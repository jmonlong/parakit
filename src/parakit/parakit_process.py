import parakit.parakit_io as pkio
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
import pyfaidx


def getRegionsFromConfig(config):
    c1 = config['c1'].split(':')
    c1[1] = [int(xx) for xx in c1[1].split('-')]
    c2 = config['c2'].split(':')
    c2[1] = [int(xx) for xx in c2[1].split('-')]
    # make sure c1 is upstream of c2 in the reference genome
    if c1[1][0] > c2[1][0]:
        ctemp = c2
        c2 = c1
        c1 = ctemp
    # region to extract
    reg_s = c1[1][0] - config['flank_size']
    reg_e = c2[1][1] + config['flank_size']
    return (c1, c2, reg_s, reg_e)


def prepareRefSeqsForMc(config):
    print('Extracting reference sequences.')
    ref_fa = pyfaidx.Fasta(config['ref_fa'])
    c1, c2, reg_s, reg_e = getRegionsFromConfig(config)
    # write files for the SD-ome graph
    recs = []
    seq = str(ref_fa[c1[0]][reg_s:reg_e])
    recs.append(SeqRecord(MutableSeq(seq.upper()), id='ref', description=""))
    SeqIO.write(recs, "seqs/ref_full.fa", "fasta")
    recs = []
    seq = str(ref_fa[c1[0]][reg_s:c2[1][0]]) + \
        str(ref_fa[c1[0]][c2[1][1]:reg_e])
    recs.append(SeqRecord(MutableSeq(seq.upper()),
                          id='ref_noc2', description=""))
    SeqIO.write(recs, "seqs/ref_noc2.fa", "fasta")
    recs = []
    seq = str(ref_fa[c1[0]][c1[1][0]:c2[1][1]])
    recs.append(SeqRecord(MutableSeq(seq.upper()),
                          id='ref_c1c2', description=""))
    SeqIO.write(recs, "seqs/ref_c1c2.fa", "fasta")
    recs = []
    seq = str(ref_fa[c1[0]][c1[1][0]:c1[1][1]])
    recs.append(SeqRecord(MutableSeq(seq.upper()),
                          id='c1_ref', description=""))
    SeqIO.write(recs, "seqs/c1_ref.fa", "fasta")
    recs = []
    seq = str(ref_fa[c1[0]][c2[1][0]:c2[1][1]])
    recs.append(SeqRecord(MutableSeq(seq.upper()),
                          id='c2_ref', description=""))
    SeqIO.write(recs, "seqs/c2_ref.fa", "fasta")


def prepareHprcSeqsForMc(config):
    hprc_seqs_for_mc = []
    coord_inf = open(config['hprc_coords'], 'rt')
    for line in coord_inf:
        [samp, coord, label] = line.rstrip().split('\t')
        out_fa = 'seqs/' + label + '.fa'
        if os.path.isfile(out_fa):
            print('{} already exists, skipping...'.format(out_fa))
        else:
            print('Extracting HPRC sequence {} at {}.'.format(label, coord))
            agc_cmd = ['agc',  'getctg', config['hprc_agc'], coord]
            agc_o = subprocess.run(agc_cmd, check=True, capture_output=True)
            out_file = open(out_fa, 'wt')
            for ii, line in enumerate(agc_o.stdout.decode().split('\n')):
                if ii == 0:
                    out_file.write('>' + label + '\n')
                else:
                    out_file.write(line + '\n')
            out_file.close()
        hprc_seqs_for_mc.append([label, out_fa])
    coord_inf.close()
    return (hprc_seqs_for_mc)


def constructPgMc(config, opref, pg_gfa):
    # prepare reference sequence(s)
    prepareRefSeqsForMc(config)
    # prepare mc input info
    mc_info_fn = opref + '.for_mc.txt'
    mc_info_f = open(mc_info_fn, 'wt')
    mc_info_f.write('ref_noc2 seqs/ref_noc2.fa\n')
    mc_info_f.write('ref_c1c2 seqs/ref_c1c2.fa\n')
    mc_info_f.write('c1_ref seqs/c1_ref.fa\n')
    mc_info_f.write('c2_ref seqs/c2_ref.fa\n')
    # prepare local sequences from HPRC
    if 'hprc_agc' in config and 'hprc_coords' in config:
        hprc_seqs_for_mc = prepareHprcSeqsForMc(config)
        for hseq in hprc_seqs_for_mc:
            mc_info_f.write('{} {}\n'.format(hseq[0], hseq[1]))
    mc_info_f.close()
    # prepare script with cactus command
    mc_sh_fn = opref + '.for_mc.sh'
    mc_js_fn = opref + '.for_mc.js'
    mc_outdir_fn = opref + '.pg'
    mc_sh_f = open(mc_sh_fn, 'wt')
    mc_sh_f.write('cactus-pangenome ' + mc_js_fn + ' ' + mc_info_fn +
                  ' --outDir ' + mc_outdir_fn + ' --outName mc_pg' +
                  ' --reference ref_noc2 --gfa\n')
    mc_sh_f.close()
    # get USER id to make sure the file permission are correct with docker
    id_o = subprocess.run(['id', '-u', os.getenv('USER')],
                          check=True, capture_output=True)
    mc_cmd = ['docker', 'run', '-it', '-v', os.getcwd() + ':/app',
              '-w', '/app',
              '-u', id_o.stdout.decode().rstrip(),
              'quay.io/comparative-genomics-toolkit/cactus:v2.6.7',
              'sh', mc_sh_fn]
    if os.path.isfile(mc_outdir_fn + '/mc_pg.gfa.gz') or \
       os.path.isfile(mc_outdir_fn + '/mc_pg.gfa'):
        print("Skipping Cactus-Minigraph.")
    else:
        subprocess.run(mc_cmd, check=True)
    # add reference path by aligning with GraphAligner
    if not os.path.isfile(mc_outdir_fn + '/mc_pg.gfa'):
        subprocess.run(['gunzip', mc_outdir_fn + '/mc_pg.gfa.gz'],
                       check=True)
    ga_cmd = ['docker', 'run', '-it', '-v', os.getcwd() + ':/app',
              '-w', '/app',
              '-u', id_o.stdout.decode().rstrip(),
              'quay.io/biocontainers/graphaligner:1.0.17b--h21ec9f0_2',
              'GraphAligner', '-g', mc_outdir_fn + '/mc_pg.gfa',
              '-f', 'seqs/ref_full.fa', '-a', mc_outdir_fn + '/ref.gaf',
              '-x', 'vg']
    outf = open(mc_outdir_fn + '/mc_pg.pg', 'w')
    subprocess.run(ga_cmd, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(['vg', 'convert', '-p', mc_outdir_fn + '/mc_pg.gfa'],
                   check=True, stdout=outf)
    outf.close()
    outf = open(mc_outdir_fn + '/mc_pg_with_ref_raw.pg', 'w')
    subprocess.run(['vg', 'augment', '-F',
                    '-B', mc_outdir_fn + '/mc_pg.pg',
                    mc_outdir_fn + '/ref.gaf'], check=True, stdout=outf)
    outf.close()
    outf = open(mc_outdir_fn + '/mc_pg_with_ref.pg', 'w')
    subprocess.run(['vg', 'mod', '-O',
                    mc_outdir_fn + '/mc_pg_with_ref_raw.pg'],
                   check=True, stdout=outf)
    outf.close()
    outf = open(pg_gfa, 'w')
    subprocess.run(['vg', 'view', mc_outdir_fn + '/mc_pg_with_ref.pg'],
                   check=True, stdout=outf)
    outf.close()
    # remove temporary files
    for ff in [mc_sh_fn, mc_info_fn]:
        os.remove(ff)
    return ({'refname': 'ref'})


def constructPgPggb(config, opref, pg_gfa):
    # use PGGB    
    # docker run -it -v `pwd`:/app -w /app -u `id -u $USER` ghcr.io/pangenome/pggb:latest
    # cat seqs/ref.fa seqs/ref.b.fa > in.fa
    # samtools faidx in.fa
    # pggb -i in.fa -o pggb.out -n 2 -c 2 -t 4
    cmd = ['docker run -it -v `pwd`:/app -w /app -u `id -u $USER` ghcr.io/pangenome/pggb:latest']
    subprocess.run(cmd, check=True)
    # move/copy final GFA to pg_gfa
    # refname might be different when using pggb?
    return ({'refname': ''})


def extractReads(config, in_reads, out_fq):
    if os.path.isfile(out_fq):
        print("{} already exists. Using those"
              " extracted reads".format(out_fq))
        return ()
    c1, c2, reg_s, reg_e = getRegionsFromConfig(config)
    reg = '{}:{}-{}'.format(c1[0], reg_s, reg_e)
    # run samtools
    temp_bam = out_fq + '.bam'
    subprocess.run(['samtools', 'view', '-h', '-O', 'bam',
                    '-o', temp_bam, in_reads, reg], check=True)
    subprocess.run(['samtools', 'fastq', '-o', out_fq, '-0', out_fq, temp_bam],
                   check=True)
    os.remove(temp_bam)


def mapReads(in_fq, pg_gfa, out_gaf):
    # get USER id to make sure the file permission are correct with docker
    id_o = subprocess.run(['id', '-u', os.getenv('USER')],
                          check=True, capture_output=True)
    # GraphAligner outputs an unzipped GAF file
    ga_gaf = out_gaf
    if out_gaf.endswith('.gz'):
        # if we want an gzipped output, make a temporary GAF
        # that will be zipped later
        ga_gaf = out_gaf + '.gaf'
    ga_cmd = ['docker', 'run', '-it', '-v', os.getcwd() + ':/app',
              '-w', '/app',
              '-u', id_o.stdout.decode().rstrip(),
              'quay.io/biocontainers/graphaligner:1.0.17b--h21ec9f0_2',
              'GraphAligner', '-g', pg_gfa, '-f', in_fq, '-a', ga_gaf,
              '-x', 'vg']
    subprocess.run(ga_cmd, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if out_gaf.endswith('.gz'):
        # if we wanted an gzipped output, zip the temporary GAF
        out_gaf_f = open(out_gaf, 'w')
        subprocess.run(['gzip', '-c', ga_gaf], check=True, stdout=out_gaf_f)
        out_gaf_f.close()
        # and delete it
        os.remove(ga_gaf)


def runRscript(config, script_r, args, annotate_gaf=''):
    # get offset from config file
    c1, c2, pos_offset, reg_e = getRegionsFromConfig(config)
    pos_offset += 1

    # visualization mode
    v_mode = 'annotate'
    if 'v' in args:
        v_mode = args.v

    # output pdf
    out_pdf = args.o
    if not args.o.endswith('.pdf'):
        out_pdf = os.path.join(out_pdf, '.pdf')

    # run R script
    rscript_cmd = ['Rscript', script_r, '-s', str(pos_offset), '-v', v_mode,
                   '-n', args.n, '-e', args.e, '-o', out_pdf]

    # add arguments if provided
    if 'm' in args:
        if v_mode == 'all_small':
            args.m = 1
        rscript_cmd += ['-m', str(args.m)]

    # check input arguments
    if v_mode in ['all', 'calls', 'all_small']:
        if args.c == '':
            print('This visualization mode requires calls input (-c).')
            exit(1)
        rscript_cmd += ['-c', args.c]
    if v_mode in ['all', 'all_small', 'calls', 'allele_support']:
        if args.r == '':
            print('This visualization mode requires reads input (-r).')
            exit(1)
        # load node info
        nodes = pkio.readNodeInfo(args.n)
        # read GAF
        reads = pkio.readGAF(args.r, nodes)
        # write read info
        reads_fn = out_pdf + '.tsv'
        outf = open(reads_fn, 'wt')
        outf.write('read\tnode\tstartpos\tendpos\treadpos\n')
        ofmt = '{}\t{}\t{}\t{}\t{}\n'
        for readn in reads.path:
            any_cspec = False
            nodes_toprint = []
            for nii, nod in enumerate(reads.path[readn]):
                if nodes[nod]['class'] != 'none':
                    any_cspec = True
                nodes_toprint.append(ofmt.format(readn,
                                                 nod,
                                                 reads.startpos[readn][nii],
                                                 reads.endpos[readn][nii],
                                                 reads.readpos[readn][nii]))
            if any_cspec:
                for tout in nodes_toprint:
                    outf.write(tout)
        outf.close()
        rscript_cmd += ['-r', reads_fn]
    if v_mode in ['all', 'all_small', 'paths']:
        if args.d == '' or args.p == '':
            print('This visualization mode requires paths input (-d and -p).')
            exit(1)
        rscript_cmd += ['-d', args.d, '-p', args.p]
    if annotate_gaf != '':
        rscript_cmd += ['-p', annotate_gaf]
    if args.l != '':
        rscript_cmd += ['-l', args.l]

    if args.t:
        subprocess.run(rscript_cmd, check=True)
    else:
        subprocess.run(rscript_cmd, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)