import argparse
import statistics
import parakit.parakit_io as pkio
import parakit.parakit_process as pkproc
import parakit.parakit_variants as pkvar
import parakit.parakit_path as pkpath
import parakit.parakit_plot as pkplot
import parakit.parakit_surject as pksurj
import json
import os


pars = argparse.ArgumentParser()
spars = pars.add_subparsers(help='sub-command help', required=True)

# construct subcommand: construct the pangenome
pars_construct = spars.add_parser('construct',
                                  help='construct pangenome')
pars_construct.add_argument('-j', help='config JSON file', required=True)
pars_construct.add_argument('-o', help='prefix for the output files',
                            default='')
pars_construct.add_argument('-t', default=4, type=int,
                            help='number of threads to use')
pars_construct.set_defaults(scmd='construct')


def scmd_construct(args):
    # assumes a 'seqs' directory exists TODO
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # check that docker is intalled TODO
    # prepare output file names
    pg_gfa = pkio.gfaFile(config, prefix=args.o, check_file=False)
    node_tsv = pkio.nodeFile(config, prefix=args.o, check_file=False)
    opref = pkio.prefixFile(config)
    refname = 'ref'
    if os.path.isfile(pg_gfa):
        print("{} exists. Skipping pangenome construction.".format(pg_gfa))
    else:
        # prepare files and run pangenome construction
        if config['method'] == 'mc':
            # use minigraph-cactus
            cons_proc = pkproc.constructPgMc(config, opref, pg_gfa)
        elif config['method'] == 'mc_collapse':
            # use minigraph-cactus experiment collapse mode
            cons_proc = pkproc.constructPgMcCollapse(config, opref, pg_gfa)
        elif config['method'] == 'pggb':
            # use PGGB
            cons_proc = pkproc.constructPgPggb(config, opref, pg_gfa,
                                               threads=args.t)
        refname = cons_proc['refname']
    # annotate the nodes in the graph
    pkio.readGFA(pg_gfa, refname=refname, out_tsv=node_tsv,
                 guess_modules=config['method'] != 'mc')
    # done
    print('Output GFA: ' + pg_gfa)
    print('Node information: ' + node_tsv)


# deconstruct: list all the variants in the pangenome
pars_deconstruct = spars.add_parser('deconstruct',
                                    help='list variants in pangenome')
pars_deconstruct.add_argument('-j', help='config JSON file', required=True)
pars_deconstruct.add_argument('-n', help='node information', default='')
pars_deconstruct.add_argument('-g', help='input GFA pangenome', default='')
pars_deconstruct.add_argument('-a', help='annotation file (e.g. from ClinVar)',
                              default='')
pars_deconstruct.add_argument('-o', help='output TSV file', default='variants.tsv')
pars_deconstruct.add_argument('-t', help='debug trace mode',
                              action='store_true')
pars_deconstruct.set_defaults(scmd='deconstruct')


def scmd_deconstruct(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # get offset from config file
    c1, c2, pos_offset, reg_e = pkproc.getRegionsFromConfig(config)
    pos_offset += 1

    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn, verbose=args.t)

    # update with edge information
    pg_gfa = pkio.gfaFile(config, fn=args.g, check_file=True)
    pkio.updateNodesSucsWithGFA(nodes, pg_gfa, verbose=args.t)

    # read annotation (optional)
    clinvar_fn = pkio.clinvarFile(args.a, config)
    if not os.path.isfile(clinvar_fn):
        clinvar_fn = None

    # decompose variants and annotate
    pkvar.decomposePangenome(nodes, annot_fn=clinvar_fn, pos_offset=pos_offset,
                             output_fn=args.o, chrom=c1[0])


# map subcommand: map reads to the pangenome
pars_map = spars.add_parser('map',
                            help='map reads to the pangenome')
pars_map.add_argument('-j', help='config JSON file', required=True)
pars_map.add_argument('-g', help='input GFA pangenome', default='')
pars_map.add_argument('-b', help='input (indexed) BAM file', required=True)
pars_map.add_argument('-o', help='output GAF file', required=True)
pars_map.add_argument('-t', help='debug trace mode', action='store_true')
pars_map.set_defaults(scmd='map')


def scmd_map(args):
    # assumes a 'seqs' directory exists TODO
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # extract local reads to fastq
    fq_fn = args.o + '.fq'
    if (args.b.endswith('.fq') or args.b.endswith('.fastq')
            or args.b.endswith('.fq.gz') or args.b.endswith('.fastq.gz')):
        print('Using input FASTQ: {}'.format(args.b))
        fq_fn = args.b
    else:
        pkproc.extractReads(config, args.b, fq_fn, trace=args.t)
    # map reads to pangenome
    pg_gfa = pkio.gfaFile(config, fn=args.g, check_file=True)
    pkproc.mapReads(fq_fn, pg_gfa, args.o)
    print('Output GAF: ' + args.o)


# call subcommand: calls variant by aggregating read support
pars_call = spars.add_parser('call',
                             help='call variants by aggregating read support')
pars_call.add_argument('-j', help='config JSON file', required=True)
pars_call.add_argument('-n', help='node information', default='')
pars_call.add_argument('-r', help='input alignments in GAF', required=True)
pars_call.add_argument('-g', help='input GFA pangenome', default='')
pars_call.add_argument('-a', help='annotation file (e.g. from ClinVar)',
                       default='')
pars_call.add_argument('-m', help='number of markers checked around the SNPs',
                       default=20, type=int)
pars_call.add_argument('-o', help='output TSV', required=True)
pars_call.add_argument('-t', help='debug trace mode', action='store_true')
pars_call.set_defaults(scmd='call')


def scmd_call(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # get offset from config file
    c1, c2, pos_offset, reg_e = pkproc.getRegionsFromConfig(config)
    pos_offset += 1

    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn, verbose=args.t)

    # update with edge information
    pg_gfa = pkio.gfaFile(config, fn=args.g, check_file=True)
    pkio.updateNodesSucsWithGFA(nodes, pg_gfa, verbose=args.t)

    # read annotation
    clinvar_fn = pkio.clinvarFile(args.a, config, check_file=True)

    # read GAF
    reads = pkio.readGAF(args.r, nodes, verbose=args.t)

    # find variants in reads and write output TSV
    pkvar.findVariants(nodes=nodes, annot_fn=clinvar_fn,
                       reads=reads, nmarkers=args.m,
                       pos_offset=pos_offset,
                       output_tsv=args.o, chrom=c1[0])


# call subcommand: calls variant by aggregating read support
pars_filtercalls = spars.add_parser('filtercalls',
                                    help='filter called variants')
pars_filtercalls.add_argument('-j', help='config JSON file', required=True)
pars_filtercalls.add_argument('-n', help='node information', default='')
pars_filtercalls.add_argument('-c', help='calls TSV', required=True)
pars_filtercalls.add_argument('-r', default=None,
                              help='region (START-END positions)')
pars_filtercalls.add_argument('-m', default=None,
                              help='the module of interest (1 or 2)')
pars_filtercalls.add_argument('-a', action='store_true',
                              help='matched with ClinVar variants')
pars_filtercalls.add_argument('-f', action='store_true',
                              help='fusion variants')
pars_filtercalls.add_argument('-o', help='output TSV', required=True)
pars_filtercalls.set_defaults(scmd='filtercalls')


def scmd_filtercalls(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # get offset from config file
    c1, c2, pos_offset, reg_e = pkproc.getRegionsFromConfig(config)
    pos_offset += 1

    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn)

    # parse region
    start_pos = None
    end_pos = None
    if args.r is not None:
        reg = args.r.replace(',', '').split('-')
        start_pos = int(reg[0])
        end_pos = int(reg[1])

    # find variants in reads and write output TSV
    pkvar.filterVariants(nodes=nodes, calls_fn=args.c,
                         module=args.m, annotated_only=args.a,
                         fusion_only=args.f,
                         start_pos=start_pos, end_pos=end_pos,
                         output_tsv=args.o)


# diplotype subcommand: find the best pair of paths by aggregating read support
pars_diplotype = spars.add_parser('diplotype',
                                  help='find diplotype best supported '
                                  'by reads')
pars_diplotype.add_argument('-j', help='config JSON file', default='')
pars_diplotype.add_argument('-n', help='node information', default='')
pars_diplotype.add_argument('-g', help='input GFA pangenome', default='')
pars_diplotype.add_argument('-r', help='input alignments in GAF',
                            required=True)
pars_diplotype.add_argument('-o', help='output TSV prefix', required=True)
pars_diplotype.add_argument('-m', help='method to use (original, default, fast, em, kcomms). Default: default')
pars_diplotype.add_argument('-t', help='debug trace mode', action='store_true')
pars_diplotype.set_defaults(scmd='diplotype')


def scmd_diplotype(args):
    # read config json file
    config = {}
    if args.j != '':
        config = json.load(open(args.j, 'rt'))

    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn, verbose=args.t)

    # update with edge information
    pg_gfa = pkio.gfaFile(config, fn=args.g, check_file=True)
    pkio.updateNodesSucsWithGFA(nodes, pg_gfa, verbose=args.t)

    # read GAF
    reads = pkio.readGAF(args.r, nodes, verbose=args.t)

    # find paths
    paths_res = pkpath.findPaths(nodes, reads, config, args)
    pkio.writePathsInfo(paths_res, nodes,
                        stats_fn=args.o + '.paths-stats.tsv',
                        info_fn=args.o + '.paths-info.tsv')


# annotate subcommand: annotate sequence(s) from a FASTA file
pars_annotate = spars.add_parser('annotate',
                                 help='annotate sequence(s) from a FASTA file')
pars_annotate.add_argument('-j', help='config JSON file', required=True)
pars_annotate.add_argument('-n', help='node information', default='')
pars_annotate.add_argument('-g', help='input GFA pangenome', default='')
pars_annotate.add_argument('-f', help='input fasta file', default='')
pars_annotate.add_argument('-r', help='input alignments in GAF', default='')
pars_annotate.add_argument('-p', help='annotate the paths in the GFA',
                           action='store_true')
pars_annotate.add_argument('-o', help='output PDF file',
                           default='parakit.viz.pdf')
pars_annotate.add_argument('-e', help='input genome element annotation TSV',
                           default='')
pars_annotate.add_argument('-a', help='annotation file (e.g. from ClinVar)',
                           default='')
pars_annotate.add_argument('-m', default=20, type=int,
                           help='number of markers checked around the SNPs')
pars_annotate.add_argument('-s', default='pangenome',
                           help='Optional. Scale for the x-axis. '
                           'Either pangenome or genome. ')
pars_annotate.add_argument('-t', help='debug trace mode', action='store_true')
pars_annotate.set_defaults(scmd='annotate')


def scmd_annotate(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # if we map fasta sequence(s) to pangenome, use this GAF file
    gaf_fn = args.f + '.' + config['method'] + '.gaf.gz'
    # pangenome GFA
    args.g = pkio.gfaFile(config, fn=args.g, check_file=True)
    if args.p:
        print("Using annotating the pangenome's paths.")
    elif args.r != '':
        print('Using {}. No alignment needed.'.format(args.r))
    elif args.f != '':
        if os.path.isfile(gaf_fn):
            print(gaf_fn + ' exists. Skipping alignment to the pangenome')
            args.r = gaf_fn
        else:
            if not os.path.isfile(args.f):
                print('Cannot find file: ' + args.f)
                exit(1)
            print('Aligning {} to pangenome...'.format(args.f))
            pkproc.mapReads(args.f, args.g, gaf_fn)
            print('Output GAF: ' + gaf_fn)
            args.r = gaf_fn
    else:
        print('Either -f, -r, or -p must be used.')
        exit(1)
    # visualize using R script
    script_path = os.path.join(os.path.dirname(__file__),
                               'parakit.viz.R')
    # update/guess paths before
    args.n = pkio.nodeFile(config, fn=args.n, check_file=True)
    args.e = pkio.geneFile(args.e, config, check_file=True)
    # if annotation provided, call variants
    if args.a != '':
        # get offset from config file
        c1, c2, pos_offset, reg_e = pkproc.getRegionsFromConfig(config)
        pos_offset += 1
        nodes = pkio.readNodeInfo(args.n, verbose=args.t)
        # update with edge information
        pg_gfa = pkio.gfaFile(config, fn=args.g, check_file=True)
        pkio.updateNodesSucsWithGFA(nodes, pg_gfa, verbose=args.t)
        # read annotation
        clinvar_fn = pkio.clinvarFile(args.a, config, check_file=True)
        # read GAF
        reads = pkio.readGAF(args.r, nodes, verbose=args.t)
        # variants in reads and write output TSV
        pkvar.findVariants(nodes=nodes, annot_fn=clinvar_fn,
                           reads=reads, nmarkers=args.m,
                           pos_offset=pos_offset, min_support=1,
                           output_tsv=args.o + '.calls.tsv')
    # run the R script to make the graphs
    pkproc.runRscript(script_path, args)


# viz subcommand: make different graphs of the results
pars_viz = spars.add_parser('viz', help='visualize the results')
pars_viz.add_argument('-v', help='visualization mode: allele_support, calls, '
                      'diplotype, all, all_small',
                      default='all_small')
pars_viz.add_argument('-j', help='config JSON file', required=True)
pars_viz.add_argument('-n', help='node information', default='')
pars_viz.add_argument('-e', help='input genome element annotation TSV',
                      default='')
pars_viz.add_argument('-r', default='', help='input alignments in GAF')
pars_viz.add_argument('-c', help='calls TSV', default='')
pars_viz.add_argument('--dip', help='diplotype index', default=1, type=int)
pars_viz.add_argument('-d', help='diplotype paths, sorted', default='')
pars_viz.add_argument('-p', help='haplotype paths information', default='')
pars_viz.add_argument('-l', help='a label to use as title of the graphs.',
                      default='')
pars_viz.add_argument('-m', help='maximum number of supporting reads to '
                      'show in graph', default=3)
pars_viz.add_argument('-o', help='output PDF file', default='parakit.viz.pdf')
pars_viz.add_argument('-S', help='Optional. R script to run instead of '
                      'the one provided.', default='')
pars_viz.add_argument('-s', default='pangenome',
                      help='Optional. Scale for the x-axis. '
                      'Either pangenome or genome. ')
pars_viz.add_argument('-t', help='debug trace mode', action='store_true')
pars_viz.set_defaults(scmd='viz')


def scmd_viz(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # visualize using R script
    script_path = os.path.join(os.path.dirname(__file__),
                               'parakit.viz.R')
    if args.S != '':
        script_path = args.S
    # update/guess paths before
    args.n = pkio.nodeFile(config, fn=args.n, check_file=True)
    args.e = pkio.geneFile(args.e, config, check_file=True)
    # run the R script to make the graphs
    pkproc.runRscript(script_path, args)
    return (True)


# viz subcommand: make different graphs of the results
pars_plot = spars.add_parser('plot', help='plot the results')
pars_plot.add_argument('-j', help='config JSON file', required=True)
pars_plot.add_argument('-n', help='node information', default='')
pars_plot.add_argument('-g', help='input GFA pangenome', default='')
pars_plot.add_argument('-e', help='input genome element annotation TSV',
                       default='')
pars_plot.add_argument('-r', default='', help='input alignments in GAF')
pars_plot.add_argument('-d', help='diplotype paths, sorted', default='')
pars_plot.add_argument('-p', help='haplotype paths information', default='')
pars_plot.add_argument('-m', help='plot mode', default='subreads')
pars_plot.add_argument('-o', help='output PDF file', default='parakit.viz.pdf')
pars_plot.set_defaults(scmd='plot')


def scmd_plot(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # get offset from config file
    c1, c2, pos_offset, reg_e = pkproc.getRegionsFromConfig(config)
    pos_offset += 1
    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn)
    # update with edge information
    pg_gfa = pkio.gfaFile(config, fn=args.g, check_file=True)
    pkio.updateNodesSucsWithGFA(nodes, pg_gfa)
    gene_fn = pkio.geneFile(args.e, config, check_file=True)
    genes = pkio.readGeneAnnotation(gene_fn, nodes,
                                    pos_offset, config['genes'])
    # read GAF
    reads = None
    if args.r != '':
        reads = pkio.readGAF(args.r, nodes, verbose=True)
    # read diplotype
    dip_paths = None
    if args.d != '' and args.p != '':
        dip_paths = pkio.readDiplotype(args.d, args.p, verbose=True)
    # if only nodes and reads -> packed subread plot
    if reads is not None and dip_paths is None:
        if args.m == 'subreads':
            pkplot.plotPackedSubreads(nodes, reads, genes, args.o)
        elif args.m == 'longreads':
            pkplot.plotLongReads(nodes, reads, genes, args.o)
    # if nodes, reads and a diplotype -> reads on haplotypes plot
    if reads is not None and dip_paths is not None:
        pkplot.plotDiplotypeReads(nodes, reads, dip_paths, args.o)
    return (True)


# copy subcommand: estimate the copy number of the modules
pars_copy = spars.add_parser('copy',
                             help='estimate the copy number of the modules')
pars_copy.add_argument('-j', help='config JSON file', required=True)
pars_copy.add_argument('-n', help='node information', default='')
pars_copy.add_argument('-r', default='', help='input alignments in GAF')
pars_copy.add_argument('-o', help='Output TSV file prefix',
                       default='parakit.cn')
pars_copy.set_defaults(scmd='copy')


def scmd_copy(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # update/guess paths before
    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn)
    # read GAF
    reads = pkio.readGAF(args.r, nodes, verbose=True)
    pkvar.estimateCopyNumber(nodes, reads, out_tsv_prefix=args.o)
    return (True)


# gafstats subcommand: estimate the copy number of the modules
pars_gafstats = spars.add_parser('gafstats',
                                 help='compute basis stats from a GAF')
pars_gafstats.add_argument('-r', required=True, help='input alignments in GAF')
pars_gafstats.add_argument('-j', help='config JSON file', required=True)
pars_gafstats.add_argument('-n', help='node information', default='')
pars_gafstats.add_argument('-o', help='Optional output TSV file', default='')
pars_gafstats.set_defaults(scmd='gafstats')


def scmd_gafstats(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # update/guess paths before
    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn)
    # read GAF
    stats = pkio.readGAFstats(args.r, nodes)
    # print summary stats
    median_rlen = statistics.median(stats['length'])
    print('Median read length: {} bp'.format(median_rlen))
    print('Maximum read length: {} bp'.format(max(stats['length'])))
    # N50
    rlens = list(stats['length'])
    rlens.sort()
    rlen_csum = 0
    rlen_sum = sum(rlens)
    for rlen in rlens:
        rlen_csum += rlen
        if rlen_csum > rlen_sum * .5:
            print('N50 read length: {} bp'.format(rlen))
            break
    # alignment stats
    prop_aln = []
    prop_match = []
    for ii in range(len(stats['length'])):
        prop_aln.append(float(stats['aligned'][ii]) / stats['length'][ii])
        prop_match.append(float(stats['matches'][ii]) / stats['length'][ii])
    prop_aln = round(100 * statistics.median(prop_aln), 2)
    print('Alignment blocks are ~{}% of the read length '
          '(median)'.format(prop_aln))
    prop_match = round(100 * statistics.median(prop_match), 2)
    print('Median ~{}% of a read matches the pangenome'.format(prop_match))
    # write output TSV if output file provided
    if args.o != '':
        outf = open(args.o, 'wt')
        outf.write('\t'.join(list(stats.keys())) + '\n')
        for ii in range(len(stats['length'])):
            rec = []
            for metric in list(stats.keys()):
                rec.append(str(stats[metric][ii]))
            outf.write('\t'.join(rec) + '\n')
        outf.close()
    return (True)


# surject subcommand: surject alignments to the linear reference
pars_surject = spars.add_parser('surject',
                                help='project reads to the linear reference')
pars_surject.add_argument('-r', required=True, help='input alignments in GAF')
pars_surject.add_argument('-j', help='config JSON file', required=True)
pars_surject.add_argument('-g', help='input GFA pangenome', default='')
pars_surject.add_argument('-f', help='input fastq reads', default='')
pars_surject.add_argument('-n', help='node information', default='')
pars_surject.add_argument('-d', help='diplotype paths, sorted', default='')
pars_surject.add_argument('-p', help='haplotype paths information', default='')
pars_surject.add_argument('-o', help='output BAM file', default='')
pars_surject.add_argument('-t', help='debug trace mode', action='store_true')
pars_surject.set_defaults(scmd='surject')


def scmd_surject(args):
    # default output name
    if args.o == '':
        args.o = 'surjected_subreads.bam'
    # try to guess the fastq file?
    if args.f == '':
        fq_fn = args.r + '.fq'
        if os.path.isfile(fq_fn):
            print("Guessing that the fastq is " + fq_fn)
            args.f = fq_fn
        else:
            print("Couldn't guess the fastq file. Please provide it with -f.")
            exit(1)
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # update/guess paths before
    # load node info
    node_fn = pkio.nodeFile(config, fn=args.n, check_file=True)
    nodes = pkio.readNodeInfo(node_fn)
    # update with edge information
    pg_gfa = pkio.gfaFile(config, fn=args.g, check_file=True)
    pkio.updateNodesSucsWithGFA(nodes, pg_gfa, verbose=args.t)
    # read GAF
    reads = pkio.readGAF(args.r, nodes, verbose=args.t)
    # read diplotype
    dip_paths = None
    if args.d != '' and args.p != '':
        dip_paths = pkio.readDiplotype(args.d, args.p, verbose=args.t)
    # identify subreads of interest
    if dip_paths is None:
        sel_sreads = pksurj.selectSubreads(reads, nodes, module='c2')
    else:
        sel_sreads = pksurj.selectSubreadsOnDip(reads, nodes, dip_paths)
    # extract them as fastqs
    sreads_fq = args.o + '.fq'
    pkproc.extractSubreads(sel_sreads, args.f, sreads_fq)
    # map them to the linear reference
    sreads_sam = pkproc.mapMinimap2Local(sreads_fq, config['ref_fa'],
                                         config['c2'])
    os.remove(sreads_fq)
    # annotate with read name, and potentially haplotype/module info
    pkproc.annotateSam(sreads_sam, sel_sreads, args.o, verbose=args.t)
    return (True)


def main():
    args = pars.parse_args()

    if args.scmd == 'construct':
        scmd_construct(args)
    elif args.scmd == 'deconstruct':
        scmd_deconstruct(args)
    elif args.scmd == 'map':
        scmd_map(args)
    elif args.scmd == 'call':
        scmd_call(args)
    elif args.scmd == 'filtercalls':
        scmd_filtercalls(args)
    elif args.scmd == 'diplotype':
        scmd_diplotype(args)
    elif args.scmd == 'viz':
        scmd_viz(args)
    elif args.scmd == 'plot':
        scmd_plot(args)
    elif args.scmd == 'copy':
        scmd_copy(args)
    elif args.scmd == 'gafstats':
        scmd_gafstats(args)
    elif args.scmd == 'annotate':
        scmd_annotate(args)
    elif args.scmd == 'surject':
        scmd_surject(args)
