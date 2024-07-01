import argparse
import parakit.parakit_io as pkio
import parakit.parakit_process as pkproc
import parakit.parakit_variants as pkvar
import parakit.parakit_path as pkpath
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
pars_construct.set_defaults(scmd='construct')

# map subcommand: map reads to the pangenome
pars_map = spars.add_parser('map',
                            help='map reads to the pangenome')
pars_map.add_argument('-j', help='config JSON file', required=True)
pars_map.add_argument('-g', help='input GFA pangenome', required=True)
pars_map.add_argument('-b', help='input (indexed) BAM file', required=True)
pars_map.add_argument('-o', help='output GAF file', required=True)
pars_map.set_defaults(scmd='map')

# call subcommand: variant calls by aggregating read support
pars_call = spars.add_parser('call',
                             help='call variants by aggregating read support')
pars_call.add_argument('-j', help='config JSON file', required=True)
pars_call.add_argument('-n', help='node information', required=True)
pars_call.add_argument('-r', help='input alignments in GAF', required=True)
pars_call.add_argument('-a', help='annotation file (e.g. from ClinVar)',
                       required=True)
pars_call.add_argument('-m', help='number of markers checked around the SNPs',
                       default=10, type=int)
pars_call.add_argument('-o', help='output TSV', required=True)
pars_call.set_defaults(scmd='call')

# paths subcommand: variant pathss by aggregating read support
pars_paths = spars.add_parser('paths',
                              help='find paths best supported by reads')
pars_paths.add_argument('-n', help='node information', required=True)
pars_paths.add_argument('-g', help='input GFA pangenome', required=True)
pars_paths.add_argument('-r', help='input alignments in GAF', required=True)
pars_paths.add_argument('-o', help='output TSV prefix', required=True)
pars_paths.add_argument('-c', default=4, type=float,
                        help='minimum read support for subread clustering.')
pars_paths.set_defaults(scmd='paths')

# paths subcommand: variant pathss by aggregating read support
pars_viz = spars.add_parser('viz', help='visualize the results')
pars_viz.add_argument('-v', help='visualization mode: allele_support',
                      default='allele_support')
pars_viz.add_argument('-j', help='config JSON file', required=True)
pars_viz.add_argument('-n', help='node information', required=True)
pars_viz.add_argument('-e', help='input genome element annotation TSV',
                      required=True)
pars_viz.add_argument('-r', default='', help='input alignments in GAF')
pars_viz.add_argument('-c', help='calls TSV', default='')
pars_viz.add_argument('-d', help='diplotype paths, sorted', default='')
pars_viz.add_argument('-p', help='haplotype paths information', default='')
pars_viz.add_argument('-l', help='a label to use as title of the graphs.',
                      default='')
pars_viz.add_argument('-m', help='maximum number of supporting reads to '
                      'show in graph', default=3)
pars_viz.add_argument('-o', help='output PDF file',
                      default='parakit.viz.pdf')
pars_viz.add_argument('-s', help='Optional. R script to run instead of '
                      'the one provided.', default='')
pars_viz.add_argument('-t', help='debug trace mode', action='store_true')
pars_viz.set_defaults(scmd='viz')


def scmd_construct(args):
    # assumes a 'seqs' directory exists TODO
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # check that docker is intalled TODO
    # make up a default ouput prefix if needed
    opref = args.o
    if opref == '':
        opref = 'parakit.{}'.format(config['method'])
    # prepare names of output files
    pg_gfa = opref + '.pg.gfa'
    node_tsv = opref + '.node_info.tsv'
    # prepare files and run pangenome construction
    if config['method'] == 'mc':
        # use minigraph-cactus
        cons_proc = pkproc.constructPgMc(config, opref, pg_gfa)
    elif config['method'] == 'pggb':
        # use PGGB
        cons_proc = pkproc.constructPgPggb(config, opref, pg_gfa)
    refname = cons_proc['refname']
    # annotate the nodes in the graph
    pkio.readGFA(pg_gfa, refname=refname, out_tsv=node_tsv)
    # done
    print('Output GFA: ' + pg_gfa)
    print('Node  information: ' + node_tsv)


def scmd_map(args):
    # assumes a 'seqs' directory exists TODO
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # extract local reads to fastq
    fq_fn = args.o + '.fq'
    pkproc.extractReads(config, args.b, fq_fn)
    # map reads to pangenome
    pkproc.mapReads(fq_fn, args.g, args.o)
    print('Output GAF: ' + args.o)


def scmd_call(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    # get offset from config file
    c1, c2, pos_offset, reg_e = pkproc.getRegionsFromConfig(config)
    pos_offset += 1

    # load node info
    nodes = pkio.readNodeInfo(args.n)

    # read annotation
    anno = pkio.readVariantAnnotation(args.a, nodes, pos_offset)
    print("   Matched " + str(len(anno['pos'])) + " input variants: ")
    for varid in anno['pos']:
        print('      n' + str(anno['pos'][varid]) + '\t' + varid)

    # read GAF
    reads = pkio.readGAF(args.r, nodes)

    var_calls = pkvar.findVariants(nodes, anno, reads,
                                   nmarkers=args.m,
                                   pos_offset=pos_offset)
    var_reads_f = var_calls['var_reads_f']
    var_reads_ref = var_calls['var_reads_ref']
    fus_reads_f = var_calls['fus_reads_f']
    read_sum = var_calls['read_sum']
    cand_func_reads = var_calls['cand_func_reads']

    # print read support summary
    print("Covered variants and their supporting reads:\n")
    # order variants by position
    var_names = sorted(list(var_reads_f), key=lambda k: anno['pos'][k])
    if len(fus_reads_f) > 0:
        var_names += list(fus_reads_f)
    for_tsv = ['\t'.join(['read', 'variant', 'node', 'allele'])]
    for read in list(read_sum.keys()) + list(cand_func_reads.keys()):
        to_print = read + '\t'
        for varid in var_names:
            for_tsv_r = '\t'.join([read, varid, anno['node'][varid]])
            if read not in read_sum or varid not in read_sum[read]:
                # check if the read support a reference node for this variant
                if varid in var_reads_ref and read in var_reads_ref[varid]:
                    # if so, print ----
                    to_print += '-' * len(varid) + '\t'
                    for_tsv.append(for_tsv_r + '\tref')
                else:
                    # if not, print a gap
                    to_print += ' ' * len(varid) + '\t'
                    for_tsv.append(for_tsv_r + '\tna')
            else:
                # the read goes through the variant, print the variant name
                to_print += varid + '\t'
                for_tsv.append(for_tsv_r + '\talt')
        print(to_print)
    print()

    # write TSV output
    print('Writing summary in TSV ' + args.o + '...')
    with open(args.o, 'wt') as outf:
        outf.write('\n'.join(for_tsv) + '\n')


def scmd_paths(args):
    # load node info
    nodes = pkio.readNodeInfo(args.n)
    pkio.updateNodesSucsWithGFA(nodes, args.g)

    # read GAF
    reads = pkio.readGAF(args.r, nodes)

    paths_res = pkpath.findPaths(nodes, reads, args)
    escores_r = paths_res['escores']
    paths = paths_res['paths']

    # write ranked list
    outf_rk = open(args.o + '.paths-stats.tsv', 'wt')
    heads_rk = ['hap1', 'hap2', 'cov_cor', 'aln_score',
                'cov_cor_adj', 'aln_score_adj', 'aln_long_prop']
    ofmt_rk = '\t'.join(['{}'] * len(heads_rk)) + '\n'
    outf_rk.write('\t'.join(heads_rk) + '\n')
    for esc in escores_r:
        outf_rk.write(ofmt_rk.format(esc['hap1'],
                                     esc['hap2'],
                                     esc['cov_cor'],
                                     esc['aln_score'],
                                     esc['cov_cor_adj'],
                                     esc['aln_score_adj'],
                                     esc['aln_long_prop']))

    # write paths
    print('Write path information...')
    outf = open(args.o + '.paths-info.tsv', 'wt')
    heads = ['hap', 'node', 'ppos', 'class', 'rpos_min', 'rpos_max']
    ofmt = '\t'.join(['{}'] * len(heads)) + '\n'
    outf.write('\t'.join(heads) + '\n')
    for mode in paths:
        for ppos, nod in enumerate(paths[mode]):
            outf.write(ofmt.format(mode,
                                   nod,
                                   ppos,
                                   nodes[nod]['class'],
                                   nodes[nod]['rpos_min'],
                                   nodes[nod]['rpos_max']))
    outf.close()


def scmd_viz(args):
    # read config json file
    config = json.load(open(args.j, 'rt'))
    script_path = os.path.join(os.path.dirname(__file__),
                               'parakit.viz.R')
    if args.s != '':
        script_path = args.s
    pkproc.runRscript(config, script_path, args)
    return (True)


def main():
    args = pars.parse_args()

    if args.scmd == 'construct':
        scmd_construct(args)
    elif args.scmd == 'map':
        scmd_map(args)
    elif args.scmd == 'call':
        scmd_call(args)
    elif args.scmd == 'paths':
        scmd_paths(args)
    elif args.scmd == 'viz':
        scmd_viz(args)
