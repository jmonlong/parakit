import gzip
import statistics
import parakit.parakit_class as pkc
import os


def readNodeInfo(filen, verbose=False):
    if verbose:
        print('Reading ' + filen + '...')
    nodes = {}
    node_sum = {}
    cyc_nodes = ['', '']
    with open(filen, 'rt') as inf:
        heads = next(inf).rstrip().split('\t')
        for line in inf:
            line = line.rstrip().split('\t')
            noden = line[heads.index('node')]
            nodes[noden] = {}
            nodes[noden]['size'] = int(line[heads.index('size')])
            nodes[noden]['c1'] = int(line[heads.index('c1')])
            nodes[noden]['c2'] = int(line[heads.index('c2')])
            nodes[noden]['ref'] = int(line[heads.index('ref')])
            nodes[noden]['rpos_min'] = int(line[heads.index('rpos_min')])
            nodes[noden]['rpos_max'] = int(line[heads.index('rpos_max')])
            node_class = line[heads.index('class')]
            nodes[noden]['class'] = node_class
            nodes[noden]['seq'] = line[heads.index('seq')]
            # tally the number of nodes in each class (for verbose mode)
            if node_class not in node_sum:
                node_sum[node_class] = 0
            node_sum[node_class] += 1
            # and cycling nodes
            if node_class == 'cyc_l':
                cyc_nodes[1] = noden
            if node_class == 'cyc_r':
                cyc_nodes[0] = noden
    if verbose:
        vout = ['cycle: ' + '-'.join(cyc_nodes)]
        for nc in node_sum:
            vout.append('{}: {}'.format(nc, node_sum[nc]))
        print('\t' + ', '.join(vout))
    return nodes


def parseCg(cg):
    nums = [str(ii) for ii in range(10)]
    num_c = ''
    res = []
    for cc in cg:
        if cc in nums:
            num_c += cc
        else:
            res.append([int(num_c), cc])
            num_c = ''
    return res


def parsePath(path):
    cur_node = ''
    cur_strand = ''
    res = []
    for ii in range(len(path)):
        if path[ii] in ['>', '<']:
            if cur_node != '':
                res.append([cur_node, cur_strand])
                cur_node = ''
            cur_strand = path[ii]
        else:
            cur_node += path[ii]
    if cur_node != '':
        res.append([cur_node, cur_strand])
    return (res)


# 'nodes' is used to know the node size
def readGAF(filen, nodes, verbose=False):
    if verbose:
        print('Reading ' + filen + '...')
    # list cycling nodes
    cyc_nodes = []
    for nod in nodes:
        if nodes[nod]['class'] == 'cyc_r':
            cyc_nodes.append(nod)
    # guess if input file is gzipped
    if filen.endswith('.gz'):
        inf_gaf = gzip.open(filen, 'rt')
    else:
        inf_gaf = open(filen, 'rt')
    # rpaths = {}
    # rpos = {}
    supp_aln_cpt = {}
    reads_tr = pkc.Reads()
    nflipped = 0
    nsupp_aln = 0
    nunmp = 0
    for line in inf_gaf:
        line = line.rstrip().split('\t')
        # if read was seen before, skip supp aln
        readn = line[0]
        if reads_tr.hasRead(readn):
            supp_aln_cpt[readn] += 1
            readn = '{}_sa{}'.format(readn, supp_aln_cpt[readn])
            nsupp_aln += 1
            continue
        else:
            supp_aln_cpt[readn] = 0
        path = parsePath(line[5])
        # find CIGAR string
        cg = ''
        for tag in line[12:]:
            tag = tag.split(':')
            if tag[0] == 'cg':
                cg = tag[2]
                break
        if cg == '':
            nunmp += 1
            continue
        # parse CIGAR string
        n_idx = 0
        node_pos = int(line[7])
        read_pos = int(line[2])
        cg = parseCg(cg)
        path_cov = []
        path_node_startpos = []
        path_node_endpos = []
        path_read_pos = []
        cur_endpos = node_pos
        for cc in cg:
            for bp in range(cc[0]):
                if cc[1] == '=':
                    # node is actually covered
                    if len(path_cov) == 0 or path_cov[-1][0] != path[n_idx][0]:
                        if len(path_cov) != 0:
                            if path_cov[-1][1] == '>':
                                path_node_endpos.append(cur_endpos)
                            else:
                                nsize = nodes[path_cov[-1][0]]['size']
                                path_node_endpos.append(nsize - cur_endpos)
                        path_cov.append(path[n_idx])
                        if path_cov[-1][1] == '>':
                            path_node_startpos.append(node_pos)
                        else:
                            nsize = nodes[path_cov[-1][0]]['size']
                            path_node_startpos.append(nsize - node_pos)
                        path_read_pos.append(read_pos)
                        cur_endpos = node_pos
                if cc[1] in ['I', 'X', '=']:
                    read_pos += 1
                if cc[1] in ['D', 'X', '=']:
                    if cc[1] in ['=']:
                        cur_endpos = node_pos
                    node_pos += 1
                    if node_pos >= nodes[path[n_idx][0]]['size']:
                        n_idx += 1
                        node_pos = 0
        if path_cov[-1][1] == '>':
            path_node_endpos.append(cur_endpos)
        else:
            nsize = nodes[path_cov[-1][0]]['size']
            path_node_endpos.append(nsize - cur_endpos)
        # flip read if mostly traversing in reverse
        if line[5].count('>') < line[5].count('<'):
            nflipped += 1
            path_cov.reverse()
            path_read_pos.reverse()
            path_node_startpos.reverse()
            path_node_endpos.reverse()
            for ii in range(len(path_node_startpos)):
                nsize = nodes[path_cov[ii][0]]['size']
                path_node_startpos[ii] = nsize - path_node_startpos[ii]
                path_node_endpos[ii] = nsize - path_node_endpos[ii]
        path_cov_node = [pc[0] for pc in path_cov]
        reads_tr.addRead(readn, path_cov_node, cyc_nodes=cyc_nodes,
                         startpos=path_node_startpos,
                         endpos=path_node_endpos,
                         readpos=path_read_pos)
    inf_gaf.close()
    if verbose:
        print('\t{} reads parsed.'.format(reads_tr.nReads()))
        print('\t{} reads were flipped to traverse the '
              'pangenome (mostly) forward.'.format(nflipped))
        print('\t{} supplementary alignment were filtered.'.format(nsupp_aln))
        print('\t{} unmapped reads.'.format(nunmp))
    return (reads_tr)


def readGAFstats(filen, nodes):
    # find flanking reference nodes (to filter reads only touching those)
    flank_nodes = []
    for nod in nodes:
        if nodes[nod]['class'] == 'ref':
            flank_nodes.append(nod)
    # prepare result object
    res = {'read': [], 'length': [], 'aligned': [], 'matches': [],
           'dv': [], 'id': []}
    # guess if input file is gzipped
    if filen.endswith('.gz'):
        inf_gaf = gzip.open(filen, 'rt')
    else:
        inf_gaf = open(filen, 'rt')
    for line in inf_gaf:
        line = line.rstrip().split('\t')
        # parse path and filter if only touching flanks
        path = parsePath(line[5])
        keep_read = False
        for node_orient in path:
            if node_orient[0] not in flank_nodes:
                # at least one touched node is not a flanking region
                keep_read = True
                break
        if not keep_read:
            continue
        # if read was seen before, skip supp aln
        res['read'].append(line[0])
        res['length'].append(int(line[1]))
        # alignment info
        res['aligned'].append(int(line[10]))
        res['matches'].append(int(line[9]))
        # find dv and/or id fields
        f_dv = 'NA'
        f_id = 'NA'
        for tag in line[12:]:
            tag = tag.split(':')
            if tag[0] == 'dv':
                f_dv = tag[2]
            elif tag[0] == 'id':
                f_id = tag[2]
        res['id'].append(f_id)
        res['dv'].append(f_dv)
    return (res)


# 'nodes' is used to know the node size
def readGFAasReads(filen, nodes, verbose=False):
    if verbose:
        print('Reading ' + filen + '...')
    inf_gfa = open(filen, 'rt')
    reads_tr = pkc.Reads()
    for line in inf_gfa:
        line = line.rstrip().split('\t')
        if line[0] == 'P':
            path = line[2].replace('+', '').replace('-', '').split(',')
            # flip path if mostly traversing in reverse
            if line[2].count('+') < line[2].count('-'):
                path.reverse()
        if line[0] == 'W':
            path = line[6].replace('<', '>').split('>')[1:]
            if line[6].count('>') < line[6].count('<'):
                path.reverse()
        if line[0] in ['W', 'P']:
            startpos = []
            endpos = []
            readpos = []
            cur_readpos = 0
            for node in path:
                readpos.append(cur_readpos)
                startpos.append(0)
                endpos.append(nodes[node]['size'])
                cur_readpos += nodes[node]['size']
            reads_tr.addRead(line[1], path, startpos=startpos,
                             endpos=endpos, readpos=readpos)
    inf_gfa.close()
    return (reads_tr)


def updateNodesSucsWithGFA(nodes, filen, verbose=False):
    if verbose:
        print('Reading ' + filen + '...')
    # read GFA and update edges
    inf = open(filen, 'rt')
    for line in inf:
        line = line.rstrip().split('\t')
        if line[0] == 'L':
            if line[1] in nodes:
                if 'sucs' not in nodes[line[1]]:
                    nodes[line[1]]['sucs'] = {}
                nodes[line[1]]['sucs'][line[3]] = True
    # init sucs for all nodes
    for nod in nodes:
        if 'sucs' not in nodes[nod]:
            if verbose:
                print('\tAdding new successor to node ', nod)
            nodes[nod]['sucs'] = {}
    inf.close()


def readGFA(gfa_fn, min_fc=3, refname='grch38', out_tsv='',
            guess_modules=False):
    # node info
    # for each node ID:
    #   size
    #   ref/c1/c2: number of times a ref/c1/c2 path traverses
    #   rpos_min/rpos_max: minimum/maximum position on the reference path
    ninfo = {}

    # read GFA and parse info
    # save paths
    inf = open(gfa_fn, 'rt')
    paths = {}
    for line in inf:
        line = line.rstrip().split('\t')
        if line[0] == 'P':
            path = line[2].replace('+', '').replace('-', '').split(',')
            # flip path if mostly traversing in reverse
            if line[2].count('+') < line[2].count('-'):
                path.reverse()
            paths[line[1]] = path
        if line[0] == 'W':
            path = line[6].replace('<', '>').split('>')[1:]
            if line[6].count('>') < line[6].count('<'):
                path.reverse()
            paths[line[3]] = path
        # save some node information like size
        # also init, future fields
        if line[0] == 'S':
            ninfo[line[1]] = {'size': len(line[2]),
                              'ref': 0, 'c1': 0, 'c2': 0,
                              'rnode': -1,
                              'rpos_min': -1, 'rpos_max': -1,
                              'class': 'none'}
            if len(line[2]) < 10:
                ninfo[line[1]]['seq'] = line[2]
            else:
                ninfo[line[1]]['seq'] = line[2][:10] + '+'
    inf.close()

    # reference path
    ref = paths[refname]

    # update node info with path cover
    refpos = 0
    for nod in ref:
        ninfo[nod]['ref'] += 1
        # also minimum/maximum position on reference
        if ninfo[nod]['rpos_min'] == -1:
            ninfo[nod]['rpos_min'] = refpos
        else:
            ninfo[nod]['rpos_min'] = min(refpos, ninfo[nod]['rpos_min'])
        ninfo[nod]['rpos_max'] = max(refpos, ninfo[nod]['rpos_max'])
        refpos += ninfo[nod]['size']

    # flag nodes contributing to the cycle edge
    # (end of copy 2 to beginning of copy 1)
    #
    #        _______________________________________
    #       |                                       |
    #  - r1 - c1 - c2 - c3 - c4 - c5 - b1 - b2 - b3 - r2 -
    #                                |              |
    #                                ----------------
    #
    # Positions on the reference path:
    #     1    2    3    4    5    6    7    8    9
    #         10   11   12   13   14                  15
    # "cycling" edge: cyc_r c5, cyc_l c1
    # find the size of the largest cycle (followed by the reference path)
    #   presumably that's the cycle we're interested in, not the smaller
    #   cycles that might exist
    cur_max_jump = 0
    for ii, nod in enumerate(ref):
        jump = ninfo[nod]['rpos_max'] - ninfo[nod]['rpos_min']
        if jump > cur_max_jump:
            cur_max_jump = jump
    # find the first/last nodes with that jump
    cyc_end = ''
    for nod in ref:
        jump = ninfo[nod]['rpos_max'] - ninfo[nod]['rpos_min']
        if jump > .9 * cur_max_jump and cyc_end == '':
            # first node with that jump
            cyc_end = nod
        if jump > .9 * cur_max_jump:
            cyc_start = nod
    ninfo[cyc_end]['class'] = 'cyc_l'
    ninfo[cyc_start]['class'] = 'cyc_r'
    print('Guessing that the cycling happens between {} and '
          '{}'.format(cyc_start, cyc_end))

    # approximate ref position on non-reference paths as
    # the position of nearest ref node
    rpos_min = {}
    rpos_max = {}
    for pname in paths:
        cur_min = 0
        cur_max = 0
        for nod in paths[pname]:
            if ninfo[nod]['rpos_min'] != -1:
                # if current node has a ref position
                # update current value
                cur_min = ninfo[nod]['rpos_min'] + ninfo[nod]['size']
            else:
                # if not save the current value
                if nod not in rpos_min:
                    rpos_min[nod] = []
                rpos_min[nod].append(cur_min)
            # same for maximum ref position
            if ninfo[nod]['rpos_max'] != -1:
                cur_max = ninfo[nod]['rpos_max'] + ninfo[nod]['size']
            else:
                if nod not in rpos_max:
                    rpos_max[nod] = []
                rpos_max[nod].append(cur_max)
    # use the median of the approximated positions for each node
    for nod in ninfo:
        if ninfo[nod]['rpos_min'] == -1 and \
           ninfo[nod]['rpos_max'] == -1 and \
           nod in rpos_min:
            ninfo[nod]['rpos_min'] = int(statistics.median(rpos_min[nod]))
            ninfo[nod]['rpos_max'] = int(statistics.median(rpos_max[nod]))

    # count how many time a node is traversed by modules 1 or 2
    if guess_modules:
        # here we first need to guess paths corresponding to the modules
        # split each full path in "module" subpaths
        modules = {}
        for pname in paths:
            modules[pname] = []
            mod_path = []
            for nod in paths[pname]:
                if nod == cyc_end:
                    # start a new module
                    mod_path = [nod]
                elif nod == cyc_start:
                    # end a module
                    modules[pname].append(mod_path)
                    mod_path = []
                elif len(mod_path) > 0:
                    # currently recording a module
                    mod_path.append(nod)
        # update the node info with each module subpath
        for pname in modules:
            if pname == refname:
                # print(modules[pname])
                # easy for the reference, first module is module 1
                # second module is module 2
                for ii, path in enumerate(modules[pname]):
                    if ii > 1:
                        print("Warning: more than two modules in the"
                              " reference?")
                        continue
                    # annotate the traversed nodes appropriately
                    for nod in path:
                        ninfo[nod]['c' + str(ii+1)] += 1
            else:
                for path in modules[pname]:
                    # compare the module to the two reference subpaths
                    # to guess which module it is
                    m1_m = 0
                    m2_m = 0
                    for nod in path:
                        if nod in modules[refname][0]:
                            m1_m += 1
                        if nod in modules[refname][1]:
                            m2_m += 1
                    modn = 'c1'
                    if m2_m > m1_m:
                        modn = 'c2'
                    # annotate the traversed nodes appropriately
                    for nod in path:
                        ninfo[nod][modn] += 1
    else:
        for pname in paths:
            if pname != refname:
                # it's an allele of module 1 or module 2
                modn = pname.split('_')[0]
                if modn not in ['c1', 'c2']:
                    continue
                # annotate the traversed nodes appropriately
                for nod in paths[pname]:
                    ninfo[nod][modn] += 1

    # flag nodes if they are specific to copy 1 or copy 2
    for pos, noden in enumerate(ninfo):
        if ninfo[noden]['class'] != 'none':
            continue
        if ninfo[noden]['c1'] > min_fc * ninfo[noden]['c2']:
            ninfo[noden]['class'] = 'c1'
        if ninfo[noden]['c2'] > min_fc * ninfo[noden]['c1']:
            ninfo[noden]['class'] = 'c2'
        if ninfo[noden]['ref'] > 0 and ninfo[noden]['c1'] == 0 \
           and ninfo[noden]['c2'] == 0:
            ninfo[noden]['class'] = 'ref'

    # pick a reference node for variant sites
    rnodes = {}
    cur_rnode = -1
    for nod in ref:
        if ninfo[nod]['class'] in ['ref', 'none']:
            cur_rnode = nod
        elif nod not in rnodes:
            rnodes[nod] = cur_rnode
    for pname in paths:
        cur_rnode = -1
        for nod in paths[pname]:
            if ninfo[nod]['class'] in ['ref', 'none']:
                cur_rnode = nod
            elif nod not in rnodes:
                rnodes[nod] = cur_rnode
    for nod in ninfo:
        if nod in rnodes:
            ninfo[nod]['rnode'] = rnodes[nod]
        else:
            ninfo[nod]['rnode'] = nod

    # write node info as TSV
    if out_tsv != '':
        outf = open(out_tsv, 'wt')
        outf.write('node\tsize\tseq\tref\tc1\tc2\trpos_min\trpos_max'
                   '\trnode\tclass\n')
        for nod in ninfo:
            ni = ninfo[nod]
            outfmt = '\t'.join(['{}'] * 10) + '\n'
            outf.write(outfmt.format(nod,
                                     ni['size'],
                                     ni['seq'],
                                     ni['ref'],
                                     ni['c1'],
                                     ni['c2'],
                                     ni['rpos_min'],
                                     ni['rpos_max'],
                                     ni['rnode'],
                                     ni['class']))
        outf.close()


def readVariantAnnotation(filen, nodes, offset):
    # for an edge, we'll save the clinvar ID(s) and if it's the ref/alt allele
    # there might be multiple clinvar records sharing an
    # edge (the ref allele edge) so save a list of clinvar info for each edge,
    # e.g. for some:
    # var_edges['nod1']['nod2'] = [{cvid: 'id1', 'ref': True},
    #                              {cvid: 'id2', 'ref': True}]
    var_edges = {}
    # list potential starting nodes for the "variant" edges,
    # i.e. on the reference path and with multiple successors
    snodes = []
    for nod in nodes:
        if 'sucs' not in nodes[nod]:
            continue
        if nodes[nod]['ref'] > 1 and len(nodes[nod]['sucs']) > 1:
            snodes.append(nod)
    # read the input tsv file
    inf = open(filen, 'rt')
    heads = next(inf).rstrip().split('\t')
    for line in inf:
        # parse the line and get some info
        line = line.rstrip().split('\t')
        pos = int(line[heads.index('start')]) - offset
        ref = line[heads.index('ref')]
        alt = line[heads.index('alt')]
        cvid = '{}_{}_{}'.format(line[heads.index('id')],
                                 line[heads.index('nuc.change')],
                                 line[heads.index('prot.change')])
        # find a potential starting node as the one that ends exactly at pos
        snode = ''
        for nod in snodes:
            # 'pos' should be the position of the base at the end of that node
            if nodes[nod]['rpos_max'] + nodes[nod]['size'] == pos:
                snode = nod
        # skip if we didn't find a suitable starting node for the edge
        if snode == '':
            continue
        # check the alleles spelled by that node's successors
        next_nodes = list(nodes[snode]['sucs'].keys())
        for nnod in next_nodes:
            vedge = ''  # is it a ref or alt edge
            # if no padding, looking for exact match in the SNV bubble
            # note: currently doesn't work for MNPs with no padding base
            if ref[0] != alt[0]:
                if ref == nodes[nnod]['seq']:
                    # reference allele of SNV
                    vedge = 'ref'
                elif alt == nodes[nnod]['seq']:
                    # alternate allele of SNV
                    vedge = 'alt'
            else:
                # padding, either insertion or deletion
                # looking for exact match of the long allele
                # assume the short allele is the other one
                # (if on the reference path)
                if alt[1:] == nodes[nnod]['seq'] and nodes[nnod]['ref'] < 2:
                    # alternate allele of insertion
                    vedge = 'alt'
                elif ref[1:] == nodes[nnod]['seq'] and nodes[nnod]['ref'] > 1:
                    print('node {}: ref {}'.format(nnod, nodes[nnod]['ref']))
                    # reference allele of deletion
                    vedge = 'ref'
                elif nodes[nnod]['ref'] > 1:
                    # the sequence doesn't match but if it's a node on the
                    # reference path it might just be the "short" allele with
                    # the deleted sequence (alt) or without the inserted
                    # sequence (ref)
                    if len(ref) < len(alt):
                        # insertion so likely the reference allele
                        vedge = 'ref'
                    else:
                        # deletion so likely the alt allele
                        vedge = 'alt'
            # if we found a match, save it to var_edges
            if vedge == 'ref' or vedge == 'alt':
                if snode not in var_edges:
                    var_edges[snode] = {}
                if nnod not in var_edges[snode]:
                    var_edges[snode][nnod] = []
                var_edges[snode][nnod].append({'cvid': cvid,
                                               'ref': vedge == 'ref'})
    inf.close()
    # find variants with both a ref and alt edge
    esum = {}
    for nod1 in var_edges:
        for nod2 in var_edges[nod1]:
            for vcv in var_edges[nod1][nod2]:
                cvid = vcv['cvid']
                ved = '{}-{}'.format(nod1, nod2)
                if cvid not in esum:
                    esum[cvid] = {'ref': '-', 'alt': '-'}
                if vcv['ref']:
                    esum[cvid]['ref'] = ved
                else:
                    esum[cvid]['alt'] = ved
    matched_variants = []
    for varid in esum:
        if esum[varid]['ref'] != '-' and esum[varid]['alt'] != '-':
            matched_variants.append(varid)
            print('      {}\t{}\t{}'.format(varid,
                                            esum[varid]['ref'],
                                            esum[varid]['alt']))
    print("Matched {} input variants.".format(len(matched_variants)))
    # only return info about those variants
    var_edges_final = {}
    for nod1 in var_edges:
        var_edges_final[nod1] = {}
        for nod2 in var_edges[nod1]:
            var_edges_final[nod1][nod2] = []
            for vcv in var_edges[nod1][nod2]:
                if vcv['cvid'] in matched_variants:
                    var_edges_final[nod1][nod2].append(vcv)
    return (var_edges_final)


def gfaFile(config, fn='', prefix='', check_file=True):
    if fn == '':
        if prefix != '':
            fn = prefix + '.pg.gfa'
        elif 'label' in config:
            fn = config['label'] + '.pg.gfa'
            if not os.path.isfile(fn) and check_file:
                fn = 'parakit.{}.pg.gfa'.format(config['method'])
        else:
            fn = 'parakit.{}.pg.gfa'.format(config['method'])
    if not os.path.isfile(fn) and check_file:
        print('Cannot find/guess GFA file: ' + fn)
    return (fn)


def nodeFile(config, fn='', prefix='', check_file=True):
    if fn == '':
        if prefix != '':
            fn = prefix + '.node_info.tsv'
        elif 'label' in config:
            fn = config['label'] + '.node_info.tsv'
            if not os.path.isfile(fn) and check_file:
                fn = 'parakit.{}.{}'.format(config['method'], 'node_info.tsv')
        else:
            fn = 'parakit.{}.{}'.format(config['method'], 'node_info.tsv')
    if not os.path.isfile(fn) and check_file:
        print('Cannot find/guess node file: ' + fn)
    return (fn)


def clinvarFile(fn, config, check_file=True):
    if fn == '' and 'clinvar' in config:
        fn = config['clinvar']
    if not os.path.isfile(fn) and check_file:
        print('Cannot find ClinVar annotation file: ' + fn +
              '\nEither provide with -a or in the config json '
              '("clinvar" field)')
    return (fn)


def geneFile(fn, config, check_file=True):
    if fn == '' and 'gene_annot' in config:
        fn = config['gene_annot']
    if not os.path.isfile(fn) and check_file:
        print('Cannot find gene annotation file: ' + fn +
              '\nEither provide with -e or in the config json '
              '("gene_annot" field)')
    return (fn)


def prefixFile(config):
    fn = 'parakit.{}'.format(config['method'])
    if 'label' in config:
        fn = config['label']
    return (fn)


def writePathsInfo(paths_res, nodes, stats_fn, info_fn):
    escores_r = paths_res['escores']
    paths = paths_res['paths']
    # write ranked list
    outf_rk = open(stats_fn, 'wt')
    heads_rk = ['hap1', 'hap2', 'cov_cor', 'cov_dev', 'aln_score',
                'cov_cor_adj', 'cov_dev_adj', 'aln_score_adj', 'aln_long_prop']
    ofmt_rk = '\t'.join(['{}'] * len(heads_rk)) + '\n'
    outf_rk.write('\t'.join(heads_rk) + '\n')
    for esc in escores_r:
        outf_rk.write(ofmt_rk.format(esc['hap1'],
                                     esc['hap2'],
                                     esc['cov_cor'],
                                     esc['cov_dev'],
                                     esc['aln_score'],
                                     esc['cov_cor_adj'],
                                     esc['cov_dev_adj'],
                                     esc['aln_score_adj'],
                                     esc['aln_long_prop']))

    # write paths
    outf = open(info_fn, 'wt')
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
