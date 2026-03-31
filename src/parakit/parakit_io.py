import gzip
import statistics
import parakit.parakit_class as pkc
import os


def readNodeInfo(filen, verbose=False):
    """Loads the node information from a file

    Reads the TSV prepared after the pangenome construction and returns an
    object with information for each node. For example, size, sequence, is
    it on the reference path or informative (assigned to a particular module).

    Args:
        filen : the path to the TSV file
        verbose : print some messages? Default is False

    Returns: dict
                node name -> dict
                                size node size (bp)
                                c1 number of module 1 traversals
                                c2 number of module 2 traversals
                                ref number of reference traversal
                                rpos_min minimum position on reference path
                                rpos_max maximum position on reference path
                                rnode assigned reference node proxy
                                class either none, c1, c2 or ref, cyc_l, cyc_r
                                seq beginning of the node sequence
    """
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
            nodes[noden]['rnode'] = line[heads.index('rnode')]
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


def parseCs(cs):
    nums = [str(ii) for ii in range(10)]
    num_c = ''
    type_c = ''
    res = []
    # we'll consider exact matches and substitutions together, like the M in
    # the CIGAR string.
    # With var_size, we keep track of the size of this M segment
    var_size = 0
    for cc in cs:
        if cc in [':', '*']:
            if type_c != 'M':
                # we were not in a M segment, start new segment
                if type_c != '':
                    # if not the first segment, save and start new one
                    num_c = 0 if num_c == '' else int(num_c)
                    res.append([num_c + var_size, type_c])
                type_c = 'M'
                num_c = ''
                var_size = 0
            else:
                # we're extending a M segment
                if num_c != '':
                    # if we were extending a ":[0-9]+", update the size
                    var_size += int(num_c)
                num_c = ''
            # if we're starting a substitution, we know it will be of size 1
            if cc == '*':
                var_size += 1
        elif cc == '+':
            if type_c != '':
                # if not the first segment, save and start new one
                num_c = 0 if num_c == '' else int(num_c)
                res.append([num_c + var_size, type_c])
            type_c = 'I'
            num_c = ''
            var_size = 0
        elif cc == '-':
            if type_c != '':
                # if not the first segment, save and start new one
                num_c = 0 if num_c == '' else int(num_c)
                res.append([num_c + var_size, type_c])
            type_c = 'D'
            num_c = ''
            var_size = 0
        elif cc in nums:
            # for numbers, update num_c
            num_c += cc
        elif cc in ['A', 'T', 'G', 'C', 'N']:
            # nucleotides mean we're in an insertion or deletion segment
            if type_c in ['I', 'D']:
                var_size += 1
        else:
            print('W: unexpected character in cs tag', cc)
    num_c = 0 if num_c == '' else int(num_c)
    res.append([num_c + var_size, type_c])
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


def readGAF(filen, nodes, verbose=False, min_read_len=0):
    """Read GAF file and load the read alignments

    Process each GAF record, including the CIGAR string, to figure out
    which nodes are traversed (on at least one base) in the
    graph. Reads are flipped so that they (mostly) traverse the
    pangenome in the forward orientation. The nodes input is used to
    know the node size.

    Args:
        filen : path to the GAF file
        nodes : the dict with node information (inc. size)
        verbose : print some informative messages? Default: False

    Returns: a Reads object

    """
    if verbose:
        print('Reading ' + filen + '...')
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
            nsupp_aln += 1
            continue
        else:
            supp_aln_cpt[readn] = 0
        # skip if read is too short
        if int(line[1]) < min_read_len:
            continue
        path = parsePath(line[5])
        # find CIGAR string
        cg = ''
        for tag in line[12:]:
            tag = tag.split(':')
            if tag[0] == 'cg':
                cg = parseCg(tag[2])
                break
        if cg == '':
            # look for cs tag instead
            for tag in line[12:]:
                tag_s = tag.split(':')
                if tag_s[0] == 'cs':
                    cg = parseCs(tag[5:])
                    break
            if cg == '':
                nunmp += 1
                continue
        # parse CIGAR string
        n_idx = 0
        node_pos = int(line[7])
        read_pos = int(line[2])
        path_cov = []
        path_node_startpos = []
        path_node_endpos = []
        path_read_pos = []
        cur_endpos = node_pos
        for cc in cg:
            for bp in range(cc[0]):
                if cc[1] in ['=', 'M']:
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
                if cc[1] in ['I', 'X', '=', 'M']:
                    read_pos += 1
                if cc[1] in ['D', 'X', '=', 'M']:
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
                spos = path_node_startpos[ii]
                path_node_startpos[ii] = path_node_endpos[ii]
                path_node_endpos[ii] = spos
        path_cov_node = [pc[0] for pc in path_cov]
        reads_tr.addRead(readn, path_cov_node,
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
    """Update the node dict with information about next nodes in graph

    Reads a GFA file and saves successors of each nodes. Adds a new
    'sucs' key in the node dict. It's a dict (used as a set) with the
    names of the successor nodes.

    Args:
        nodes : dict with nodes information
    """
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
    cycled = False
    for nod in ref:
        ninfo[nod]['ref'] += 1
        # also minimum/maximum position on reference
        if ninfo[nod]['rpos_min'] == -1 and not cycled:
            # first time the reference traverses this node
            ninfo[nod]['rpos_min'] = refpos
        else:
            # second time, fill rpos_max
            ninfo[nod]['rpos_max'] = refpos
            cycled = True
        refpos += ninfo[nod]['size']
    # assign position on the flanks
    ii = 0
    while ninfo[ref[ii]]['rpos_max'] == -1:
        if ninfo[ref[ii]]['rpos_min'] != 1:
            ninfo[ref[ii]]['rpos_max'] = ninfo[ref[ii]]['rpos_min']
        else:
            break
        ii += 1
    ii = len(ref) - 1
    while ninfo[ref[ii]]['rpos_min'] == -1:
        if ninfo[ref[ii]]['rpos_max'] != 1:
            ninfo[ref[ii]]['rpos_min'] = ninfo[ref[ii]]['rpos_max']
        else:
            break
        ii += -1

    # approximate ref position on non-reference paths as
    # the position of nearest ref node
    rpos_min = {}
    rpos_max = {}
    for pname in paths:
        cur_min = 0
        cur_max = None
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
            elif cur_max is not None:
                if nod not in rpos_max:
                    rpos_max[nod] = []
                rpos_max[nod].append(cur_max)
    # use the median of the approximated positions for each node
    for nod in ninfo:
        if ninfo[nod]['rpos_min'] == -1 and \
           nod in rpos_min:
            if nod == '1279':
                print(ninfo[nod])
            ninfo[nod]['rpos_min'] = int(statistics.median(rpos_min[nod]))
        if ninfo[nod]['rpos_max'] == -1 and \
           nod in rpos_max:
            ninfo[nod]['rpos_max'] = int(statistics.median(rpos_max[nod]))

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


def gfaFile(config, fn='', prefix='', check_file=True):
    """Returns the path to the most appropriate GFA

    Checks, in that order, the provided filename, the provided prefix, the
    config file.

    Args:
        config : dict with the configuration values (e.g. from input JSON)
        fn : filename
        prefix : label or prefix for this locus/experiment
        check_file : check if the file exists?
    """
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
    """Returns the path to the most appropriate node info TSV

    Checks, in that order, the provided filename, the provided prefix, the
    config file.

    Args:
        config : dict with the configuration values (e.g. from input JSON)
        fn : filename
        prefix : label or prefix for this locus/experiment
        check_file : check if the file exists?
    """
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
    """Returns the path to the most appropriate ClinVar annotation TSV

    Checks, in that order, the provided filename, the config file.

    Args:
        config : dict with the configuration values (e.g. from input JSON)
        fn : filename
        check_file : check if the file exists?
    """
    if fn == '' and 'clinvar' in config:
        fn = config['clinvar']
    if not os.path.isfile(fn) and check_file:
        print('Cannot find ClinVar annotation file: ' + fn +
              '\nEither provide with -a or in the config json '
              '("clinvar" field)')
    return (fn)


def geneFile(fn, config, check_file=True):
    """Returns the path to the most appropriate gene annotation file

    Checks, in that order, the provided filename, the config file.

    Args:
        config : dict with the configuration values (e.g. from input JSON)
        fn : filename
        check_file : check if the file exists?
    """
    if fn == '' and 'gene_annot' in config:
        fn = config['gene_annot']
    if not os.path.isfile(fn) and check_file:
        print('Cannot find gene annotation file: ' + fn +
              '\nEither provide with -e or in the config json '
              '("gene_annot" field)')
    return (fn)


def prefixFile(config):
    """Prepare a prefix for a project based on the config file

    If label is not defined in the config JSON, makes sure we use a
    consistent default.

    Args:
        config : dict with the configuration values
    """
    fn = 'parakit.{}'.format(config['method'])
    if 'label' in config:
        fn = config['label']
    return (fn)


def writePathsInfo(paths_res, nodes, stats_fn, info_fn):
    """Write inferred path information

    Add some node information to a list of predicted paths and
    evaluation metrics, and write them in two files: a "stats" file
    with the diplotype pair evaluation, and a "info" file with the
    nodes traversed by each path.

    Args:
        paths_res : dict with escores and paths produced by findPaths
        nodes : dict with node information
        stats_fn: output file name for the "stats" file (diplotype evaluation)
        info_fn: output file name for the "info" file with path definition.

    """
    escores_r = paths_res['escores']
    paths = paths_res['paths']
    # write ranked list
    outf_rk = open(stats_fn, 'wt')
    heads_rk = ['hap1', 'hap2', 'cov_cor', 'cov_dev', 'hap_ll', 'aln_score',
                'cov_cor_adj', 'cov_dev_adj', 'hap_ll_adj',
                'aln_score_adj', 'aln_long_prop']
    ofmt_rk = '\t'.join(['{}'] * len(heads_rk)) + '\n'
    outf_rk.write('\t'.join(heads_rk) + '\n')
    for esc in escores_r:
        outf_rk.write(ofmt_rk.format(esc['hap1'],
                                     esc['hap2'],
                                     esc['cov_cor'],
                                     esc['cov_dev'],
                                     esc['hap_ll'],
                                     esc['aln_score'],
                                     esc['cov_cor_adj'],
                                     esc['cov_dev_adj'],
                                     esc['hap_ll_adj'],
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


def readDiplotype(stats_fn, info_fn, verbose=False):
    # read the stats file and extract the names of the best diplotype
    hap_names = []
    inf = open(stats_fn, 'rt')
    for line in inf:
        line = line.rstrip().split('\t')
        if line[0] != 'hap1':
            hap_names.append(line[0])
            hap_names.append(line[1])
            break
    inf.close()
    # read info and save path of the two haplotypes
    dip_paths = {}
    for hn in hap_names:
        dip_paths[hn] = []
    inf = open(info_fn, 'rt')
    for line in inf:
        line = line.rstrip().split('\t')
        if line[0] in hap_names:
            dip_paths[line[0]].append(line[1])
    inf.close()
    # print info
    if verbose:
        for hn in hap_names:
            print('{}: {} nodes.'.format(hn, len(dip_paths[hn])))
    return dip_paths


def readGeneAnnotation(gene_fn, nodes, offset, gene_list=[]):
    # read TSV file, keep only 'gene' and 'exon' line
    inf = open(gene_fn, 'rt')
    header = None
    elts = []
    for line in inf:
        line = line.rstrip().split('\t')
        if header is None:
            header = {}
            for idx, col in enumerate(line):
                header[col] = idx
            continue
        # skip genes not in gene list
        if len(gene_list) > 0 and line[header['gene_name']] not in gene_list:
            continue
        # keep if exon or gene record
        if line[header['type']] in ['exon', 'gene']:
            elt = {}
            for col in header:
                if col in ['start', 'end']:
                    line[header[col]] = int(line[header[col]])
                elt[col] = line[header[col]]
            elts.append(elt)
    inf.close()
    # add start/end nodes
    for node in nodes:
        rpmin = nodes[node]['rpos_min'] + offset
        rpmax = nodes[node]['rpos_max'] + offset
        nsize = nodes[node]['size']
        for elt in elts:
            if (elt['start'] >= rpmin and
                    elt['start'] < rpmin + nsize):
                elt['node_start'] = node
                elt['mod_start'] = 'c1'
            elif (elt['start'] >= rpmax and
                    elt['start'] < rpmax + nsize):
                elt['node_start'] = node
                elt['mod_start'] = 'c2'
            if (elt['end'] >= rpmin and
                    elt['end'] < rpmin + nsize):
                elt['node_end'] = node
                elt['mod_end'] = 'c1'
            elif (elt['end'] >= rpmax and
                    elt['end'] < rpmax + nsize):
                elt['node_end'] = node
                elt['mod_end'] = 'c2'
    return elts
