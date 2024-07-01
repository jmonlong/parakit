import gzip
import statistics
import parakit.parakit_class as pkc


def readNodeInfo(filen, verbose=True):
    if verbose:
        print('Reading ' + filen + '...')
    nodes = {}
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
            nodes[noden]['class'] = line[heads.index('class')]
            nodes[noden]['seq'] = line[heads.index('seq')]
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
def readGAF(filen, nodes, verbose=True):
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
    for line in inf_gaf:
        line = line.rstrip().split('\t')
        # if read was seen before, skip supp aln
        readn = line[0]
        if reads_tr.hasRead(readn):
            supp_aln_cpt[readn] += 1
            readn = '{}_sa{}'.format(readn, supp_aln_cpt[readn])
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
            path_cov.reverse()
            path_node_startpos.reverse()
            path_node_endpos.reverse()
            path_read_pos.reverse()
        path_cov_node = [pc[0] for pc in path_cov]
        reads_tr.addRead(readn, path_cov_node, cyc_nodes=cyc_nodes,
                         startpos=path_node_startpos,
                         endpos=path_node_endpos,
                         readpos=path_read_pos)
    inf_gaf.close()
    return (reads_tr)


def updateNodesSucsWithGFA(nodes, filen, verbose=True):
    if verbose:
        print('Reading ' + filen + '...')
    # init sucs for all nodes
    # for nod in nodes:
    #     nodes[nod]['sucs'] = {}
    # read GFA and update edges
    inf = open(filen, 'rt')
    for line in inf:
        line = line.rstrip().split('\t')
        if line[0] == 'L':
            if line[1] in nodes:
                if 'sucs' not in nodes[line[1]]:
                    nodes[line[1]]['sucs'] = {}
                nodes[line[1]]['sucs'][line[3]] = True
    inf.close()


def readGFA(gfa_fn, min_fc=3, refname='grch38', out_tsv=''):
    # node info
    # for each node ID:
    #   size
    #   ref/c1/c2: number of times a ref/c1/c2 path traverses
    #   rpos_min/rpos_max: minimum/maximum position on the reference path
    ninfo = {}
    # reference path
    ref = []
    # copy 1 (c1) and copy 2 (c2) paths
    c1 = []
    c2 = []

    # read GFA and parse info
    inf = open(gfa_fn, 'rt')
    for line in inf:
        line = line.rstrip().split('\t')
        if line[0] == 'P' and line[1] == refname:
            # reference path was added by alignment+augment, hence is a P line
            ref = line[2].replace('+', '').replace('-', '').split(',')
        if line[0] == 'W':
            # otherwise it's an allele of copy 1 or copy 2
            path = line[6].replace('<', '>').split('>')[1:]
            # flip path if mostly traversing in reverse
            if line[6].count('>') < line[6].count('<'):
                path.reverse()
            # add the path to either c1 or c2 path list
            seqn = line[3].split('_')
            if seqn[0] == 'c1':
                c1.append(path)
            elif seqn[0] == 'c2':
                c2.append(path)
        # save some node information like size
        # also init, future fields
        if line[0] == 'S':
            ninfo[line[1]] = {'size': len(line[2]),
                              'ref': 0, 'c1': 0, 'c2': 0,
                              'rpos_min': -1, 'rpos_max': -1,
                              'class': 'none'}
            if len(line[2]) < 10:
                ninfo[line[1]]['seq'] = line[2]
            else:
                ninfo[line[1]]['seq'] = line[2][:10] + '+'
    inf.close()

    # update node info with path cover
    for path in c1:
        for nod in path:
            ninfo[nod]['c1'] += 1
    for path in c2:
        for nod in path:
            ninfo[nod]['c2'] += 1
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

    # flag nodes if they are specific to copy 1 or copy 2
    for pos, noden in enumerate(ninfo):
        if ninfo[noden]['c1'] > min_fc * ninfo[noden]['c2']:
            ninfo[noden]['class'] = 'c1'
        if ninfo[noden]['c2'] > min_fc * ninfo[noden]['c1']:
            ninfo[noden]['class'] = 'c2'
        if ninfo[noden]['ref'] > 0 and ninfo[noden]['c1'] == 0 \
           and ninfo[noden]['c2'] == 0:
            ninfo[noden]['class'] = 'ref'

    # flag nodes contributing to the cycle edge
    # (end of copy 2 to beginning of copy 1)
    cyc_edge_l = ''
    cyc_edge_r = ''
    for nod in ref:
        if ninfo[nod]['rpos_max'] - ninfo[nod]['rpos_min'] > 0:
            cyc_edge_r = nod
            if cyc_edge_l == '':
                cyc_edge_l = nod
    ninfo[cyc_edge_l]['class'] = 'cyc_l'
    ninfo[cyc_edge_r]['class'] = 'cyc_r'

    # approximate ref position on non-reference paths as
    # the position of nearest ref node
    rpos_min = {}
    rpos_max = {}
    for path in c1 + c2:
        cur_min = 0
        cur_max = 0
        for nod in path:
            if ninfo[nod]['rpos_min'] != -1:
                # if current node has a ref position
                # update current value
                if nod not in rpos_min:
                    rpos_min[nod] = []
                rpos_min[nod].append(cur_min)
                cur_min = ninfo[nod]['rpos_min'] + ninfo[nod]['size']
            else:
                # if not save the current value
                if nod not in rpos_min:
                    rpos_min[nod] = []
                rpos_min[nod].append(cur_min)
            # same for maximum ref position
            if ninfo[nod]['rpos_max'] != -1:
                if nod not in rpos_max:
                    rpos_max[nod] = []
                rpos_max[nod].append(cur_max)
                cur_max = ninfo[nod]['rpos_max'] + ninfo[nod]['size']
            else:
                if nod not in rpos_max:
                    rpos_max[nod] = []
                rpos_max[nod].append(cur_max)
    # use the median of the approximated positions for each node
    for nod in ninfo:
        if ninfo[nod]['rpos_min'] == ninfo[nod]['rpos_max'] and \
           nod in rpos_min:
            ninfo[nod]['rpos_min'] = int(statistics.median(rpos_min[nod]))
            ninfo[nod]['rpos_max'] = int(statistics.median(rpos_max[nod]))

    # write node info as TSV
    if out_tsv != '':
        outf = open(out_tsv, 'wt')
        outf.write('node\tsize\tseq\tref\tc1\tc2\trpos_min\trpos_max\tclass\n')
        for nod in ninfo:
            ni = ninfo[nod]
            outfmt = '\t'.join(['{}'] * 9) + '\n'
            outf.write(outfmt.format(nod,
                                     ni['size'],
                                     ni['seq'],
                                     ni['ref'],
                                     ni['c1'],
                                     ni['c2'],
                                     ni['rpos_min'],
                                     ni['rpos_max'],
                                     ni['class']))
        outf.close()


def readVariantAnnotation(filen, nodes, offset):
    vars_pos = {}
    vars_node = {}
    with open(filen, 'rt') as inf:
        heads = next(inf).rstrip().split('\t')
        for line in inf:
            line = line.rstrip().split('\t')
            pos = int(line[heads.index('start')]) - offset
            for nod in nodes:
                if nodes[nod]['rpos_max'] == pos:
                    cvid = '{}_{}_{}'.format(line[heads.index('id')],
                                             line[heads.index('nuc.change')],
                                             line[heads.index('prot.change')])
                    if nodes[nod]['class'] == 'c1' and \
                       nodes[nod]['seq'] == line[heads.index('alt')]:
                        vars_pos[cvid] = pos
                        vars_node[cvid] = nod
                        if 'clinvar' not in nodes[nod]:
                            nodes[nod]['clinvar'] = []
                        nodes[nod]['clinvar'].append(cvid)
                    elif (nodes[nod]['class'] == 'c2') or \
                         (nodes[nod]['seq'] == line[heads.index('ref')]):
                        if 'clinvar_ref' not in nodes[nod]:
                            nodes[nod]['clinvar_ref'] = []
                        nodes[nod]['clinvar_ref'].append(cvid)
    return ({'pos': vars_pos, 'node': vars_node})
