import parakit.parakit_class as pkc
import statistics as stat
import math


def findPaths(nodes, reads, args):
    # upate nodes with read info
    for readn in reads.path:
        path = reads.path[readn]
        for pos, noden in enumerate(path):
            if 'reads' not in nodes[noden]:
                nodes[noden]['reads'] = []
            nodes[noden]['reads'].append([readn, pos])
            if pos < len(path) - 1:
                nodes[noden]['sucs'][path[pos+1]] = True
    # init path list
    paths = {}
    if args.c == 0:
        print('Error: -c must be greater than 0')
        exit(1)
    paths_cls = clusterSubreads(nodes, reads, min_read_support=args.c)
    for walk in paths_cls:
        paths[walk] = paths_cls[walk]

    best_paths = list(paths.keys())

    # evaluate pairs of path
    print('Evaluate quality of pairs of the ' +
          str(len(best_paths)) + ' best paths...')
    # precompute alignment stats for each path first
    aln_score = {}
    for mii in range(len(best_paths)):
        pathn = best_paths[mii]
        aln_score[pathn] = pathReadGraphAlign(paths[pathn],
                                              reads.path, nodes)
    # evaluate each pair
    escores = []
    max_cor = 0
    max_aln = 0
    for mii in range(len(best_paths)):
        for mjj in range(mii, len(best_paths)):
            mode1 = best_paths[mii]
            mode2 = best_paths[mjj]
            path1 = paths[mode1]
            path2 = paths[mode2]
            node_cov = pathNodeCoverage(path1 + path2, nodes)
            read_aln1 = aln_score[mode1]
            read_aln2 = aln_score[mode2]
            esc = evaluatePaths(node_cov, nodes, reads.path,
                                read_aln1, read_aln2)
            esc['hap1'] = mode1
            esc['hap2'] = mode2
            max_cor = max(max_cor, esc['cov_cor'])
            max_aln = max(max_aln, esc['aln_score'])
            escores.append(esc)
    # adjust scores
    for esc in escores:
        esc['cov_cor_adj'] = esc['cov_cor'] / max_cor
        esc['aln_score_adj'] = esc['aln_score'] / max_aln
    # rank hap pairs
    escores_r = sorted(escores,
                       key=lambda k: k['cov_cor_adj'] + k['aln_score_adj'],
                       reverse=True)

    return ({'escores': escores_r, 'paths': paths})


def clusterSubreads(nodes, reads, min_read_support=3, max_cycles=3):
    # init subreads
    sreads = pkc.Subreads()
    sreads.splitReads(reads, nodes)
    # init list of subreads clusters
    sreads_list = [sreads]
    sreads_list_final = []
    # process the list
    while len(sreads_list) > 0:
        csreads = sreads_list.pop(0)
        # look at read coverage and find markers
        csreads.computeCoverage()
        csreads.findMarkers(min_read_support=min_read_support)
        # if some supported markers
        if csreads.nbMarkers() > 0:
            # split the reads in two
            csreads.biClusterReads()
            sreads_list.append(csreads.subsetByCluster(0))
            sreads_list.append(csreads.subsetByCluster(1))
        else:
            # otherwise, save this cluster as a candidate allele
            sreads_list_final.append(csreads)
        # iterate until no clear markers suggesting two alleles
    # enumerate alleles
    res = sreads.enumerateAlleles(sreads_list_final, max_cycles=4,
                                  min_read_support=2*min_read_support)
    return (res)


def pathReadGraphAlign(path, reads, nodes={}, max_node_gap=10):
    # record the node position on the path for each alignment candidate
    path_pos = {}
    for pii, nod in enumerate(path):
        if nod in path_pos:
            path_pos[nod].append(pii)
        else:
            path_pos[nod] = [pii]
    read_scores = {}
    for readn in reads:
        read = reads[readn]
        alns = []
        for rii, nod in enumerate(read):
            # node in read but not in path: skip
            if nod not in path_pos:
                continue
                # node in read but not in path, so record as "-1"
                # if len(alns) == 0:
                #     alns = [-1]
                # else:
                #     for aln in alns:
                #         aln.append(-1)
            # node in both: get the position in path
            for pii in path_pos[nod]:
                if len(alns) == 0:
                    alns.append([pii])
                else:
                    aln_extended = False
                    for aln in alns:
                        if abs(pii - aln[-1]) < max_node_gap:
                            # extend alignment if new position
                            # in path is not too far
                            aln.append(pii)
                            aln_extended = True
                    # if none of the current aln were extended, start a new one
                    if not aln_extended:
                        alns.append([pii])
        # keep alignment with most nodes in common between read and path
        best_score = 0
        for aln in alns:
            aln_score = 0
            for ppi in aln:
                if len(nodes) == 0 or nodes[path[ppi]]['class'] == 'none':
                    aln_score += 1
                else:
                    aln_score += 10
            if aln_score > best_score:
                best_score = aln_score
        max_score = 0
        for nod in read:
            if len(nodes) == 0 or nodes[nod]['class'] == 'none':
                max_score += 1
            else:
                max_score += 10
        read_scores[readn] = round(best_score / max_score, 3)
    return read_scores


def pathNodeCoverage(path, nodes):
    # compute the node coverage in path and reads
    node_cov = {}
    for nod in nodes:
        # if nodes[nod]['class'] != 'none':
        readc = 0
        if 'reads' in nodes[nod]:
            readc = len(nodes[nod]['reads'])
        node_cov[nod] = {'reads': readc, 'path': 0}
    for nod in path:
        # if nodes[nod]['class'] != 'none':
        node_cov[nod]['path'] += 1
    return node_cov


def evaluatePaths(node_cov, nodes, reads, read_aln, read_aln2=[]):
    # correlation between coverage on the predicted path and the reads
    read_c = []
    path_c = []
    for nod in node_cov:
        if nodes[nod]['size'] < 10 and node_cov[nod]['path'] > 0:
            read_c.append(node_cov[nod]['reads'])
            path_c.append(node_cov[nod]['path'])
    if stat.variance(read_c) > 0 and stat.variance(path_c) > 0:
        cov_cor = stat.correlation(read_c, path_c)
    elif stat.variance(read_c) == 0 and stat.variance(path_c) == 0:
        cov_cor = 1
    else:
        cov_cor = .5
    # read alignment on the path
    # find longest reads
    longest_reads = sorted(list(reads.keys()), key=lambda k: -len(reads[k]))
    longest_read_l = len(reads[longest_reads[0]])
    # compute alignment score
    aln_score = 0
    w_sum = 0
    for readn in reads:
        r_score = read_aln[readn]
        if len(read_aln2) > 0:
            r_score = max(read_aln2[readn], r_score)
        read_w = math.exp(20*len(reads[readn])/longest_read_l)
        aln_score += read_w * r_score
        w_sum += read_w
    aln_score = float(aln_score) / w_sum
    # count how many of the top longest reads align well
    rii = 0
    rii_tot = min(len(longest_reads), 20)
    nsupp_reads = 0
    while rii < rii_tot:
        r_score = read_aln[longest_reads[rii]]
        if len(read_aln2) > 0:
            r_score = max(read_aln2[longest_reads[rii]], r_score)
        if r_score > .99:
            nsupp_reads += 1
        rii += 1
    # return scores
    return ({'cov_cor': cov_cor, 'aln_score': aln_score,
             'aln_long_prop': float(nsupp_reads)/rii_tot})
