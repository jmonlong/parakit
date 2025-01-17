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
            # record which read and at which position it traverses this node
            nodes[noden]['reads'].append([readn, pos])
            # just in case, also make sure edges followed
            # by the read are recorded
            if pos < len(path) - 1:
                nodes[noden]['sucs'][path[pos+1]] = True
    if args.c == 0:
        print('Error: minimum read support -c must be greater than 0')
        exit(1)
    # enumerate path candidates by clustering (sub)reads
    paths = clusterSubreads(nodes, reads, min_read_support=args.c)
    # if too many paths, rerun with more stringent read support
    while len(paths) > 500:
        args.c += 1
        print("Too many paths ({}). Rerunning with min read support: "
              "{}".format(len(paths), args.c))
        paths = clusterSubreads(nodes, reads, min_read_support=args.c)

    # potentially first select the best paths (in case there are too
    # many paths and we don't want to consider all pairs)
    # for now, we keep them all
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
    # sort the reads by length
    longest_reads = sorted(list(reads.path.keys()),
                           key=lambda k: -len(reads.path[k]))
    # precompute read coverage on each node
    read_cov = {}
    for nod in nodes:
        readc = 0
        if 'reads' in nodes[nod]:
            readc = len(nodes[nod]['reads'])
        read_cov[nod] = readc
    # precompute path coverage on each node
    path_cov = {}
    for pathn in best_paths:
        path_cov[pathn] = {}
        for nod in paths[pathn]:
            if nod not in path_cov[pathn]:
                path_cov[pathn][nod] = 1
            else:
                path_cov[pathn][nod] += 1
    # evaluate each pair
    escores = []
    # keep track of some min/max values to adjust scores later
    max_cor = 0
    max_cov_dev = 0
    min_cov_dev = math.inf
    max_aln = 0
    for mii in range(len(best_paths)):
        for mjj in range(mii, len(best_paths)):
            mode1 = best_paths[mii]
            mode2 = best_paths[mjj]
            esc = evaluatePaths(read_cov, path_cov[mode1], path_cov[mode2],
                                nodes, reads.path,
                                longest_reads,
                                aln_score[mode1], aln_score[mode2])
            esc['hap1'] = mode1
            esc['hap2'] = mode2
            escores.append(esc)
            # update the extreme values used to scale score later
            max_cor = max(max_cor, esc['cov_cor'])
            max_cov_dev = max(max_cov_dev, esc['cov_dev'])
            min_cov_dev = min(min_cov_dev, esc['cov_dev'])
            max_aln = max(max_aln, esc['aln_score'])
    # adjust scores to the [0,1] range
    cov_dev_a = (max_cov_dev - min_cov_dev)
    for esc in escores:
        esc['cov_cor_adj'] = esc['cov_cor'] / max_cor
        esc['cov_dev_adj'] = (max_cov_dev - esc['cov_dev']) / cov_dev_a
        esc['aln_score_adj'] = esc['aln_score'] / max_aln
    # rank hap pairs
    escores_r = sorted(escores,
                       key=lambda k: k['cov_dev_adj'] + k['aln_score_adj'],
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
        # to make more candidates and lower the risk of clustering too much,
        # we also save intermediate clusters
        # to be used for path enumeration later.
        sreads_list_final.append(csreads)
        # look at read coverage and find markers
        csreads.computeCoverage()
        csreads.findMarkers(min_read_support=min_read_support)
        # if some supported markers
        if csreads.nbMarkers() > 0:
            # split the reads in two
            csreads.biClusterReads()
            sreads_list.append(csreads.subsetByCluster(0))
            sreads_list.append(csreads.subsetByCluster(1))
        # iterate until no clear markers suggesting two alleles
    # enumerate alleles
    res = sreads.enumerateAlleles(sreads_list_final, max_cycles=4,
                                  min_read_support=min_read_support)
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
        if max_score == 0:
            read_scores[readn] = 0
        else:
            read_scores[readn] = round(best_score / max_score, 3)
    return read_scores


def evaluatePaths(read_cov, path_cov_1, path_cov_2,
                  nodes, reads,
                  longest_reads, read_aln, read_aln2=[]):
    # correlation between coverage on the predicted path and the reads
    read_c = []
    path_c = []
    # add nodes in both paths or unique to first path
    for nod in path_cov_1:
        if nodes[nod]['size'] < 10:
            read_c.append(read_cov[nod])
            if nod in path_cov_2:
                path_c.append(path_cov_1[nod] + path_cov_2[nod])
            else:
                path_c.append(path_cov_1[nod])
    # add nodes unique to second path
    for nod in path_cov_2:
        if nodes[nod]['size'] < 10 and nod not in path_cov_1:
            read_c.append(read_cov[nod])
            path_c.append(path_cov_2[nod])
    if stat.variance(read_c) > 0 and stat.variance(path_c) > 0:
        cov_cor = stat.correlation(read_c, path_c)
    elif stat.variance(read_c) == 0 and stat.variance(path_c) == 0:
        cov_cor = 1
    else:
        cov_cor = .5
    # deviation from copy-number normalized counts
    one_copy_cov = []
    for ii in range(len(read_c)):
        if path_c[ii] > 0:
            one_copy_cov.append(float(read_c[ii]) / path_c[ii])
    one_copy_cov = stat.median(one_copy_cov)
    cov_dev = []
    for ii in range(len(read_c)):
        cov_dev.append(abs(read_c[ii] - path_c[ii] * one_copy_cov))
    cov_dev = stat.mean(cov_dev)
    # read alignment on the path
    # compute alignment score
    aln_score = 0
    w_sum = 0
    for readn in reads:
        r_score = read_aln[readn]
        if len(read_aln2) > 0:
            r_score = max(read_aln2[readn], r_score)
        aln_score += r_score
        w_sum += 1
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
    return ({'cov_cor': cov_cor, 'cov_dev': cov_dev,
             'aln_score': aln_score,
             'aln_long_prop': float(nsupp_reads)/rii_tot})
