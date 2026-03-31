import parakit.parakit_path_original as pkpo
import statistics as stat
import math
import random
# import networkx as nx

DEBUG = False


def findPaths(nodes, reads, args):
    """Find paths through the pangenome from read alignment

    This is the main function that will call the different approaches to
    enumerate the best paths for the input reads. Here "path" is equivalent
    to haplotype, and in the end we'll want to enumerate a pair of them
    (a diplotype).

    Args:
        nodes : dict with nodes information (size, class, etc)
        reads : a Reads object
        args : argparse arguments given to parakit

    Returns: dict
                escores : list of dict
                                    hap1, hap2, hap_ll, cov_cor, aln_score
                paths : dict path name -> list of nodes
    """
    if args.c == 0:
        print('Error: minimum read support -c must be greater than 0')
        exit(1)
    if args.t:
        print("Enumerating candidate haplotypes..")
        print("\tMinimum read support: {}".format(args.c))
    # upate nodes with read info
    for readn in reads.path:
        path = reads.path[readn]
        for pos, noden in enumerate(path):
            if 'reads' not in nodes[noden]:
                nodes[noden]['reads'] = []
                # record which read and at which position it traverses
                # this node
            nodes[noden]['reads'].append([readn, pos])
            # just in case, also make sure edges followed
            # by the read are recorded
            if pos < len(path) - 1:
                nodes[noden]['sucs'][path[pos+1]] = True
                # two modes to cluster the reads in module-alleles by
                # default, the first one implemented iteratively
                # splits clusters in two
    args_m = args.m if args.m else '50,50'
    args_m = args_m.split(',')
    cl_o = None
    if len(args_m) == 2:
        # once at most max_cls candidate clusters are enumerated
        # they are combined into at most max_haps haplotypes
        max_cls = int(args_m[0])
        max_haps = int(args_m[1])
        # first cluster subreads
        cl_o = pkpo.biClusterSubreads(nodes, reads,
                                      min_read_support=args.c,
                                      max_cls=max_cls,
                                      max_cycles=4,
                                      max_haps=max_haps,
                                      verbose=args.t)
    else:
        # this is the EM algorithm that will test different module copy
        # number (from min_modc, to max_modc). n_em_cand repetitions will
        # be run for each copy number
        min_modc = int(args_m[0])
        max_modc = int(args_m[1])
        n_em_cand = int(args_m[2])
        cl_o = clusterSubreadsEM(nodes, reads,
                                 min_read_support=args.c,
                                 min_module_copies=min_modc,
                                 max_module_copies=max_modc,
                                 ncandidates=n_em_cand,
                                 verbose=args.t)

    # evaluate pairs of path
    if args.t:
        print('Evaluating quality of {} '
              'diplotypes...'.format(len(cl_o['hap_pairs'])))
        # precompute alignment stats for each path first
    if args.t:
        print('\tPrecomputing alignment of reads to each haplotype...')
    aln_score = {}
    for pathn in cl_o['paths']:
        aln_score[pathn] = pathReadGraphAlign(cl_o['paths'][pathn],
                                              reads.path, nodes)
    # sort the reads by length
    longest_reads = sorted(list(reads.path.keys()),
                           key=lambda k: -len(reads.path[k]))
    # precompute read coverage on each node
    if args.t:
        print('\tPrecomputing node coverage for reads and haplotypes...')
    read_cov = {}
    for nod in nodes:
        readc = 0
        if 'reads' in nodes[nod]:
            readc = len(nodes[nod]['reads'])
            read_cov[nod] = readc
            # precompute path coverage on each node
    path_cov = {}
    for pathn in cl_o['paths']:
        path_cov[pathn] = {}
        for nod in cl_o['paths'][pathn]:
            if nod not in path_cov[pathn]:
                path_cov[pathn][nod] = 1
            else:
                path_cov[pathn][nod] += 1
                # evaluate each pair
    if args.t:
        print('\tEvaluating each haplotype pair...')
    escores = []
    # keep track of some min/max values to adjust scores later
    max_cor = 0
    max_cov_dev = 0
    min_cov_dev = math.inf
    max_aln = 0
    min_hap_ll = 0
    max_hap_ll = -math.inf
    # evaluate each haplotype pair
    for hpair_names in cl_o['hap_pairs']:
        mode1 = hpair_names[0]
        mode2 = hpair_names[1]
        esc = evaluatePaths(read_cov, path_cov[mode1], path_cov[mode2],
                            nodes, longest_reads,
                            aln_score[mode1], aln_score[mode2])
        esc['hap1'] = mode1
        esc['hap2'] = mode2
        escores.append(esc)
        # update the extreme values used to scale score later
        max_cor = max(max_cor, esc['cov_cor'])
        max_cov_dev = max(max_cov_dev, esc['cov_dev'])
        min_cov_dev = min(min_cov_dev, esc['cov_dev'])
        max_aln = max(max_aln, esc['aln_score'])
        max_hap_ll = max(max_hap_ll, esc['hap_ll'])
        min_hap_ll = min(min_hap_ll, esc['hap_ll'])
        # adjust scores to the [0,1] range
    cov_dev_a = (max_cov_dev - min_cov_dev)
    hap_ll_a = (min_hap_ll - max_hap_ll)
    for esc in escores:
        esc['cov_cor_adj'] = esc['cov_cor'] / max_cor
        if cov_dev_a == 0:
            esc['cov_dev_adj'] = 1
        else:
            esc['cov_dev_adj'] = (max_cov_dev - esc['cov_dev']) / cov_dev_a
            # esc['cov_dev_adj'] = 1 - esc['cov_dev'] / max_cov_dev
        if hap_ll_a == 0:
            esc['hap_ll_adj'] = 1
            esc['aln_score_adj'] = 1
        else:
            esc['hap_ll_adj'] = (min_hap_ll - esc['hap_ll']) / hap_ll_a
            esc['aln_score_adj'] = esc['aln_score'] / max_aln
            # rank hap pairs
    escores_r = sorted(escores,
                       key=lambda k: (k['cov_dev_adj']
                                      + k['hap_ll_adj']
                                      + k['cov_cor_adj']
                                      + 3 * k['aln_score_adj']),
                       reverse=True)
    return ({'escores': escores_r, 'paths': cl_o['paths']})


def clusterSubreadsEM(nodes, reads, min_read_support=3, min_module_copies=1,
                      max_module_copies=6, ncandidates=3, verbose=False):
    """Run the EM subread clustered on a range of parameters

    Cluster subreads using the EM approach. Test multiple parameters
    for the number of clusters (from min_module_copies to
    max_module_copies).

    Args:
        nodes : dict with node information
        reads : Reads object
        min_read_support : minimum read support to define markers?
        min_module_copies : minimum number of copies N to consider
        max_module_copies : maximum number of copies N to consider
        ncandidates : number of candidates generated for each N
        verbose : show informative message. Default: False

    Returns: a dict with 'paths' -> path for each haplotype
                         'hap_pairs' -> list of list with haplotype names
    """
    module_copies = list(range(min_module_copies, max_module_copies + 1))
    # init subreads
    sreads = Subreads()
    sreads.splitReads(reads, nodes)
    # prepare markers
    marker_finder = MarkerFinder(nodes)
    marker_finder.findMarkers()
    sreads.profileSubreads(marker_finder, min_read_support=min_read_support)
    if verbose:
        print('\t\tClustering subreads...')
    # list of subreads clusters
    # try to cluster in K profiles, K being any possible module number
    paths = {}
    hap_pair_names = []
    for nmodules in module_copies:
        # run the greedy clusterer?
        if ncandidates >= 0:
            if verbose:
                print('\t\t\t{} module(s) (greedy)...'.format(nmodules))
            sreads.clearClusters()
            sreads.clusterSubreadsGreedy(nmodules)
            if verbose:
                sreads.print()
            cpaths = sreads.threadDiplotype(nodes)
            hpair_names = []
            for pathn in cpaths:
                new_pathn = '{}m_greedy_{}'.format(nmodules, pathn)
                paths[new_pathn] = cpaths[pathn]
                hpair_names.append(new_pathn)
            if len(hpair_names) > 0:
                hap_pair_names.append(hpair_names)
            if verbose:
                print('\t\t\t{} module(s) (phase block merging)'
                      '...'.format(nmodules))
            sreads.clearClusters()
            sreads.clusterSubreadsPhaseBlockMerging(nmodules)
            if verbose:
                sreads.print()
            cpaths = sreads.threadDiplotype(nodes)
            hpair_names = []
            for pathn in cpaths:
                new_pathn = '{}m_pb_{}'.format(nmodules, pathn)
                paths[new_pathn] = cpaths[pathn]
                hpair_names.append(new_pathn)
            if len(hpair_names) > 0:
                hap_pair_names.append(hpair_names)
        # skip EM if candidate number is 0
        if ncandidates == 0:
            continue
        # let's cluster the subreads multiple times to get
        # more than one candidate
        for ncand in range(ncandidates):
            if verbose:
                print('\t\t\t{} module(s) ({}/{})...'.format(nmodules,
                                                             ncand + 1,
                                                             ncandidates))
            sreads.clearClusters()
            sreads.clusterSubreads(nmodules)
            if verbose:
                sreads.print()
            cpaths = sreads.threadDiplotype(nodes)
            hpair_names = []
            for pathn in cpaths:
                new_pathn = '{}m_c{}_{}'.format(nmodules, ncand, pathn)
                paths[new_pathn] = cpaths[pathn]
                hpair_names.append(new_pathn)
            if len(hpair_names) > 0:
                hap_pair_names.append(hpair_names)
    # return results
    res = {}
    res['paths'] = paths
    res['hap_pairs'] = hap_pair_names
    return (res)


def computeHapPairLikelihood(nodes, path_cov_1, path_cov_2,
                             read_scores_1, read_scores_2,
                             prob_error=.01):
    """Compute the likelihood of read assignments to a diplotype

    For each node, compute the likelihood of the read assignment to
    each haplotype (expected 50/50 chance). The probabily of error
    represent how likely it is to see a node in the read but not in
    the haplotype.

    Args:
        nodes : dict with node info
        path_cov_1 : node coverage of first haplotype
        path_cov_2 : node coverage of second haplotype
        read_scores_1 : alignment score of reads on first haplotype
        read_scores_2 : alignment score of reads on second haplotype
        prob_error : probability of a node in the read but not in the haplotype

    Returns: the log-likelihood of the read assignment to the diplotype

    """
    # for each node, compute a likelihood of observing this path assignment
    tot_ll = 0
    for noden in nodes:
        # skip nodes with no read coverage
        if 'reads' not in nodes[noden]:
            continue
        # otherwise count number of reads assigned to each path
        read_to_path = [0, 0]
        for read_pos in nodes[noden]['reads']:
            if read_scores_1[read_pos[0]] == read_scores_2[read_pos[0]]:
                # tie, assign to both?
                read_to_path[0] += 1
                read_to_path[1] += 1
            elif read_scores_1[read_pos[0]] > read_scores_2[read_pos[0]]:
                # assign to first path
                read_to_path[0] += 1
            else:
                # assign to second path
                read_to_path[1] += 1
        tot_reads = sum(read_to_path)
        # compare with number in each path
        nc1 = path_cov_1[noden] if noden in path_cov_1 else 0
        nc2 = path_cov_2[noden] if noden in path_cov_2 else 0
        # compute log-likelihood
        ll = None
        if nc1 + nc2 == 0:
            # likelihood of mapping errors
            ll = tot_reads * math.log(prob_error)
        else:
            # binomial distribution across pair of paths
            ll = math.log(math.comb(tot_reads, read_to_path[0]))
            prob1 = nc1 / (nc1 + nc2) if nc1 > 0 else prob_error
            prob2 = nc2 / (nc1 + nc2) if nc2 > 0 else prob_error
            ll += read_to_path[0] * math.log(prob1)
            ll += read_to_path[1] * math.log(prob2)
        tot_ll += ll
    # return the sum of log-likelihood across all nodes
    return (tot_ll)


def evaluatePaths(read_cov, path_cov_1, path_cov_2, nodes, longest_reads,
                  read_aln, read_aln2=[]):
    """Compute evaluation metrics for this diplotype

    Compure coverage, alignment, and assignment likelihood for this diplotype.

    Args:
        read_cov : node coverage across all reads
        path_cov_1 : node coverage in the first haplotype
        path_cov_2 : node coverage in the second haplotype
        nodes : dict with node information
        longest_reads : list of longes reads (or sorted by decreasing length)
        read_aln : alignment score of all reads on first haplotype
        read_aln2 : alignment score of all reads on second haplotype

    Returns: a dict with multiple evaluation metrics
    """
    # prepare vectors of read and node coverage across selected nodes
    read_c = []
    path_c = []
    for nod in nodes:
        if nodes[nod]['size'] >= 5:
            # keep only small nodes, to keep it simple and focus on important
            # nodes (those that are more likely copy-number variable)
            continue
        if nod not in read_cov:
            continue
        read_c.append(read_cov[nod])
        pcov = 0
        if nod in path_cov_1:
            pcov += path_cov_1[nod]
        if nod in path_cov_2:
            pcov += path_cov_2[nod]
        path_c.append(pcov)
    # correlation between coverage on the predicted path and the reads
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
    max_aln_score = 0
    w_sum = 0
    for readn in read_aln:
        r_score = read_aln[readn]
        if len(read_aln2) > 0:
            r_score = max(read_aln2[readn], r_score)
        aln_score += r_score
        max_aln_score = max(max_aln_score, r_score)
        w_sum += 1
    aln_score = float(aln_score) / w_sum
    # count how many of the top longest reads align well
    aln_score_th = aln_score + .8 * (max_aln_score - aln_score)
    rii = 0
    rii_tot = min(len(longest_reads), len(read_aln)/10)
    nsupp_reads = 0
    while rii < rii_tot:
        r_score = read_aln[longest_reads[rii]]
        if len(read_aln2) > 0:
            r_score = max(read_aln2[longest_reads[rii]], r_score)
        if r_score > aln_score_th:
            nsupp_reads += 1
        rii += 1
    # compute likelihood of read assignment to haplotypes
    hap_ll = computeHapPairLikelihood(nodes, path_cov_1, path_cov_2,
                                      read_aln, read_aln2)
    # return scores
    return ({'cov_cor': cov_cor, 'cov_dev': cov_dev,
             'aln_score': aln_score, 'hap_ll': hap_ll,
             'aln_long_prop': float(nsupp_reads)/rii_tot})


def pathReadGraphAlign(path, reads, nodes={}, max_node_gap=10):
    """Compute alignment scores for all reads on a haplotype path

    Find alignment for each read on the path and save a score. The
    score is computed as the best score divided by the maximum
    possible score (perfect alignment). If all nodes are weighted
    equally, the score is the proportion of nodes matching between
    aligned read and the path. Informative nodes (class: c1, c2 or
    ref) can be weighted up (10 vs 1), if a 'nodes' argument is
    provided.

    Args:
        path : path as a list of nodes
        reads : read path as a list of nodes
        nodes : Optional. Dict with node information to weight nodes
        max_node_gap : maximum allowed gap (in number of nodes) allowed

    Returns: a dict: read name -> alignment score

    """
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


class Subread:
    """Subread object with basic information

    Organize basic informaiton about a subread.

    Attributes:
        path : list of nodes traversed by the subread
        name : subread name
        cluster : assigned cluster
        cluster_conf : confidence in the assigned cluster
        in_module : is it a subread within the module? Otherwise a flank

    Methods:
        assignCluster : assign to a cluster with a certain confidence
    """
    def __init__(self, name, path, read_name):
        self.name = name
        self.read_name = read_name
        self.path = path
        self.cluster = None
        self.cluster_conf = 0
        self.in_module = True

    def assignCluster(self, cluster, confidence):
        self.cluster = cluster
        self.cluster_conf = confidence

    def clearCluster(self):
        if self.in_module:
            self.cluster = None
            self.cluster_conf = 0

    def makeObs(self, marker_finder):
        """"Prepare list of position and alleles for the subread"""
        # check the marker for each edge
        pos_val = {}
        for srpos in range(len(self.path)-1):
            edge = '{}_{}'.format(self.path[srpos], self.path[srpos+1])
            marks = marker_finder.getMarkers(edge)
            for edgem in marks:
                pos_val[edgem] = marks[edgem]
        obs = Observation(pos_val)
        return obs

    def makeObsPos(self, marker_finder):
        """"Prepare list of position and alleles for the subread"""
        # check the marker for each edge
        pos_val = {}
        pos_absent = set()
        for srpos in range(len(self.path)-1):
            edge = '{}_{}'.format(self.path[srpos], self.path[srpos+1])
            marks = marker_finder.getMarkers(edge)
            for posn in marks:
                if marks[posn] == 'absent':
                    pos_absent.add(posn)
                else:
                    pos_val[posn] = marks[posn]
        # mark absent positions if not already set to a marker
        for posn in pos_absent:
            if posn not in pos_val:
                pos_val[posn] = 'absent'
        obs = Observation(pos_val)
        return obs


# class NetworkBiClusterer:
#     def __init__(self):
#         self.edge_count = {}
#         self.vertex_col = {}
#         self.edge_to_vert = {}

#     def print(self):
#         blocks = [str(bb) for bb in self.vertex_col]
#         print('NetworkBiClusterer with {} blocks'
#               ': {}.'.format(len(blocks), ', '.join(blocks)))

#     def addNode(self, elt):
#         if elt not in self.vertex_col:
#             self.vertex_col[elt] = len(self.vertex_col)

#     def addEdge(self, elt1, elt2):
#         ename = '{}_{}'.format(elt1, elt2)
#         if elt1 > elt2:
#             ename = '{}_{}'.format(elt2, elt1)
#         if ename not in self.edge_to_vert:
#             self.edge_to_vert[ename] = [elt1, elt2]
#             self.edge_count[ename] = 0
#         self.edge_count[ename] += 1
#         if elt1 not in self.vertex_col:
#             self.vertex_col[elt1] = len(self.vertex_col)
#         if elt2 not in self.vertex_col:
#             self.vertex_col[elt2] = len(self.vertex_col)

#     def biCluster(self):
#         ncolors = len(self.vertex_col)
#         # sort edges by counts
#         edge_count_s = sorted(self.edge_count, reverse=True)
#         # iteratively merge
#         cur_edge_ii = 0
#         while ncolors > 2 and cur_edge_ii < len(edge_count_s):
#             # combine next most supported edge
#             elts = self.edge_to_vert[edge_count_s[cur_edge_ii]]
#             if self.vertex_col[elts[0]] != self.vertex_col[elts[1]]:
#                 # change the color of the second element
#                 col_mod = self.vertex_col[elts[1]]
#                 new_col = self.vertex_col[elts[0]]
#                 for elt in self.vertex_col:
#                     if self.vertex_col[elt] == col_mod:
#                         self.vertex_col[elt] = new_col
#                 ncolors += -1
#             cur_edge_ii += 1
#         # TODO ideally recompute the support after each merge
#         # make two lists
#         colors = set()
#         for elt in self.vertex_col:
#             colors.add(self.vertex_col[elt])
#         res = []
#         for col in colors:
#             col_elts = set()
#             for elt in self.vertex_col:
#                 if self.vertex_col[elt] == col:
#                     col_elts.add(elt)
#             res.append(col_elts)
#         return res


class Subreads:
    """Set of subreads that can be clustered

    A set of Subread objects for which we remember the original read
    order and can cluster to generate candidate haplotypes.

    Attributes:
        sreads : dict subread name -> Subread object
        mod_adj : record adjacency between subreads

    Methods:
        splitReads : import full reads and make Subreads
        profileSubreads : prepare the profile for all subreads
        clusterSubreads : cluster the subreads into N groups
        threadDiplotype : enumerate pairs of haplotypes
    """
    def __init__(self):
        self.sreads = {}
        self.mod_adj = {}
        self.pos_alleles = {}
        self.read_prof = {}
        self.read_to_sub = {}

    def print(self, long=False):
        read_to_sub = {}
        nparts = {'flankr': 0, 'flankl': 0, 'buffer': 0, 'module': 0}
        ncls = {}
        for sreadn in self.sreads:
            cl = self.sreads[sreadn].cluster
            read_name = self.sreads[sreadn].read_name
            if read_name not in read_to_sub:
                read_to_sub[read_name] = []
            read_to_sub[read_name].append(sreadn)
            if self.sreads[sreadn].in_module:
                nparts['module'] += 1
                if cl not in ncls:
                    ncls[cl] = 0
                ncls[cl] += 1
            else:
                nparts[cl] += 1
        parts_s = []
        for pp in nparts:
            if pp == 'module':
                cls_s = ['{}:{}'.format(cl, ncls[cl]) for cl in ncls]
                cls_s = ', '.join(cls_s)
                parts_s.append('{} {} ({})'.format(nparts[pp], pp, cls_s))
            else:
                parts_s.append('{} {}'.format(nparts[pp], pp))
        parts_s = ', '.join(parts_s)
        print('Subreads with {} reads split in {} parts'
              ': {}.'.format(len(read_to_sub), len(self.sreads), parts_s))
        if long:
            for readn in read_to_sub:
                print("\t" + readn)
                sreads = read_to_sub[readn]
                sreads.sort()
                for sreadn in sreads:
                    sread = self.sreads[sreadn]
                    print('\t\t{}: {} ({})'.format(sreadn, sread.cluster,
                                                   sread.cluster_conf))

    def clearClusters(self):
        """Clear the cluster assignment for subreads"""
        for sreadn in self.sreads:
            self.sreads[sreadn].clearCluster()

    def splitReads(self, reads, nodes):
        # find cycle's boundaries and start node
        cyc_l = None
        cyc_r = None
        start_node = None
        for node in nodes:
            if nodes[node]['class'] == 'cyc_l':
                cyc_l = node
            elif nodes[node]['class'] == 'cyc_r':
                cyc_r = node
            elif nodes[node]['class'] == 'ref':
                if start_node is None:
                    start_node = node
                else:
                    if nodes[node]['rpos_min'] < nodes[start_node]['rpos_min']:
                        start_node = node
        # group the nodes into flankl, flankr buffer and module
        # follow the reference path
        node = start_node
        node_group = {}
        cur_group = 'flankl'
        buffer_done = False
        while node is not None:
            if node == cyc_l:
                if cur_group == 'buffer':
                    buffer_done = True
                cur_group = 'module'
            elif node == cyc_r:
                if buffer_done:
                    cur_group = 'flankr'
                else:
                    cur_group = 'buffer'
            else:
                node_group[node] = cur_group
            nnode = None
            nn_pos = None
            for nn in nodes[node]['sucs']:
                # just consider nodes on the reference
                if nodes[nn]['ref'] == 0:
                    continue
                # skip if the rpos_min is not increasing except when the
                # buffer meets the left cycle bound
                if (nodes[nn]['rpos_min'] <= nodes[node]['rpos_min'] and
                        (cur_group != 'buffer' or nn != cyc_l)):
                    continue
                # update the position information
                if nn_pos is None:
                    nn_pos = nodes[nn]['rpos_min']
                    nnode = nn
                elif nn_pos > nodes[nn]['rpos_min'] and not buffer_done:
                    # prioritize reference edges "closer" to previous position
                    # we want to loop once
                    nn_pos = nodes[nn]['rpos_min']
                    nnode = nn
                elif nn_pos < nodes[nn]['rpos_min'] and buffer_done:
                    # prioritize reference edges "farther" to previous position
                    # because we've looped once already
                    nn_pos = nodes[nn]['rpos_min']
                    nnode = nn
            node = nnode
        # process each read
        for readn in reads.path:
            # skip if read with no informative nodes
            any_inf_nodes = False
            for nod in reads.path[readn]:
                if nodes[nod]['class'] in ['c1', 'c2']:
                    any_inf_nodes = True
                    break
            if not any_inf_nodes:
                continue
            # split the reads at node involved in the cycle
            subreads = [[]]
            subreads_t = []
            for nod in reads.path[readn]:
                if len(subreads_t) < len(subreads) and nod in node_group:
                    subreads_t.append(node_group[nod])
                if nod == cyc_l:
                    if len(subreads[-1]) > 0:
                        subreads.append([])
                if nod == cyc_r:
                    subreads[-1].append(nod)
                    subreads.append([])
                else:
                    subreads[-1].append(nod)
            # skip if just one subread and not in the region of interest
            if len(subreads_t) == 1 and subreads_t[0] != 'module':
                continue
            # otherwise, save the subreads
            for sbi, sr_path in enumerate(subreads):
                sreadn = '{}_{}'.format(readn, sbi)
                subread = Subread(sreadn, sr_path, readn)
                if subreads_t[sbi] != 'module':
                    subread.in_module = False
                    subread.assignCluster(subreads_t[sbi], 1)
                self.sreads[sreadn] = subread
                # record adjacency with next read
                if readn not in self.read_to_sub:
                    self.read_to_sub[readn] = []
                self.read_to_sub[readn].append(sreadn)
                if len(subreads) > sbi + 1:
                    self.mod_adj[sreadn] = '{}_{}'.format(readn, sbi + 1)

    def profileSubreads(self, marker_finder, min_read_support=3):
        # prepare marker profiles of each subread
        self.read_prof = {}
        count_alleles = {}
        for sreadn in self.sreads:
            # skip if not within the module
            if not self.sreads[sreadn].in_module:
                continue
            obs = self.sreads[sreadn].makeObs(marker_finder)
            self.read_prof[sreadn] = obs
            for pos in obs.pos_val:
                if pos not in count_alleles:
                    count_alleles[pos] = {}
                if obs.pos_val[pos] not in count_alleles[pos]:
                    count_alleles[pos][obs.pos_val[pos]] = 0
                count_alleles[pos][obs.pos_val[pos]] += 1
        # find positions that are not supported by enough reads
        self.pos_alleles = {}
        pos_to_rm = []
        for pos in count_alleles:
            supp_als = 0
            for al in count_alleles[pos]:
                if count_alleles[pos][al] >= min_read_support:
                    supp_als += 1
            if supp_als > 1:
                if pos not in self.pos_alleles:
                    self.pos_alleles[pos] = set()
                for al in count_alleles[pos]:
                    if count_alleles[pos][al] >= min_read_support:
                        self.pos_alleles[pos].add(al)
            else:
                pos_to_rm.append(pos)
        # remove them from the profiles
        for sreadn in self.sreads:
            # skip if not within the module
            if not self.sreads[sreadn].in_module:
                continue
            self.read_prof[sreadn].removePositions(pos_to_rm)

    def clusterSubreads(self, n_clusters):
        # run EM to create profiles
        best_em = None
        # init the number of differences with a high number
        best_em_ndiff = None
        # try a couple of attempts and keep the best
        # should be enough to try ~25 attempts with ~40 iterations
        for attempt in range(20):
            # init the EM profiles with 10 states
            em = EM(self.pos_alleles, n_clusters)
            # iterate
            for iter in range(40):
                em.updateWithObs(self.read_prof, init=iter == 0)
            em.consensusWithObs(self.read_prof)
            ndiffs = em.computeReadDiff(self.read_prof)
            ndiff = ndiffs['ndiffs']
            if best_em_ndiff is None or ndiff < best_em_ndiff:
                best_em_ndiff = ndiff
                best_em = em
            attempt += 1
        # ndiffs = best_em.computeReadDiff(self.read_prof)
        # cn = best_em.estimateCopyNumber(ndiffs)
        # print('Copy number: {}'.format(cn))
        # for key in ndiffs:
        #     print('{}: {}'.format(key, ndiffs[key]))
        # assign clusters based on best profile
        for sreadn in self.read_prof:
            ob = self.read_prof[sreadn]
            # find the best profile
            ob.findBestEMProfile(best_em)
            self.sreads[sreadn].assignCluster(ob.profile, ob.profile_conf)

    def clusterSubreadsEMsupp(self, n_clusters):
        # run EM to create profiles
        best_em = None
        best_em_ndiff = None
        # try a couple of attempts and keep the best
        # should be enough to try ~25 attempts with ~40 iterations
        for attempt in range(50):
            # init the EM profiles with 10 states
            em = EM(self.pos_alleles, n_clusters)
            # iterate
            for iter in range(100):
                em.updateWithObs(self.read_prof, init=iter == 0)
            em.consensusWithObs(self.read_prof)
            ndiffs = em.computeReadDiff(self.read_prof)
            ndiff = ndiffs['ndiffs']
            if best_em_ndiff is None or ndiff < best_em_ndiff:
                best_em_ndiff = ndiff
                best_em = em
            attempt += 1
        # remove states that are not covered enough
        ndiffs = best_em.computeReadDiff(self.read_prof)
        best_em.keepSuppStates(ndiffs)
        # assign clusters based on best profile
        for sreadn in self.read_prof:
            ob = self.read_prof[sreadn]
            # find the best profile
            ob.findBestEMProfile(best_em)
            self.sreads[sreadn].assignCluster(ob.profile, ob.profile_conf)

    def clusterSubreadsEMtarget(self, n_clusters):
        # run EM to create profiles
        best_em = None
        # init the number of differences with a high number
        best_em_ndiff = None
        # try a couple of attempts and keep the best
        # should be enough to try ~25 attempts with ~40 iterations
        attempt = 0
        while attempt < 10:
            # init the EM profiles with 10 states
            em = EM(self.pos_alleles, 8)
            # iterate
            for iter in range(20):
                em.updateWithObs(self.read_prof, init=iter == 0)
            em.consensusWithObs(self.read_prof)
            ndiffs = em.computeReadDiff(self.read_prof)
            # record copy number estimate
            cn = em.estimateCopyNumber(ndiffs)
            if cn != n_clusters:
                continue
            # remove a state until the proprer number
            em.keepSuppStates(ndiffs)
            for iter in range(20):
                em.updateWithObs(self.read_prof)
            em.consensusWithObs(self.read_prof)
            ndiffs = em.computeReadDiff(self.read_prof)
            ndiff = ndiffs['ndiffs']
            if best_em_ndiff is None or ndiff < best_em_ndiff:
                best_em_ndiff = ndiff
                best_em = em
            attempt += 1
        ndiffs = best_em.computeReadDiff(self.read_prof)
        # assign clusters based on best profile
        for sreadn in self.read_prof:
            ob = self.read_prof[sreadn]
            # find the best profile
            ob.findBestEMProfile(best_em)
            self.sreads[sreadn].assignCluster(ob.profile, ob.profile_conf)

    def clusterSubreadsGreedy(self, n_clusters):
        # start with one profile
        clusterer = GreedyMarkerBiClusterer(self.read_prof)
        # iteratively split the profile with most assigned reads
        while clusterer.nstates < n_clusters:
            clusterer.splitProfile(self.read_prof)
            # polish the profiles with the EM approach?
            for ii in range(10):
                clusterer.updateWithObs(self.read_prof)
            clusterer.consensusWithObs(self.read_prof)
        # assign clusters based on best profile
        for sreadn in self.read_prof:
            ob = self.read_prof[sreadn]
            # find the best profile
            ob.findBestEMProfile(clusterer)
            self.sreads[sreadn].assignCluster(ob.profile, ob.profile_conf)

    def clusterSubreadsPhaseBlockMerging(self, n_clusters):
        # init the phase blocks
        pos_to_pb = {}
        pblocks = []
        for pos in self.pos_alleles:
            if pos in pos_to_pb:
                continue
            pos_to_pb[pos] = len(pblocks)
            pb = PhaseBlock(n_clusters, pos, self.read_prof)
            pblocks.append(pb)
        if DEBUG:
            print('{} phase blocks'.format(len(pblocks)))
        # define an order to combine the blocks
        # compute the coverage for each position pair
        pair_cov = {}
        tot_cov = 0
        for sr in self.read_prof:
            ob = self.read_prof[sr]
            for p1 in ob.pos_val:
                for p2 in ob.pos_val:
                    if p1 <= p2:
                        continue
                    if (p1, p2) not in pair_cov:
                        pair_cov[(p1, p2)] = 0
                    pair_cov[(p1, p2)] += 1
                    tot_cov += 1
        # start with mean coverage? highest coverage might be enriched
        # in small/bad reads
        cov_m = tot_cov / len(pair_cov)
        pairs_sorted = sorted(pair_cov, key=lambda k: abs(cov_m - pair_cov[k]),
                              reverse=True)
        # pairs_sorted = sorted(pair_cov, key=lambda k: pair_cov[k],
        #                       reverse=True)
        # random.shuffle(pairs_sorted)
        if DEBUG:
            print('Top PB pairs:')
            ii = 0
            for p1, p2 in pairs_sorted:
                print('{} - {}'.format(p1, p2))
                ii += 1
                if ii > 10:
                    break
            print('{} phase blocks'.format(len(pblocks)))
        # combine the phase blocks
        for p1, p2 in pairs_sorted:
            idx1 = pos_to_pb[p1]
            idx2 = pos_to_pb[p2]
            if idx1 == idx2:
                # already combined these positions
                continue
            # otherwise combine the second with the first
            pblocks[idx1].combinePhaseBlock(pblocks[idx2], self.read_prof)
            for pos in pblocks[idx2].pbprof[0]:
                pos_to_pb[pos] = idx1
            # adjust profile probabilities
            pblocks[idx1].adjustProfileProbs(self.read_prof)
        # export final profile
        pb = pos_to_pb[list(pos_to_pb.keys())[0]]
        pb = pblocks[pb]
        # pb.makeProfile()
        # TODO check that all blocks are merged?
        # assign clusters based on best profile
        for sreadn in self.read_prof:
            ob = self.read_prof[sreadn]
            # find the best profile
            mprob = pb.computeReadProb(ob)
            if DEBUG:
                print('{} probs: {}'.format(sreadn, mprob))
            best_prob = None
            best_prof = None
            for ii in range(len(mprob)):
                if best_prob is None or best_prob < mprob[ii]:
                    best_prob = mprob[ii]
                    best_prof = ii
            ob.profile = best_prof
            mprob.sort(reverse=True)
            ob.profile_conf = mprob[0] - mprob[1]
            self.sreads[sreadn].assignCluster(ob.profile, ob.profile_conf)

    def makeConsensus(self, nodes):
        """Make consensus for each cluster"""
        # list clusters
        cls = set()
        for sreadn in self.sreads:
            # skip if not within the module
            if not self.sreads[sreadn].in_module:
                continue
            cls.add(self.sreads[sreadn].cluster)
        # find the first cycle boundary
        start_node = None
        for node in nodes:
            if nodes[node]['class'] == 'cyc_l':
                start_node = node
                break
        if start_node is None:
            print("Error: couldn't find cycle boundary.")
        # find consensus for each consensus
        warn_low_cov_nodes = set()
        cons_path = {}
        for cl in cls:
            # index the edge coverage for subread in the cluster
            eidx = EdgeIndex()
            for sreadn in self.sreads:
                if self.sreads[sreadn].cluster == cl:
                    eidx.addPath(self.sreads[sreadn].path)
            # eidx.print()
            # start from the first cycling node
            path = [start_node]
            while nodes[path[-1]]['class'] != 'cyc_r':
                node = path[-1]
                nnode = eidx.getBestSuccessor(node)
                if nnode is None:
                    warn_low_cov_nodes.add(node)
                    # pick most frequent successor in the graph
                    freq_suc = -1
                    for nn in nodes[node]['sucs']:
                        freq = nodes[nn]['c1'] + nodes[nn]['c2']
                        if freq > freq_suc:
                            nnode = nn
                    if nnode is None:
                        print("Error: no successor during consensus")
                        exit
                path.append(nnode)
            cons_path[cl] = path
        if len(warn_low_cov_nodes) > 0:
            print("Warning: no coverage during consensus "
                  "(for {} nodes).".format(len(warn_low_cov_nodes)))
        return cons_path

    def threadDiplotype(self, nodes):
        """Thread clusters of subreads into a pair of haplotypes

        Use the module adjacency (how reads were split into subreads)
        to connect clusters. Find the two best paths/haplotypes and
        return a list of their path.

        Args:
            nodes : dict with node information

        Returns: a dict: haplotype name -> path (list of nodes)

        """
        # make consensus for each cluster
        cl_paths = self.makeConsensus(nodes)
        # prepare the consensus for the flanks and buffer?
        fl_l_path = makeFlankConsensus(nodes, 'flankl')
        fl_r_path = makeFlankConsensus(nodes, 'flankr')
        buff_path = makeFlankConsensus(nodes, 'buffer')
        # create a module diplotyper
        md = ModuleDiplotyper(include_buffer=len(buff_path) > 0)
        for readn in self.read_to_sub:
            md.addRead(self.read_to_sub[readn], self.sreads)
        md.prepareDiplotypes()
        if DEBUG:
            md.print()
        best_dip = md.findBestDiplotype()
        # prepare full haplotype paths
        res = {}
        for hap in best_dip:
            hap = best_dip[hap]
            hapn = '_'.join([str(hh) for hh in hap])
            path = []
            for mod in hap:
                if mod == 'flankl':
                    path += fl_l_path
                elif mod == 'flankr':
                    path += fl_r_path
                elif mod == 'buffer':
                    path += buff_path
                else:
                    path += cl_paths[mod]
            res[hapn] = path
        return res


def makeFlankConsensus(nodes, flank_type='flankl'):
    """Prepare the path for a flanking region

    Uses the node information to enumerate the path of a flanking or
    buffer region. It basically starts at the first node or a cycle
    boundary and follows "successors" on the reference path.

    Args:
        nodes : a dict with node information
        flank_type : which flank: "flankl", "flankr", or "buffer"

    Returns: a path as a list of nodes
    """
    if flank_type not in ['flankl', 'flankr', 'buffer']:
        print('Error: wrong flank type in makeFlankConsensus.')
        exit
    # find cycle's boundaries and start node
    cyc_l = None
    cyc_r = None
    start_node = None
    for node in nodes:
        if nodes[node]['class'] == 'cyc_l':
            cyc_l = node
        elif nodes[node]['class'] == 'cyc_r':
            cyc_r = node
        elif nodes[node]['class'] == 'ref':
            if start_node is None:
                start_node = node
            else:
                if nodes[node]['rpos_min'] < nodes[start_node]['rpos_min']:
                    start_node = node
    # start the path and define end node
    path = []
    if flank_type == 'flankl':
        path.append(start_node)
    else:
        path.append(cyc_r)
    # follow the reference path
    while True:
        # stop if we're done
        if flank_type in ['flankl', 'buffer'] and path[-1] == cyc_l:
            break
        elif len(nodes[path[-1]]['sucs']) == 0:
            break
        # find next node
        node = path[-1]
        nn_pos = None
        nnode = None
        for nn in nodes[node]['sucs']:
            # just consider nodes on the reference
            if nodes[nn]['ref'] == 0:
                continue
            # skip if the rpos_min is not increasing except when the
            # buffer meets the left cycle bound
            if (nodes[nn]['rpos_min'] <= nodes[node]['rpos_min'] and
                    (flank_type != 'buffer' or nn != cyc_l)):
                continue
            # update the position information
            if nn_pos is None:
                nn_pos = nodes[nn]['rpos_min']
                nnode = nn
            elif (flank_type == 'buffer' and
                  (nn_pos > nodes[nn]['rpos_min'] or nn == cyc_l)):
                # for the buffer we want the most upstream ref node
                nn_pos = nodes[nn]['rpos_min']
                nnode = nn
            elif flank_type == 'flankr' and nn_pos < nodes[nn]['rpos_min']:
                # for the right flank we want the most downstream ref node
                nn_pos = nodes[nn]['rpos_min']
                nnode = nn
        path.append(nnode)
    # remove cycle's boundary
    if path[0] in [cyc_l, cyc_r]:
        path = path[1:]
    if path[-1] in [cyc_l, cyc_r]:
        path = path[:-1]
    return path


class EdgeIndex:
    """Simple class to tally edge counts

    Helpher class to manage a dictionnary edge -> count.

    Attributes:
        edge_count : the dict edge -> count

    Methods:
        addPath : increment the counts from a path (node list)
        getBestSuccessor : return most supported next node
    """

    def __init__(self):
        self.edge_count = {}

    def addPath(self, path):
        for pos in range(len(path)-1):
            if path[pos] not in self.edge_count:
                self.edge_count[path[pos]] = {}
            if path[pos+1] not in self.edge_count[path[pos]]:
                self.edge_count[path[pos]][path[pos+1]] = 0
            self.edge_count[path[pos]][path[pos+1]] += 1

    def getBestSuccessor(self, node):
        if node not in self.edge_count:
            return None
        best_suc = None
        best_suc_count = 0
        for suc in self.edge_count[node]:
            if self.edge_count[node][suc] > best_suc_count:
                best_suc = suc
                best_suc_count = self.edge_count[node][suc]
        return best_suc

    def print(self, long=False):
        min_count = None
        max_count = None
        for e1 in self.edge_count:
            for e2 in self.edge_count[e1]:
                if min_count is None or min_count > self.edge_count[e1][e2]:
                    min_count = self.edge_count[e1][e2]
                if max_count is None or max_count < self.edge_count[e1][e2]:
                    max_count = self.edge_count[e1][e2]
        print('EdgeIndex with {} edges with counts in '
              '[{}, {}].'.format(len(self.edge_count), min_count, max_count))
        if long:
            for edge in self.edge_count:
                print('\t{}: {}'.format(edge, self.edge_count[edge]))


class Observation:
    """A series of position and markers observed in a read

    Represent a read or subread and the corresponding
    marker/allele. This Observation can be used to assign the read to
    a profile (also made up of positions and allele/values).

    Attributes:
        pos_val : dict position -> value (~allele)
        profile : the assigned profile
        profile_conf : our confidence in the assigned profile

    Methods:
        findBestEMProfile : which profile from an EM best match the Observation
        removePositions : remove a list of positions from the Observation

    """

    def __init__(self, position_values):
        self.pos_val = position_values
        self.profile = None
        self.profile_conf = 0
        self.profile_mismatch = 0

    def print(self, long=True):
        if long:
            for pos in self.pos_val:
                print('\t{}: {}'.format(pos, self.pos_val[pos]))
        else:
            out_h = ''
            out_v = ''
            prev_size = 1
            for pos in self.pos_val:
                out_h += ' ' * (4 - prev_size) + str(pos)
                if self.pos_val[pos]:
                    out_v += ' ' * (4 - 1) + '+'
                else:
                    out_v += ' ' * (4 - 1) + '-'
                prev_size = len(str(pos))
            print(out_h)
            print(out_v)

    def findBestEMProfile(self, em, include_none=False):
        best_state = []
        best_count = -1
        best_mm = 0
        all_counts = []
        for state in range(em.nstates):
            match_count = 0
            mismatch_count = 0
            for pos in self.pos_val:
                if em.prof[pos][state] == self.pos_val[pos]:
                    match_count += 1
                elif include_none and em.prof[pos][state] is None:
                    match_count += 1
                else:
                    mismatch_count += 1
            all_counts.append(match_count)
            if match_count > best_count:
                best_count = match_count
                best_mm = mismatch_count
                best_state = [state]
            elif match_count == best_count:
                best_state.append(state)
        # assign to that profile
        if len(best_state) > 1:
            random.shuffle(best_state)
        self.profile = best_state[0]
        self.profile_mismatch = best_mm
        # compute a confidence metric as the difference between best
        # profile and second best profile
        if len(all_counts) < 2:
            self.profile_conf = 1
        else:
            all_counts.sort(reverse=True)
            self.profile_conf = all_counts[0] - all_counts[1]

    def removePositions(self, pos_to_rm):
        """Remove specified positions from the profile"""
        for pos in pos_to_rm:
            if pos in self.pos_val:
                del (self.pos_val[pos])

    def size(self):
        return len(self.pos_val)


class EM:
    """EM profile builder

    Builde N profiles using an EM approach: iteratively assign subread
    to best profile and update profile.

    Attributes:
        prof : the profiles as a dict of list. Use as prof[position][state]
        nstates : number of profiles to consider
        marker_values : a dict pos -> set of the values each marker can take

    Methods:
        constructor : init the profiles with random values
        updateWithObs : run a round of direct profile updates from Observations
        consensusWithObs : assign Observations and update consensus profile
        compareWithTruth : helper function to benchmark with a truth profile

    """

    def __init__(self, marker_values, nstates_per_pos):
        self.prof = {}
        self.nstates = nstates_per_pos
        self.marker_values = marker_values
        for pos in marker_values:
            vals = []
            vals_space = list(self.marker_values[pos])
            for state in range(self.nstates):
                # init with a random value
                vals.append(random.sample(vals_space, 1)[0])
            self.prof[pos] = vals

    def print(self, digits=4):
        for pos in self.marker_values:
            tp = []
            for val in self.prof[pos]:
                tp.append('{}'.format(val))
            print('S{}:\t'.format(pos) + '  '.join(tp))

    def updateWithObs(self, obs, long_reads_first=False, init=False):
        obs_names = list(obs.keys())
        if long_reads_first:
            obs_names = sorted(list(obs.keys()), key=lambda k: -obs[k].size())
        else:
            random.shuffle(obs_names)
        # loop over observations and integrate them one by one
        for obsn in obs_names:
            ob = obs[obsn]
            # find the best profile
            ob.findBestEMProfile(self, include_none=init)
            # overwrite profile
            for pos in ob.pos_val:
                self.prof[pos][ob.profile] = ob.pos_val[pos]

    def consensusWithObs(self, obs):
        # loop over observations, assign, and tally
        counts = {}
        for obn in obs:
            ob = obs[obn]
            # find the best profile
            ob.findBestEMProfile(self)
            # count each type of marker at each position of this profile
            for pos in ob.pos_val:
                oname = '{}_{}'.format(pos, ob.profile)
                if oname not in counts:
                    counts[oname] = {}
                if ob.pos_val[pos] not in counts[oname]:
                    counts[oname][ob.pos_val[pos]] = 0
                counts[oname][ob.pos_val[pos]] += 1
        # pick the most supported value for each profile position
        for pp in self.prof:
            for state in range(self.nstates):
                top_val = None
                top_val_c = 0
                oname = '{}_{}'.format(pp, state)
                if oname not in counts:
                    continue
                for val in counts[oname]:
                    if top_val is None or counts[oname][val] > top_val_c:
                        top_val_c = counts[oname][val]
                        top_val = val
                if top_val is not None:
                    self.prof[pp][state] = top_val

    def computeReadDiff(self, obs):
        # loop over observations, assign, and tally
        ndiffs = 0
        cl_ndiffs = {}
        supp = {}
        for prof in range(self.nstates):
            cl_ndiffs[prof] = {'ndiffs': 0, 'nmatch': 0, 'conf': 0}
            supp[prof] = {}
        for obn in obs:
            ob = obs[obn]
            # find the best profile
            ob.findBestEMProfile(self)
            # if ob.profile_conf == 0:
            #     continue
            cl_ndiffs[ob.profile]['conf'] += ob.profile_conf
            # count each type of marker at each position of this profile
            for pos in ob.pos_val:
                if ob.pos_val[pos] != self.prof[pos][ob.profile]:
                    ndiffs += 1
                    cl_ndiffs[ob.profile]['ndiffs'] += 1
                else:
                    if pos not in supp[ob.profile]:
                        supp[ob.profile][pos] = 0
                    supp[ob.profile][pos] += 1
                    cl_ndiffs[ob.profile]['nmatch'] += 1
        # compute support for each position in the profile
        for state in range(self.nstates):
            supp_v = []
            supp_c = 0
            for pos in self.prof:
                if pos in supp[state]:
                    supp_v.append(supp[state][pos])
                    supp_c += 1
                else:
                    supp_v.append(0)
            cl_ndiffs[state]['nsupp'] = supp_c
            cl_ndiffs[state]['nsupp_mean'] = sum(supp_v) / len(supp_v)
            ntot = cl_ndiffs[state]['ndiffs'] + cl_ndiffs[state]['nmatch']
            if ntot > 0:
                mprop = round(cl_ndiffs[state]['nmatch'] / ntot, 4)
                cl_ndiffs[state]['match_prop'] = mprop
            else:
                cl_ndiffs[state]['match_prop'] = 0
        cl_ndiffs['ndiffs'] = ndiffs
        return (cl_ndiffs)

    def removeOneState(self, ndiffs):
        # find the state with least support
        supp_min = None
        supp2_min = None
        st_to_rm = None
        for state in range(self.nstates):
            if (supp_min is None or supp_min > ndiffs[state]['nsupp']
                    or (supp_min == ndiffs[state]['nsupp']
                        and supp2_min > ndiffs[state]['nmatch'])):
                supp_min = ndiffs[state]['nsupp']
                supp2_min = ndiffs[state]['nmatch']
                st_to_rm = state
        # make a new profile without this state
        prof = {}
        for pos in self.prof:
            prof[pos] = []
            for state in range(self.nstates):
                if state == st_to_rm:
                    continue
                prof[pos].append(self.prof[pos][state])
        # update state number
        self.nstates += -1

    def keepSuppStates(self, ndiffs):
        # find the state with least support
        max_supp = None
        for state in range(self.nstates):
            if max_supp is None or max_supp < ndiffs[state]['nsupp']:
                max_supp = ndiffs[state]['nsupp']
        # how many states with about the maximum support
        st_to_keep = []
        for state in range(self.nstates):
            if 0.9 * max_supp < ndiffs[state]['nsupp']:
                st_to_keep.append(state)
        # make a new profile without this state
        prof = {}
        for pos in self.prof:
            prof[pos] = []
            for state in range(self.nstates):
                if state == st_to_keep:
                    prof[pos].append(self.prof[pos][state])
        # update state number
        self.nstates = len(st_to_keep)

    def estimateCopyNumber(self, ndiffs):
        # supps = []
        # find the state with least support
        max_supp = None
        for state in range(self.nstates):
            if max_supp is None or max_supp < ndiffs[state]['nsupp']:
                max_supp = ndiffs[state]['nsupp']
            # supps.append(ndiffs[state]['nsupp'])
        # supps.sort()
        # print(supps)
        # how many states with about the maximum support
        nsupp = 0
        for state in range(self.nstates):
            if 0.9 * max_supp < ndiffs[state]['nsupp']:
                nsupp += 1
        return nsupp

    def compareWithTruth(self, hap_truth):
        # compare each profile to each true haplotype
        best_truth = [-1] * self.nstates
        best_truth_match = [0] * self.nstates
        for prof_ii in range(self.nstates):
            for truth_ii in range(len(hap_truth)):
                nmatch = 0
                for pos in self.marker_values:
                    if self.prof[pos][prof_ii] == hap_truth[truth_ii][pos]:
                        nmatch += 1
                if nmatch > best_truth_match[prof_ii]:
                    best_truth_match[prof_ii] = nmatch
                    best_truth[prof_ii] = truth_ii
        # any true haplotype not found
        missing_truth = []
        for truth_ii in range(len(hap_truth)):
            if truth_ii not in best_truth:
                missing_truth.append(truth_ii)
        # print profile next to their most similar true haplotype
        tp = 'read\t'
        for prof_ii in range(self.nstates):
            tp += '  {}/{}'.format(prof_ii, best_truth[prof_ii])
        tp += '\t'
        for truth_ii in missing_truth:
            tp += '  {}'.format(truth_ii)
        print(tp)
        print('---')
        for pp in self.marker_values:
            tp = str(pp) + '\t'
            for p_ii in range(self.nstates):
                match_s = 'x'
                if self.prof[pp][p_ii] == hap_truth[best_truth[p_ii]][pp]:
                    match_s = '='
                tp += '  {}{}{}'.format(self.prof[pp][p_ii], match_s,
                                        hap_truth[best_truth[p_ii]][pp])
            tp += '\t'
            for truth_ii in missing_truth:
                tp += '  {}'.format(hap_truth[truth_ii][pp])
            print(tp)


class Segment:
    """A phased segment object used to split profiles

    A segment is just a pair of alleles for a set of positions. It's
    used in the profile splitter where we start with single-position
    segments and iteratively combine pairs based on the "phasing"
    information from the reads.

    Attributes:
        al1 : dict pos->allele for each position for the first haplotype
        al2 : dict pos->allele for each position for the second haplotype

    Methods:
        constructor : init with one position and two alleles
        combineSegment : combine another Segment with this one
    """

    def __init__(self, pos, al1, al2):
        self.al1 = {}
        self.al1[pos] = al1
        self.al2 = {}
        self.al2[pos] = al2

    def print(self):
        pos_sum = []
        for pos in self.al1:
            if len(pos_sum) > 3:
                break
            pos_sum.append(pos)
        pos_sum = ', '.join(pos_sum)
        print('Segment with {} markers ({}...).'.format(len(self.al1),
                                                        pos_sum))

    def combineSegment(self, seg_o, read_prof, read_list=[]):
        # find out if we should combine al1-al1/al2-al2 (keep phase)
        # or al1-al2/al2-al1 (flip phase)
        flip_phase_supp = 0
        keep_phase_supp = 0
        for sr in read_prof:
            # skip if not a read we want
            if len(read_list) > 0 and sr not in read_list:
                continue
            ob = read_prof[sr]
            # which allele is better support in this segment?
            seg_al1 = None
            a1_supp = 0
            a2_supp = 0
            for pos in ob.pos_val:
                if pos in self.al1 and ob.pos_val[pos] == self.al1[pos]:
                    a1_supp += 1
                if pos in self.al2 and ob.pos_val[pos] == self.al2[pos]:
                    a2_supp += 1
            if a1_supp == a2_supp:
                # skip non-informative reads
                continue
            elif a1_supp > a2_supp:
                seg_al1 = True
            else:
                seg_al1 = False
            # which allele is better support in this segment?
            seg_al2 = None
            a1_supp = 0
            a2_supp = 0
            for pos in ob.pos_val:
                if pos in seg_o.al1 and ob.pos_val[pos] == seg_o.al1[pos]:
                    a1_supp += 1
                if pos in seg_o.al2 and ob.pos_val[pos] == seg_o.al2[pos]:
                    a2_supp += 1
            if a1_supp == a2_supp:
                # skip non-informative reads
                continue
            elif a1_supp > a2_supp:
                seg_al2 = True
            else:
                seg_al2 = False
            # tally the support for keeping/flipping the phases
            if seg_al1 == seg_al2:
                keep_phase_supp += 1
            else:
                flip_phase_supp += 1
        # now we know how to combine the alleles
        for pos in seg_o.al1:
            if keep_phase_supp > flip_phase_supp:
                self.al1[pos] = seg_o.al1[pos]
                self.al2[pos] = seg_o.al2[pos]
            else:
                self.al1[pos] = seg_o.al2[pos]
                self.al2[pos] = seg_o.al1[pos]


class ProfileSplitter:
    """Runs the profile-splitting approach for a set of reads

    Uses the profile from a list of reads to init a set of
    Segments. The Segments are then combined iteratively, starting
    with the most supported pairs of segments first.

    Attributes:
        mark_co_cov : coverage of a pair of markers
        pos_to_seg : dict position to Segment currently holding that position
        segs : list of segments

    Methods:
        constructor : build the coverage matrix and init the segments
        split : iteratively combine the segments
        getFinalSegment : return the remaining Segment with the phased profiles
    """

    def __init__(self, read_prof, read_list=[]):
        # how many reads cover two positions: mark_co_cov["edge1|edge2"]
        self.mark_co_cov = {}
        # start by counting all allele pairs across reads
        al_p_cov = {}
        # also count the allele support for each allele at a pos
        al_cov = {}
        # check each read
        for sr in read_prof:
            # skip if not a read we want
            if len(read_list) > 0 and sr not in read_list:
                continue
            # otherwise get the Observation profile
            obs = read_prof[sr]
            # record coverage for every position pairs
            for p1 in obs.pos_val:
                # update the allele coverage
                if p1 not in al_cov:
                    al_cov[p1] = {}
                if obs.pos_val[p1] not in al_cov[p1]:
                    al_cov[p1][obs.pos_val[p1]] = 0
                al_cov[p1][obs.pos_val[p1]] += 1
                # check for other positions to update allele pair coverage
                for p2 in obs.pos_val:
                    if p1 >= p2:
                        # record pairs once, as (p, P) with p the
                        # "smallest" position
                        continue
                    # increment the coverage dict
                    if p1 not in al_p_cov:
                        al_p_cov[p1] = {}
                    if p2 not in al_p_cov[p1]:
                        al_p_cov[p1][p2] = {}
                    al_pair = '{}+{}'.format(obs.pos_val[p1],
                                             obs.pos_val[p2])
                    if al_pair not in al_p_cov[p1][p2]:
                        al_p_cov[p1][p2][al_pair] = 0
                    al_p_cov[p1][p2][al_pair] += 1
        # keep the coverage of the second most supported allele pair
        for p1 in al_p_cov:
            for p2 in al_p_cov[p1]:
                cov = []
                for al_pair in al_p_cov[p1][p2]:
                    cov.append(al_p_cov[p1][p2][al_pair])
                if len(cov) > 1:
                    cov.sort(reverse=True)
                    pos_pair = '{}|{}'.format(p1, p2)
                    self.mark_co_cov[pos_pair] = cov[1]
        # create segments with top 2 most supported allele for each position
        self.segs = []
        self.pos_to_seg = {}
        # start with one segment per position
        for pos in al_cov:
            # skip if only one allele
            if len(al_cov[pos]) < 2:
                continue
            # sort the alleles by coverage
            al_sorted = sorted(al_cov[pos], key=lambda k: -al_cov[pos][k])
            # skip is low support for the second allele
            if al_cov[pos][al_sorted[1]] < 2:
                continue
            seg = Segment(pos, al_sorted[0], al_sorted[1])
            self.pos_to_seg[pos] = len(self.segs)
            self.segs.append(seg)

    def split(self, read_prof, read_list=[]):
        # find the next most supported position pair to merge
        pairs_to_merge = sorted(self.mark_co_cov,
                                key=lambda k: -self.mark_co_cov[k])
        for pos_pair in pairs_to_merge:
            # try to combine those two positions
            pos_pair_v = pos_pair.split('|')
            p1 = pos_pair_v[0]
            p2 = pos_pair_v[1]
            # skip if they are already in the same segment
            if (p1 not in self.pos_to_seg or p2 not in self.pos_to_seg or
                    self.pos_to_seg[p1] == self.pos_to_seg[p2]):
                continue
            # otherwise combine segments
            seg2 = self.segs[self.pos_to_seg[p2]]
            self.segs[self.pos_to_seg[p1]].combineSegment(seg2, read_prof,
                                                          read_list=read_list)
            # update the segments for the combined positions
            for p_seg2 in seg2.al1:
                # positions from second segment now point to first segment
                self.pos_to_seg[p_seg2] = self.pos_to_seg[p1]

    def getFinalSegment(self):
        # we should now have only one segment
        seg_idx = None
        for pos in self.pos_to_seg:
            if seg_idx is None:
                seg_idx = self.pos_to_seg[pos]
            elif seg_idx != self.pos_to_seg[pos]:
                print("Warning: ProfileSplitter failed to "
                      "combine all segments")
        # return segment
        return self.segs[seg_idx]


class GreedyMarkerBiClusterer:
    """Runs the iterative profile-splitting approach to cluster reads

    Starts with a single profile. Then iteratively split the profile
    with the highest read coverage. The profile will be split using
    the helper class ProfileSplitter.

    Attributes:
        nstates : the number of profiles currently
        prof : the profiles as a dict of list. Use as prof[position][state]

    Methods:
        constructor : init the profile
        splitProfile : assign reads to profiles and split the most covered one
    """

    def __init__(self, read_prof):
        # init the one profile/state
        self.nstates = 1
        self.prof = {}
        for sr in read_prof:
            for pos in read_prof[sr].pos_val:
                if pos not in self.prof:
                    # add this position to the single profile
                    self.prof[pos] = [read_prof[sr].pos_val[pos]]
                    # note: the init value shouldn't matter here

    def print(self):
        for pos in self.prof:
            tp = []
            for val in self.prof[pos]:
                tp.append('{}'.format(val))
            print('S{}:\t'.format(pos) + '  '.join(tp))

    def splitProfile(self, read_prof):
        # assign reads to each profile
        prof_to_reads = []
        prof_mismatch_c = [0] * self.nstates
        for prof in range(self.nstates):
            prof_to_reads.append([])
        for sr in read_prof:
            read_prof[sr].findBestEMProfile(self)
            prof = read_prof[sr].profile
            prof_to_reads[prof].append(sr)
            prof_mismatch_c[prof] += read_prof[sr].profile_mismatch
        # find profile with most reads/bases
        big_prof = None
        # big_cov = 0
        big_mm = 0
        for prof, reads in enumerate(prof_to_reads):
            if prof_mismatch_c[prof] > big_mm:
                big_mm = prof_mismatch_c[prof]
                big_prof = prof
        # make two profiles for those reads
        ps = ProfileSplitter(read_prof, prof_to_reads[big_prof])
        ps.split(read_prof, prof_to_reads[big_prof])
        # prepare a new profile, init with the profile that was split
        for pos in self.prof:
            self.prof[pos].append(self.prof[pos][big_prof])
        # replace profile and add new one
        seg = ps.getFinalSegment()
        for pos in seg.al1:
            self.prof[pos][big_prof] = seg.al1[pos]
            self.prof[pos][self.nstates] = seg.al2[pos]
        self.nstates += 1

    def consensusWithObs(self, obs):
        # loop over observations, assign, and tally
        counts = {}
        for obn in obs:
            ob = obs[obn]
            # find the best profile
            ob.findBestEMProfile(self)
            if ob.profile_conf == 0:
                continue
            # count each type of marker at each position of this profile
            for pos in ob.pos_val:
                oname = '{}_{}'.format(pos, ob.profile)
                if oname not in counts:
                    counts[oname] = {}
                if ob.pos_val[pos] not in counts[oname]:
                    counts[oname][ob.pos_val[pos]] = 0
                counts[oname][ob.pos_val[pos]] += 1
        # pick the most supported value for each profile position
        for pp in self.prof:
            for state in range(self.nstates):
                top_val = None
                top_val_c = 0
                oname = '{}_{}'.format(pp, state)
                if oname not in counts:
                    continue
                for val in counts[oname]:
                    if top_val is None or counts[oname][val] > top_val_c:
                        top_val_c = counts[oname][val]
                        top_val = val
                if top_val is not None:
                    self.prof[pp][state] = top_val

    def updateWithObs(self, obs):
        obs_names = sorted(list(obs.keys()), key=lambda k: -obs[k].size())
        # loop over observations and integrate them one by one
        for obsn in obs_names:
            ob = obs[obsn]
            # find the best profile
            ob.findBestEMProfile(self)
            if ob.profile_conf == 0:
                continue
            # overwrite profile
            for pos in ob.pos_val:
                self.prof[pos][ob.profile] = ob.pos_val[pos]


class PhaseBlock:
    def __init__(self, ploidy, position, read_prof):
        self.ploidy = ploidy
        # precompute the allele counts
        al_c = {}
        tot_ac = 0
        for sr in read_prof:
            if position in read_prof[sr].pos_val:
                al = read_prof[sr].pos_val[position]
                if al not in al_c:
                    al_c[al] = 0
                al_c[al] += 1
                tot_ac += 1
        # expected read coverage for one copy
        exp_cov = tot_ac / ploidy
        # init the P profiles
        self.pbprof = []
        for pp in range(ploidy):
            prof = {}
            prof[position] = {}
            exp_cov_prof = exp_cov
            for al in al_c:
                if exp_cov_prof <= 0:
                    # there is no more "coverage" left to distribute
                    prof[position][al] = 0
                    continue
                if al_c[al] == 0:
                    # "coverage" already distributed to this allele
                    prof[position][al] = 0
                elif al_c[al] > exp_cov_prof:
                    # distribute some coverage to this allele
                    prof[position][al] = exp_cov_prof / exp_cov
                    # this allele can still receive some coverage in
                    # another profile
                    al_c[al] += -round(exp_cov_prof)
                    # but there is no coverage left to distribute in
                    # this profile
                    exp_cov_prof = 0
                else:
                    # distribute all possible coverage to this allele
                    prof[position][al] = al_c[al] / exp_cov
                    # update how much coverage is left for this profile
                    exp_cov_prof += -al_c[al]
                    # we've used up this allele
                    al_c[al] = 0
            self.pbprof.append(prof)

    def combinePhaseBlock(self, o_block, read_prof):
        # compute match probability for each pair of profiles between
        # the two phase blocks add position to this block
        prof_pair_cost = [[0] * self.ploidy for p in range(self.ploidy)]
        for sr in read_prof:
            # compute probability of matching each profile in this block
            mprob1 = self.computeReadProb(read_prof[sr])
            # compute probability of matching each profile in this other block
            mprob2 = o_block.computeReadProb(read_prof[sr])
            for i1 in range(self.ploidy):
                for i2 in range(self.ploidy):
                    # define a cost for this read and profile pair
                    # weight by the read size/information content
                    prof_pair_cost[i1][i2] += (1 - mprob1[i1] * mprob2[i2]) * len(read_prof[sr].pos_val)
                    # prof_pair_cost[i1][i2] += 1 - mprob1[i1] * mprob2[i2]
        # find the strongest match between profiles
        # best_m = matchBipartite(prof_pair_cost)
        best_m = matchBipartiteOpt(prof_pair_cost)
        # combine other profile to this phase block
        for pp in range(self.ploidy):
            # other profile to combine
            pp2 = best_m[pp]
            for pos in o_block.pbprof[pp2]:
                self.pbprof[pp][pos] = o_block.pbprof[pp2][pos]
        # clean/empty the other phase block? risky except if we make
        # sure to copy the other block

    def computeReadProb(self, ob):
        """Compute the probability of this observation matching each profile"""
        # compute first probability for each state
        probs = []
        tot_probs = 0
        for state in range(self.ploidy):
            # multiply probabilities at each position
            prob_state = 1
            for pos in ob.pos_val:
                # skip if not a position in this block
                if pos not in self.pbprof[state]:
                    continue
                # check for this allele at that position
                al = ob.pos_val[pos]
                if al not in self.pbprof[state][pos]:
                    # don't think this should happen.
                    print('Warning: allele {} not modeled at position {} '
                          'in PhaseBlock'.format(al, pos))
                prob_state = prob_state * self.pbprof[state][pos][al]
            probs.append(prob_state)
            tot_probs += prob_state
        # rescale probability to [0-1]
        probs_rs = []
        for prob in probs:
            if tot_probs > 0:
                probs_rs.append(prob/tot_probs)
            else:
                probs_rs.append(1/self.ploidy)
        return probs_rs

    def adjustProfileProbs(self, read_prof):
        """Match reads and adjust probabilities"""
        # init new profile with 0s
        pbprof = []
        pos_s = set()
        for state in range(self.ploidy):
            prof = {}
            for pos in self.pbprof[state]:
                prof[pos] = {}
                for al in self.pbprof[state][pos]:
                    prof[pos][al] = 0
                pos_s.add(pos)
            pbprof.append(prof)
        for sr in read_prof:
            # compute probability of matching each profile in this block
            mprob = self.computeReadProb(read_prof[sr])
            for pos in read_prof[sr].pos_val:
                if pos not in pos_s:
                    continue
                for state in range(self.ploidy):
                    al = read_prof[sr].pos_val[pos]
                    pbprof[state][pos][al] += mprob[state]
        # normalize to "probabilities" TODO double-check
        for state in range(self.ploidy):
            for pos in pbprof[state]:
                tot_al = 0
                for al in pbprof[state][pos]:
                    tot_al += pbprof[state][pos][al]
                if DEBUG and tot_al == 0:
                    print('Warning: no allele support in state {} ({} positions), position {}'.format(state, len(pbprof[state]), pos))
                    print('Previous profile {}'.format(self.pbprof[state][pos][al]))
                for al in pbprof[state][pos]:
                    if tot_al == 0:
                        # pbprof[state][pos][al] = self.pbprof[state][pos][al]
                        pbprof[state][pos][al] = 1 / len(pbprof[state][pos])
                        # TODO is that normal?
                    else:
                        pbprof[state][pos][al] = pbprof[state][pos][al] / tot_al
        self.pbprof = pbprof

    def makeProfile(self):
        self.prof = {}
        self.nstates = self.ploidy
        for pos in self.pbprof[0]:
            vals = []
            for state in range(self.ploidy):
                best_al = None
                best_pb = None
                for al in self.pbprof[state][pos]:
                    pb = self.pbprof[state][pos][al]
                    if best_al is None or best_pb < pb:
                        best_al = al
                        best_pb = pb
                vals.append(best_al)
            self.prof[pos] = vals


class MarkerFinder:
    """Finds edge marker in the pangenoe

    Marker are edges where at least another edge in the pangenome is
    mutually exclusive. For example, the two alleles of a SNP bubble
    are made of edge markers. Similarly a deletion edge would be a
    marker mutually exclusive with all edges in the other part of the
    bubble.

    Attributes:
        e_reach_mat : reachability "matrix" between all edges
                        dict edge -> list of boolean
        eid : map an edge to the index in the reachability boolean list
        prev_edge : map an edge to a set of previous edges in the graph
        next_edge : map an edge to a set of next edges in the graph
        last_edge : last edge (exiting the module region)
        mark_me: map an edge marker to a set of mutually exclusive edges

    Methods:
        constructor : Prepare the edge structures using node information
        findMarkers : Find marker edges and their mutually exclusive markers
        getMarkers : edge name to dict: position -> marker alleles
    """
    def __init__(self, nodes):
        """Prepare a MarkerFinder object

        Go through the node information and prepare the edge graph and
        structures.

        Args:
            nodes : dict with the node information
        """
        # the objects to prepare
        self.prev_edge = {}
        self.next_edge = {}
        self.e_reach_mat = {}
        self.last_edge = None
        self.eid = {}
        self.mark_me = {}
        # start at the left cycle node and enumerate edges until right cycle
        for start_node in nodes:
            if nodes[start_node]['class'] == 'cyc_l':
                break
        node_stack = {start_node}
        node_done = set()
        while len(node_stack) > 0:
            node = node_stack.pop()
            # for each edge starting at that node, look for next edges
            for node_2 in nodes[node]['sucs']:
                # save that edge in the reachability matrix
                e1 = '{}_{}'.format(node, node_2)
                self.e_reach_mat[e1] = []
                # this is the last edge we want to consider
                if nodes[node]['class'] == 'cyc_r':
                    self.last_edge = e1
                    # no need to consider other edges with that node
                    # or record the next edges
                    continue
                # enumerate next edges
                for node_3 in nodes[node_2]['sucs']:
                    e2 = '{}_{}'.format(node_2, node_3)
                    # record that edge 1 is just before edge 2
                    if e2 not in self.prev_edge:
                        self.prev_edge[e2] = set()
                    self.prev_edge[e2].add(e1)
                    # record that edge 2 is just after edge 1
                    if e1 not in self.next_edge:
                        self.next_edge[e1] = set()
                    self.next_edge[e1].add(e2)
                # add that next node to the stack if not done already
                if node_2 not in node_done:
                    node_stack.add(node_2)
            # remember that we've processed this node
            node_done.add(node)
        # fill the reachability matrix with zeros
        n_edges = len(self.e_reach_mat)
        for edge in self.e_reach_mat:
            self.eid[edge] = len(self.eid)
            self.e_reach_mat[edge] = [False] * n_edges

    def print(self):
        nedges = 0
        nmarks = 0
        for edge in self.mark_me:
            nmarks += len(self.mark_me[edge])
            nedges += 1
        print('MarkerFinder: {} edges with {} mutually exclusive '
              'markers'.format(nedges, nmarks))

    def computeReachability(self):
        """Compute the reachability of edges

        Performs a reverse breadth-first search of the graph (starting
        from last_edge). At each step, go forward and update the
        reachable edges. When we encounter an edge whose reachability
        had already been computed, update the edge of interest using
        an "OR" operation, and stop the forward search. This fills up
        the reachability matrix e_reach_mat, propagating information
        when all the edges reachable from one edge have been recorded
        already.
        """
        # start from last edge and do a reverse breadth-first search
        edge_stack = [self.last_edge]
        edge_done = set()
        stack_ii = 0
        while stack_ii < len(edge_stack):
            edge = edge_stack[stack_ii]
            if edge in edge_done:
                stack_ii += 1
                continue
            # look forward and update reachability?
            f_stack = {edge}
            while len(f_stack) > 0:
                r_edge = f_stack.pop()
                self.e_reach_mat[edge][self.eid[r_edge]] = True
                if r_edge not in edge_done:
                    # we haven't computed the reachability of this
                    # edge yet so we'll need to keep going
                    if r_edge in self.next_edge:
                        for n_edge in self.next_edge[r_edge]:
                            f_stack.add(n_edge)
                else:
                    # we've already computed that edge reachability
                    # just update the reachability with an OR operation
                    self.updateReach(edge, r_edge)
            # mark that edge as done
            edge_done.add(edge)
            # add the previous edges to the stack
            if edge in self.prev_edge:
                for p_edge in self.prev_edge[edge]:
                    if p_edge not in edge_done:
                        edge_stack.append(p_edge)
            # move to next edge in the stack
            stack_ii += 1

    def updateReach(self, e_up, e_done):
        """Helper function to update the reachability matrix

        Update reachability of e_up with everything reachable from
        e_done in the reachability matrix e_reach_mat. We do this when
        e_done is reachable from e_up and we've already found all the
        edges reachable from e_done.

        Args:
            e_up : edge to update in the matrix
            e_done : edge with complete reachability, reachable by e_up

        """
        for eid in range(len(self.e_reach_mat[e_up])):
            self.e_reach_mat[e_up][eid] = (self.e_reach_mat[e_up][eid]
                                           or self.e_reach_mat[e_done][eid])

    def findMarkers(self):
        """Prepare all markers and their mutually exclusive markers

        Go over the reachability matrix and check every edge pair to
        test if they are mutually exclusive. Save that information for
        edges with at least one mutually exclusive edge in the mark_me
        dict.

        """
        # compute the reachability of every edges
        self.computeReachability()
        # prepare the dict mapping edges to their mutually exclusive edges
        self.mark_me = {}
        for edge in self.eid:
            self.mark_me[edge] = set()
        # check every pair (not very optimized...)
        for e1 in self.eid:
            for e2 in self.eid:
                # to avoid double-checking?
                if self.eid[e1] < self.eid[e2]:
                    continue
                # check if they're mutually exclusive
                if (not self.e_reach_mat[e1][self.eid[e2]]
                        and not self.e_reach_mat[e2][self.eid[e1]]):
                    # mutually exclusive
                    self.mark_me[e1].add(e2)
                    self.mark_me[e2].add(e1)
        # remove edges that have not mutually exclusive edges
        edges_to_rm = set()
        for edge in self.mark_me:
            if len(self.mark_me[edge]) == 0:
                edges_to_rm.add(edge)
        for edge in edges_to_rm:
            del (self.mark_me[edge])

    def getMarkers(self, edge):
        """Find all positions and allele markers corresponding to that edge

        If this edge is a marker, output it as "present". Also outputs
        the "absent" allele for positions of markers that are mutually
        exclusive with that edge.

        Args:
            edge : the edge name ("node1_node2")

        Returns: dict with edge -> allele ("present" or "absent")
        """
        res = {}
        # return nothing if that edge is not a marker
        if edge not in self.mark_me:
            return res
        # this edge's starting node
        res[edge] = 'present'
        for medge in self.mark_me[edge]:
            res[medge] = 'absent'
        return res

    def getMarkersPos(self, edge):
        """Find all positions and allele markers corresponding to that edge

        If this edge is a marker, outputs its position (starting node)
        and allele (edge). Also outputs the "absent" allele for
        positions of markers that are mutually exclusive with that
        edge.

        Args:
            edge : the edge name ("node1_node2")

        Returns: dict with position (node) -> allele (edge name or "absent")
        """
        res = {}
        # return nothing if that edge is not a marker
        if edge not in self.mark_me:
            return res
        # this edge's starting node
        edge_pos = edge.split('_')[0]
        res[edge_pos] = edge
        for medge in self.mark_me[edge]:
            medge_pos = medge.split('_')[0]
            if edge_pos != medge_pos:
                # return absent for other positions
                res[medge_pos] = 'absent'
        return res


class ModuleDiplotyper:
    def __init__(self, include_buffer=True):
        # list of confident module adjacencies from the reads
        # (e.g. [flankl, 2, buffer, 3], [3, flankr])
        self.sread_adj = []
        self.clusters = set()
        # list of possible diplotype
        self.dips = {}
        self.non_cl = ['flankl', 'flankr']
        if include_buffer:
            self.non_cl.append('buffer')

    def print(self):
        print('ModuleDiplotyper with {} clusters, {} informative '
              'read adjacencies, {} possible '
              'diplotypes.'.format(len(self.clusters),
                                   len(self.sread_adj),
                                   len(self.dips)))
        adj_sum = {}
        for adj in self.sread_adj:
            for adji in range(len(adj) - 1):
                adj_e = '{} - {}'.format(adj[adji], adj[adji + 1])
                if adj_e not in adj_sum:
                    adj_sum[adj_e] = 0
                adj_sum[adj_e] += 1
        print(adj_sum)

    def addRead(self, sread_list, sread_info):
        cur_adj = []
        for sreadn in sread_list:
            cl = sread_info[sreadn].cluster
            if sread_info[sreadn].cluster_conf > 0:
                cur_adj.append(cl)
            else:
                if len(cur_adj) > 1:
                    self.sread_adj.append(cur_adj)
                cur_adj = []
            if cl in self.non_cl or cl in self.clusters:
                continue
            self.clusters.add(cl)
        # potentially save the last adjacency
        if len(cur_adj) > 1:
            self.sread_adj.append(cur_adj)

    def prepareDiplotypes(self):
        # first list all possible haplotypes in the form
        # flank-mod-(buffer-mod-)*flank
        haps_prog = [['flankl']]
        haps_done = {}
        while len(haps_prog) > 0:
            hap = haps_prog.pop()
            if hap[-1] in ['flankl', 'buffer']:
                # add a cluster if not already present
                for cl in self.clusters:
                    if cl in hap:
                        continue
                    nhap = list(hap)
                    nhap.append(cl)
                    haps_prog.append(nhap)
            elif hap[-1] in self.clusters:
                # add either a buffer (if there is a buffer) and continue
                if 'buffer' in self.non_cl:
                    nhap = list(hap)
                    nhap.append('buffer')
                    haps_prog.append(nhap)
                else:
                    for cl in self.clusters:
                        if cl in hap:
                            continue
                        nhap = list(hap)
                        nhap.append(cl)
                        haps_prog.append(nhap)
                # or the right flank and stop
                hap.append('flankr')
                # save the haplotype with other using the same clusters.
                hap_name = []
                for cl in self.clusters:
                    if cl in hap:
                        hap_name.append('1')
                    else:
                        hap_name.append('0')
                hap_name = ''.join(hap_name)
                if hap_name not in haps_done:
                    haps_done[hap_name] = []
                haps_done[hap_name].append(hap)
        # list all mutually exclusive pairs
        hap_used = set()
        for hap_name in haps_done:
            # find the name of the haplotype groups using the other
            # clusters, the complement of this one
            hap_comp = []
            for clint in hap_name:
                if clint == '0':
                    hap_comp.append('1')
                else:
                    hap_comp.append('0')
            hap_comp = ''.join(hap_comp)
            # skip if we've already used that haplotype or it's not enumerated
            if hap_comp in hap_used or hap_comp not in haps_done:
                continue
            # otherwise enumerate all the diplotype pairs
            dip_ii = 0
            for h1 in haps_done[hap_name]:
                for h2 in haps_done[hap_comp]:
                    dip_name = '{}_{}_{}'.format(hap_name, hap_comp, dip_ii)
                    dip_ii += 1
                    self.dips[dip_name] = {}
                    self.dips[dip_name]['hap1'] = h1
                    self.dips[dip_name]['hap2'] = h2
            hap_used.add(hap_name)
            hap_used.add(hap_comp)

    def scoreDiplotypes(self, dip_name):
        hap1 = self.dips[dip_name]['hap1']
        hap2 = self.dips[dip_name]['hap2']
        # count number of adjacencies that match?
        score = 0
        for adj in self.sread_adj:
            if anyExactMatch(adj, hap1) or anyExactMatch(adj, hap2):
                score += 1
        return score

    def findBestDiplotype(self):
        # score all the diplotypes
        best_dip = None
        best_dip_score = None
        for dip_name in self.dips:
            score = self.scoreDiplotypes(dip_name)
            if best_dip is None or best_dip_score < score:
                best_dip = dip_name
                best_dip_score = score
        # return the best one (highest score)
        return self.dips[best_dip]


def anyExactMatch(query, target):
    # helper function to check for any exact match
    for t_pos in range(len(target) - len(query)):
        if query[0] == target[t_pos]:
            # check the other positions
            all_match = True
            for offset in range(len(query)-1):
                if query[offset + 1] != target[t_pos + offset + 1]:
                    all_match = False
                    break
            if all_match:
                return True
    return False


def matchBipartite(pair_cost):
    n = len(pair_cost)
    # sort pairs by increasing cost
    cost_d = {}
    for i1 in range(n):
        for i2 in range(n):
            cost_d[(i1, i2)] = pair_cost[i1][i2]
    pairs_sorted = sorted(cost_d, key=lambda k: cost_d[k])
    matching = {}
    matching[0] = []
    matching[1] = []
    matching['cost'] = 0
    for i1, i2 in pairs_sorted:
        if i1 in matching[0] or i2 in matching[1]:
            continue
        matching[0].append(i1)
        matching[1].append(i2)
        matching['cost'] += pair_cost[i1][i2]
        if len(matching[1]) == n:
            # we've assigned all pairs
            break
    # convert to dict with matching
    res = {}
    for ii in range(n):
        res[matching[0][ii]] = matching[1][ii]
    return res


def matchBipartiteOpt(pair_cost):
    n = len(pair_cost)
    cur_cands = []
    matching = {}
    matching[0] = []
    matching[1] = []
    matching['cost'] = 0
    cur_cands.append(matching)
    final_cands = []
    while len(cur_cands) > 0:
        matching = cur_cands.pop()
        for i1 in range(n):
            for i2 in range(n):
                if i1 in matching[0] or i2 in matching[1]:
                    continue
                new_matching = matching.copy()
                new_matching[0].append(i1)
                new_matching[1].append(i2)
                new_matching['cost'] += pair_cost[i1][i2]
                if len(new_matching[1]) == n:
                    # we've assigned all pairs
                    final_cands.append(new_matching)
                else:
                    cur_cands.append(new_matching)
    # find best matching (minimum cost)
    best_matching = None
    best_cost = None
    for matching in final_cands:
        if best_cost is None or matching['cost'] < best_cost:
            best_cost = matching['cost']
            best_matching = matching
    # convert to dict with matching
    res = {}
    for ii in range(n):
        res[best_matching[0][ii]] = best_matching[1][ii]
    return res


def mapReadsToHap(path, reads, max_node_gap=10):
    # record the node position on the path for each alignment candidate
    path_pos = {}
    for pii, nod in enumerate(path):
        if nod in path_pos:
            path_pos[nod].append(pii)
        else:
            path_pos[nod] = [pii]
    read_alns = {}
    for readn in reads:
        read = reads[readn]
        alns = []
        for rii, nod in enumerate(read):
            # node in read but not in path: skip
            if nod not in path_pos:
                continue
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
        best_aln = None
        for aln in alns:
            if best_aln is None or len(aln) > best_score:
                best_score = len(aln)
                best_aln = aln
        if best_score is None:
            continue
        # prepare output dict
        read_alns[readn] = {}
        read_alns[readn]['score'] = best_score
        # convert to alignment track
        nread_edits = 0
        hpos = []
        aln_type = []
        aln_pos = 0
        for nod in read:
            # print('node {}; aln_pos {}'.format(nod, aln_pos))
            if (nod not in path_pos or aln_pos >= len(best_aln)
                    or best_aln[aln_pos] not in path_pos[nod]):
                # either in read and not path or not the next match
                nread_edits += 1
            else:
                # match the path
                cur_hpos = best_aln[aln_pos]
                if nread_edits > 0:
                    # we need to add the read-specific edits before
                    # adding this match
                    if len(hpos) == 0:
                        # inject "mismatch" node at previous position
                        hpos.append(cur_hpos - 1)
                    else:
                        hpos.append((hpos[-1] + cur_hpos) / 2)
                    aln_type.append('mismatch')
                    nread_edits = 0
                # add the match
                hpos.append(cur_hpos)
                aln_type.append('match')
                # move to next aligned position
                aln_pos += 1
        read_alns[readn]['aln_pos'] = hpos
        read_alns[readn]['aln_type'] = aln_type
    return read_alns
