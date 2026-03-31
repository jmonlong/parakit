import networkx as nx


def biClusterSubreads(nodes, reads, min_read_support=3, max_cls=50,
                      max_cycles=4, max_haps=50, verbose=False):
    """Main function to cluster subreads by bi-clustering iteratively

    Here we create candidate sub-paths by by-clustering the reads
    iteratively until no marker/site have diverging read support (with
    minimum support). The function to do this is actually
    biClusterSubreadsIter. Here we call this function with increasing
    stringency (minimum read support) until we get a number of
    candidate sub-paths that is below max_cls.

    Args:
        nodes : dict with node informations
        reads : a Reads object
        min_read_support : minimum reads support to define markers
        max_cls : maximum number of clusters/subpaths to return.
        max_cycles : maximum number of cycles allowed
        max_haps : maximum number of haplotype to enumerate
        verbose : print some informative messages?

    Returns: a pair with: a Subreads object, a list of all subpath candidates
    """
    # start with the provided minimum support. Will be incremented if necessary
    min_rsupp = min_read_support
    # count the attempts to show a warning if too many attempts
    nattempt = 1
    # keep trying with increasing minimum read support
    while True:
        sreads, sread_cls = biClusterSubreadsIter(nodes, reads,
                                                  min_read_support=min_rsupp,
                                                  max_cls=max_cls,
                                                  verbose=verbose)
        if len(sread_cls) <= max_cls:
            # stop if we have a manageable number of subread clusters
            break
            # return (sreads, sread_cls)
        else:
            # warn if we've tried several times already
            if nattempt > 3:
                print("Warning: more than {} module-read clusters found "
                      "with multiple marker support threshold "
                      "(currently {}). The reads might be too "
                      "short/noisy for accurate "
                      "diplotyping.".format(max_cls, min_rsupp))
            # if too many subreads clusters, rerun with more
            # stringent read support
            min_rsupp += 1
            nattempt += 1
            if verbose:
                print("\tToo many subread clusters ({}). "
                      "Rerunning with min "
                      "read support: {}".format(len(sread_cls),
                                                min_rsupp))
    # then enumerate path candidates from subread clusters
    paths = sreads.enumAlleles(sread_cls, max_cycles=4,
                               max_haps=max_haps,
                               min_read_support=min_read_support,
                               verbose=verbose)
    # enumerate all pairs
    path_names = list(paths.keys())
    hap_pairs_names = []
    for mii in range(len(path_names)):
        for mjj in range(mii, len(path_names)):
            hap_pairs_names.append([path_names[mii],
                                    path_names[mjj]])
    return {'paths': paths, 'hap_pairs': hap_pairs_names}


def biClusterSubreadsIter(nodes, reads, min_read_support=3,
                          max_cls=50, verbose=False):
    """Run one iteration of the subreads bi-clustering method

    Bi-cluster reads until there are no divergent markers (with at
    least min_read_support supporting reads).

    Args:
        nodes : dict with node information
        reads : a Reads object
        min_read_support : minimum reads support to define markers
        max_cls : maximum number of clusters/subpaths to return.
        verbose : print some informative messages?

    Returns: a pair with: a Subreads object, a list of all subpath candidates

    """
    # init subreads
    sreads = Subreads()
    sreads.splitReads(reads, nodes)
    sreads.computeCoverage()
    if verbose:
        print('\t\tClustering subreads...')
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
        csreads.findMarkers(min_read_support=min_read_support)
        # if some supported markers
        if csreads.nbMarkers() > 0:
            if verbose:
                print('{} marker(s) found, '
                      'splitting'.format(csreads.nbMarkers()))
            # split the reads in two
            csreads.biClusterReads()
            sreads_list.append(csreads.subsetByCluster(0))
            sreads_list.append(csreads.subsetByCluster(1))
        # iterate until no clear markers suggesting two alleles
        if verbose and len(sreads_list) % 100 == 0 and len(sreads_list) > 0:
            print('\t\tWatchdog: {} clusters in '
                  'progress ({} finished).'.format(len(sreads_list),
                                                   len(sreads_list_final)))
        # abort if we already have too many clusters
        if len(sreads_list_final) > max_cls:
            return ((sreads, sreads_list_final))
    return (sreads, sreads_list_final)


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


class Subreads:
    """Reads split at the cycle boundaries

    For the haplotype reconstruction, we first try to reconstruct the
    different alleles of the module (the cycling part of the
    pangenome). To do that we work with subreads, i.e. the full reads
    cut in pieces at the cycle's boundaries. This class helps
    manipulating them.

    Attributes:
        cyc_edge : list of the two cycling edges
        nsucs : dict node -> default successor node (for consensus generation)
        flankn : lists the two paths corresponding to the flanks
        path : dict subread -> path through the graph
        cyc : dict with subreads that cycle at least once
        flankl : dict with subreads that overlap the left flank
        flankr : dict with subreads that overlap the right flank
        ecov : dict edge [n][m] -> subread coverage (absent means not covered)
        ecov_prev : same as ecov but from the previous iteration
        markers : dict with nodes that are differently covered by subreads
        site_markers : lists sites/edges with different subread cover (for EM)
        sparts : dict subread -> cluster (0 or 1)

    Methods:
        splitReads : split Reads and fill the Subreads object
        computeCoverage : update the edge coverage ecov
        findMarkers :
        nbMarkers :
        biClusterReads :
        subsetByCluster :
        getConsensus :
        assignReads :
        enumAlleles :
        clusterEM :
        enumAllelePair :
        buildClusterEdges :
        clustersToPaths :
        print :
    """

    def __init__(self):
        # saves the cycling "edge"
        self.cyc_edge = ['', '']
        # saves default successor nodes
        self.nsucs = {}
        # saves the reference flanking nodes
        self.flankn = [[], []]
        # saves the path of each subread through the graph (node sequence)
        self.path = {}
        # saves subreads that cycle back
        self.cyc = {}
        # saves subreads to start on the left flank
        self.flankl = {}
        # saves subreads to continue to the right flank
        self.flankr = {}
        # saves the coverage (nb subreads) of each edge in the graph
        # (absent means not covered)
        self.ecov = {}
        # same but for the previous iteration
        # (useful if not covered with current subread subset)
        self.ecov_prev = {}
        # saves the node markers where subreads differ to use for clustering
        self.markers = {}
        # saves the bi-partition for each subread
        self.sparts = {}

    def splitReads(self, reads, nodes):
        # before looking at reads, look at the nodes and save some information
        # for example, define "default" edges from each node in nsucs. If we
        # have no/limited support during the inference but want to continue
        # traversing a path, we'll use that edge
        # here, also save the edges to the next reference (used for fast flank
        # reconstruction)
        ed = {}
        ed_rev = {}
        for nod in nodes:
            # also remember the cycle's boundaries
            if nodes[nod]['class'] == 'cyc_r':
                self.cyc_edge[0] = nod
            elif nodes[nod]['class'] == 'cyc_l':
                self.cyc_edge[1] = nod
            # save default successor for each node
            if 'sucs' in nodes[nod]:
                # to save the "best" node
                bnod = ''
                for nnod in nodes[nod]['sucs']:
                    nnod_c12 = bnod_c12 = 0
                    if bnod != '':
                        nnod_c12 = nodes[nnod]['c1'] + nodes[nnod]['c2']
                        bnod_c12 = nodes[bnod]['c1'] + nodes[bnod]['c2']
                    # keep the nodes with most haplotype/ref support
                    cond = bnod == ''
                    cond = cond or (nodes[nnod]['ref'] > 0 and
                                    nnod_c12 > bnod_c12)
                    cond = cond or (nodes[nnod]['ref'] > nodes[bnod]['ref'])
                    if cond:
                        bnod = nnod
                    # save edges to reference path for flank approximation
                    if nodes[nnod]['class'] == 'ref':
                        ed[nod] = nnod
                    if nodes[nod]['class'] == 'ref':
                        ed_rev[nnod] = nod
                self.nsucs[nod] = bnod
        # find reference region for left flank by starting at cycling edge
        flank_l_path = [self.cyc_edge[1]]
        while flank_l_path[0] in ed_rev:
            onod = ed_rev[flank_l_path[0]]
            if onod in flank_l_path:
                break
            flank_l_path = [onod] + flank_l_path
        flank_r_path = [self.cyc_edge[0]]
        while flank_r_path[-1] in ed:
            onod = ed[flank_r_path[-1]]
            if onod in flank_r_path:
                break
            flank_r_path.append(ed[flank_r_path[-1]])
        self.flankn = [flank_l_path, flank_r_path]
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
            for nod in reads.path[readn]:
                if nod in self.cyc_edge:
                    subreads.append([])
                subreads[-1].append(nod)
            # guess what type of subreads we have
            rpos_flankl_end = nodes[self.cyc_edge[1]]['rpos_min']
            rpos_flankr_start = nodes[self.cyc_edge[0]]['rpos_max']
            subreads_t = []
            for sri, subread in enumerate(subreads):
                subread_t = ''
                for nod in subread:
                    if nodes[nod]['rpos_min'] < rpos_flankl_end:
                        # definitely the left flank
                        subread_t = 'flankl'
                        break
                    elif nodes[nod]['rpos_max'] > rpos_flankr_start:
                        # definitely the right flank
                        subread_t = 'flankr'
                        break
                    elif nod == self.cyc_edge[1]:
                        # definitely a subread to analyze
                        subread_t = 'subread'
                        break
                    elif (nodes[nod]['class'] == 'ref' and
                          nodes[nod]['rpos_min'] > rpos_flankl_end and
                          nodes[nod]['rpos_max'] < rpos_flankr_start):
                        # definitely within a buffer region
                        subread_t = 'buffer'
                        continue
                # if we haven't guessed yet, try to use the next subread
                # to figure out
                if subread_t == '' and len(subreads) > sri + 1:
                    # there is a next subread
                    if subreads[sri+1][0] == self.cyc_edge[0]:
                        # next subreads exit the region of interest
                        # hence this is a subread to analyze
                        subread_t = 'subread'
                    elif subreads[sri+1][0] == self.cyc_edge[1]:
                        # next subreads enters the region of interest
                        # hence this is a buffer region
                        subread_t = 'buffer'
                    else:
                        subread_t = 'unknown'
                if subread_t == '' and len(subreads) == sri + 1:
                    # there is no next subread, might be a read completely
                    # within the region of interest
                    subread_t = 'subread'
                subreads_t.append(subread_t)
            # save the subreads to analyze,
            # flag them if touching the flanks or cycling
            sreadc = 0
            sreadn = '{}_{}'.format(readn, sreadc)
            for sbi, subread in enumerate(subreads):
                if subreads_t[sbi] == 'subread':
                    self.path[sreadn] = subread
                    # check previous subread
                    if sbi > 0 and subreads_t[sbi-1] == 'flankl':
                        self.flankl[sreadn] = True
                    # check next subreads
                    nsrs = subreads_t[(sbi+1):]
                    if len(nsrs) > 0:
                        if 'subread' in nsrs:
                            self.cyc[sreadn] = True
                        elif nsrs[0] == 'flankr':
                            self.flankr[sreadn] = True
                    # prepare next subread
                    sreadc += 1
                    sreadn = '{}_{}'.format(readn, sreadc)

    def computeCoverage(self):
        for sreadn in self.path:
            for pos, nod in enumerate(self.path[sreadn][:-1]):
                nnod = self.path[sreadn][pos+1]
                if nod not in self.ecov:
                    self.ecov[nod] = {}
                if nnod not in self.ecov[nod]:
                    self.ecov[nod][nnod] = 0
                self.ecov[nod][nnod] += 1

    def findMarkers(self, min_read_support=3):
        self.markers = {}
        self.site_markers = []
        # start from within the collapsed part of the pangenome
        cnod = self.cyc_edge[1]
        while cnod != self.cyc_edge[0]:
            if cnod not in self.ecov or len(self.ecov[cnod]) == 0:
                # no edge coverage, look at the coverage from previous rounds
                # first check if we know that "previous" coverage
                if cnod not in self.ecov_prev:
                    self.ecov_prev[cnod] = {}
                if len(self.ecov_prev[cnod]) == 0:
                    # if missing, set coverage to 1
                    self.ecov_prev[cnod][self.nsucs[cnod]] = 1
                # set edge coverage to most supported edge in previous round
                best_nnod = ''
                best_supp = 0
                for nnod in self.ecov_prev[cnod]:
                    if self.ecov_prev[cnod][nnod] > best_supp:
                        best_nnod = nnod
                        best_supp = self.ecov_prev[cnod][nnod]
                # save most supported edge in current coverage
                self.ecov[cnod] = {}
                self.ecov[cnod][best_nnod] = best_supp
            # check support for each outgoing edge
            nsupp = []
            rsupp = []
            best_nnod = ''
            best_supp = 0
            for nnod in self.ecov[cnod]:
                if self.ecov[cnod][nnod] >= best_supp:
                    best_nnod = nnod
                    best_supp = self.ecov[cnod][nnod]
                if self.ecov[cnod][nnod] >= min_read_support:
                    nsupp.append(nnod)
                    rsupp.append(self.ecov[cnod][nnod])
            # save those nodes as markers if more than one supported edge
            if len(nsupp) > 1:
                rsupp.sort(reverse=True)
                for nnod in nsupp:
                    if nnod in self.markers:
                        self.markers[nnod] = max(self.markers[nnod],
                                                 rsupp[1])
                    else:
                        self.markers[nnod] = rsupp[1]
                self.site_markers.append({'prev': cnod, "next": nsupp})
            # move to next node
            cnod = best_nnod

    def nbMarkers(self):
        return (len(self.markers))

    def biClusterReads(self):
        # prepare marker signature for each read
        msig = {}
        for sreadn in self.path:
            msig[sreadn] = {}
            for nod in self.path[sreadn]:
                if nod in self.markers:
                    msig[sreadn][nod] = True
        # compare reads pairwise and build network
        G = nx.Graph()
        sread_names = list(self.path.keys())
        for sr1 in range(len(sread_names)-1):
            sreadn1 = sread_names[sr1]
            mark1 = list(msig[sreadn1].keys())
            for sr2 in range(sr1 + 1, len(sread_names)):
                sreadn2 = sread_names[sr2]
                # compute nb of markers in common
                n_com_mark = 0
                for mark in mark1:
                    if mark in msig[sreadn2]:
                        n_com_mark += 1
                if n_com_mark > 0:
                    G.add_edge(sreadn1, sreadn2, weight=n_com_mark)
        # find best cut
        sparts = nx.community.kernighan_lin_bisection(G, seed=123)
        # parse groups
        for cl in [0, 1]:
            for sreadn in sparts[cl]:
                self.sparts[sreadn] = cl

    def subsetByCluster(self, cl):
        sreads_cl = Subreads()
        sreads_cl.cyc_edge = self.cyc_edge
        for sreadn in self.path:
            if sreadn in self.sparts and self.sparts[sreadn] == cl:
                sreads_cl.path[sreadn] = self.path[sreadn]
        sreads_cl.ecov_prev = self.ecov
        sreads_cl.nsucs = self.nsucs
        sreads_cl.computeCoverage()
        return (sreads_cl)

    def getConsensus(self):
        # start from within the collapsed part of the pangenome
        cnod = self.cyc_edge[1]
        cons_path = [cnod]
        while cnod != self.cyc_edge[0]:
            if cnod not in self.ecov or len(self.ecov[cnod]) == 0:
                if cnod not in self.ecov_prev:
                    self.ecov_prev[cnod] = {}
                if len(self.ecov_prev[cnod]) == 0:
                    self.ecov_prev[cnod][self.nsucs[cnod]] = 1
                # no coverage, look at the coverage from previous round
                best_nnod = ''
                best_supp = 0
                for nnod in self.ecov_prev[cnod]:
                    if self.ecov_prev[cnod][nnod] > best_supp:
                        best_nnod = nnod
                        best_supp = self.ecov_prev[cnod][nnod]
                # save most supported edge in current coverage
                self.ecov[cnod] = {}
                self.ecov[cnod][best_nnod] = best_supp
            # check support for each outgoing edge
            best_nnod = ''
            best_supp = 0
            for nnod in self.ecov[cnod]:
                if self.ecov[cnod][nnod] >= best_supp:
                    best_nnod = nnod
                    best_supp = self.ecov[cnod][nnod]
            # move to next node
            cnod = best_nnod
            cons_path.append(cnod)
        return (cons_path)

    def assignReads(self, cl_paths):
        cl_assign = {}
        rscores = {}
        for cl_path in cl_paths:
            rscores_cl = pathReadGraphAlign(cl_path, self.path)
            for sreadn in rscores_cl:
                if sreadn not in rscores:
                    rscores[sreadn] = []
                rscores[sreadn].append(rscores_cl[sreadn])
        for sreadn in rscores:
            cl_assign[sreadn] = []
            max_rscore = 0
            for rscore in rscores[sreadn]:
                max_rscore = max(rscore, max_rscore)
            for cl_ii, rscore in enumerate(rscores[sreadn]):
                if rscore > max_rscore * .98:
                    cl_assign[sreadn].append(cl_ii)
        return (cl_assign)

    def enumAlleles(self, cluster_list, max_cycles=3, max_haps=500,
                    min_read_support=3, verbose=False):
        # print the cluster sizes
        if verbose:
            print('\t\tEnumerating haplotypes from {} '
                  'clusters...'.format(len(cluster_list)))
        # make a consensus path for each cluster
        cl_paths = []
        for cl in cluster_list:
            cl_paths.append(cl.getConsensus())
        # assign each subread to best cluster (duplicate ties)
        cl_assign = self.assignReads(cl_paths)
        # link subread clusters using read splitting information
        block_e = self.buildClusterEdges(cl_assign)
        # enumerate cluster sequence
        cur_paths = [{'path': ['flankl'], 'support': []}]
        final_paths = []
        while len(cur_paths) > 0 and len(final_paths) < 20000:
            cpath = cur_paths.pop(0)
            path = cpath['path']
            if path[-1] == 'flankr':
                final_paths.append(cpath)
                continue
            if path[-1] in block_e:
                for cl in block_e[path[-1]]:
                    # TODO improve threshold for block edges?
                    # if block_e[path[-1]][cl] >= min_read_support and \
                    #    path.count(cl) <= max_cycles:
                    if block_e[path[-1]][cl] >= min_read_support and \
                       len(path) < max_cycles + 1:
                        npath = {}
                        npath['path'] = path + [cl]
                        npath['support'] = cpath['support']
                        npath['support'] += [block_e[path[-1]][cl]]
                        cur_paths.insert(0, npath)
            if verbose and len(cur_paths) % 5000 == 0 and len(cur_paths) > 0:
                print('\t\tWatchdog: {} haplotype candidate in '
                      'progress ({} finished).'.format(len(cur_paths),
                                                       len(final_paths)))
        # enumerate (node) paths
        return (self.clustersToPaths(final_paths, cl_paths, max_haps=max_haps,
                                     verbose=verbose))

    def buildClusterEdges(self, cl_assign):
        block_e = {'flankl': {}}
        for sreadn in cl_assign:
            # subreads starting in the "left" flank
            if sreadn in self.flankl:
                for clii in cl_assign[sreadn]:
                    if clii not in block_e['flankl']:
                        block_e['flankl'][clii] = 0
                    block_e['flankl'][clii] += 1
            # subreads ending to the "right" flank
            if sreadn in self.flankr:
                for clii in cl_assign[sreadn]:
                    if clii not in block_e:
                        block_e[clii] = {'flankr': 0}
                    block_e[clii]['flankr'] += 1
            # subreads that cycle back somewhere
            if sreadn in self.cyc:
                # find where is the next subread (if anywhere)
                sr_id = int(sreadn.split('_')[-1])
                readn = '_'.join(sreadn.split('_')[:-1])
                next_sr = readn + '_' + str(sr_id+1)
                if next_sr in cl_assign:
                    for cl_s in cl_assign[sreadn]:
                        if cl_s not in block_e:
                            block_e[cl_s] = {'flankr': 0}
                        for cl_e in cl_assign[next_sr]:
                            if cl_e not in block_e[cl_s]:
                                block_e[cl_s][cl_e] = 0
                            block_e[cl_s][cl_e] += 1
        return (block_e)

    def clustersToPaths(self, paths, cl_paths,
                        max_haps=None, verbose=False):
        paths_n = {}

        def getSupport(ii):
            return (-float(sum(paths[ii]['support'])) / len(paths[ii]))
        # sort paths by support
        sorted_idx = sorted(list(range(len(paths))), key=getSupport)
        mod_n = {}
        warn_printed = False
        for ii in sorted_idx:
            path = paths[ii]['path']
            if len(path) not in mod_n:
                mod_n[len(path)] = 0
            mod_n[len(path)] += 1
            # if too many paths, keep the top ones
            # (higest mean block-junction support)
            if max_haps is not None and mod_n[len(path)] >= max_haps:
                if verbose and not warn_printed:
                    print("\t\t\tToo many haplotypes ({}). Attempting to keep "
                          "the top {} for each module "
                          "number.".format(len(paths), max_haps))
                    warn_printed = True
                continue
            path_exp = []
            for cl in path:
                if cl == 'flankl':
                    path_exp += self.flankn[0]
                elif cl == 'flankr':
                    path_exp += self.flankn[1]
                else:
                    path_exp += cl_paths[cl]
            paths_n['_'.join([str(p) for p in path])] = path_exp
        return (paths_n)

    def print(self):
        print('{} subreads, {} markers, '
              '{} partitioned subreads'.format(len(self.path),
                                               len(self.markers),
                                               len(self.sparts)))
