import networkx as nx
import parakit.parakit_path as pkpath
import random


class ReadPosList:
    def __init__(self):
        # record reads and position(s)
        self.read_to_pos = {}
        # record if the read-position is after the read having cycled once
        self.readpos_to_cyc = {}
        # list of reads to exclude from future updates
        # (because they were removed earlier)
        self.exc_list = {}

    def addReadPos(self, read_name, read_pos, has_cycled=False):
        if read_name not in self.read_to_pos:
            self.read_to_pos[read_name] = []
        self.read_to_pos[read_name].append(read_pos)
        self.readpos_to_cyc[read_name + '_' + str(read_pos)] = has_cycled

    def addReadPosList(self, rpl):
        # add read positions
        for readn in rpl.read_to_pos:
            if readn not in self.read_to_pos:
                self.read_to_pos[readn] = []
            for read_pos in rpl.read_to_pos[readn]:
                self.read_to_pos[readn].append(read_pos)
        # add cycles
        for rp in rpl.readpos_to_cyc:
            self.readpos_to_cyc[rp] = rpl.readpos_to_cyc[rp]

    def hasCycled(self, read_name, read_pos):
        return (self.readpos_to_cyc[read_name + '_' + str(read_pos)])

    def updateReads(self, other_rp, add_reads=False, max_pos_diff=-1,
                    max_rm_reads=5, verbose=False):
        radd = []
        for readn in other_rp.read_to_pos:
            # if absent, skip except is we want to add new reads
            if readn not in self.read_to_pos or \
               len(self.read_to_pos[readn]) == 0:
                if add_reads:
                    # skip if the reads was removed multiple times
                    if readn in self.exc_list and \
                       self.exc_list[readn] > max_rm_reads:
                        if verbose:
                            print(readn + ' excluded')
                        continue
                    # find smallest position
                    min_rpos = other_rp.read_to_pos[readn][0]
                    for rpos in other_rp.read_to_pos[readn]:
                        min_rpos = min(min_rpos, rpos)
                    # add new read if it didn't loop before
                    if not other_rp.hasCycled(readn, min_rpos):
                        self.addReadPos(readn, min_rpos)
                        radd.append(readn)
                continue
            # if we don't want to check positions, increment all positions
            if max_pos_diff < 0:
                for ii in range(len(self.read_to_pos[readn])):
                    self.read_to_pos[readn][ii] += 1
                continue
            # otherwise, make sure the pos are not too far from each other
            for ii, rpos in enumerate(self.read_to_pos[readn]):
                for rpos2 in other_rp.read_to_pos[readn]:
                    if abs(rpos - rpos2) <= max_pos_diff:
                        if verbose and False:
                            vtp = '{} pos updated: {} -> {}'
                            print(vtp.format(readn,
                                             self.read_to_pos[readn][ii],
                                             rpos2))
                        self.read_to_pos[readn][ii] = rpos2
        if verbose:
            for readn in radd:
                print(readn + ' added')

    def countIntersect(self, other_rp, max_pos_diff=-1):
        ninter = 0
        for readn in self.read_to_pos:
            # sanity check that the array of positions is not empty
            if len(self.read_to_pos[readn]) == 0:
                continue
            # if absent from other list, skip
            if readn not in other_rp.read_to_pos or \
               len(other_rp.read_to_pos[readn]) == 0:
                continue
            # if we don't want to check positions, count as matched
            if max_pos_diff < 0:
                ninter += 1
                continue
            # otherwise, make sure the pos are not too far from each other
            for rpos in self.read_to_pos[readn]:
                for rpos2 in other_rp.read_to_pos[readn]:
                    if abs(rpos - rpos2) <= max_pos_diff:
                        ninter += 1
        return (ninter)

    def removeReads(self, other_rp, max_pos_diff=-1, verbose=False):
        nrem = []
        for readn in other_rp.read_to_pos:
            # if absent, skip
            if readn not in self.read_to_pos or \
               len(self.read_to_pos[readn]) == 0:
                continue
            # if we don't want to check positions, remove no matter what
            if max_pos_diff < 0:
                self.read_to_pos[readn] = []
                nrem.append(readn)
                continue
            # otherwise, make sure the positions are not to far from each other
            for rpos in self.read_to_pos[readn]:
                for rpos2 in other_rp.read_to_pos[readn]:
                    if abs(rpos - rpos2) <= max_pos_diff:
                        self.read_to_pos[readn] = []
                        nrem.append(readn)
        # keep track of which reads were removed
        # we've decided they don't support this haplotype, we don't want to
        # add them again later
        for readn in nrem:
            if readn not in self.exc_list:
                self.exc_list[readn] = 0
            self.exc_list[readn] += 1
        # print the names of removed reads?
        if verbose:
            for readn in nrem:
                print(readn + ' removed')

    def getReads(self):
        all_reads = {}
        for readn in self.read_to_pos:
            if len(self.read_to_pos[readn]) > 0:
                all_reads[readn] = True
        return (list(all_reads.keys()))

    def hasRead(self, read_name):
        return (read_name in self.read_to_pos and
                len(self.read_to_pos[read_name]) > 0)

    def getReadNumber(self):
        nreads = 0
        for readn in self.read_to_pos:
            if len(self.read_to_pos[readn]) > 0:
                nreads += 1
        return (nreads)

    def print(self):
        print('{} reads.'.format(self.getReadNumber()))


class Reads:
    def __init__(self):
        self.edge_to_readpos = {}
        self.nsuc = {}
        self.path = {}
        self.startpos = {}
        self.endpos = {}
        self.readpos = {}

    def addRead(self, read_name, path, cyc_nodes=[],
                startpos=[], endpos=[], readpos=[]):
        # save path and position in the sequenced read
        self.path[read_name] = path
        if len(startpos) > 0:
            self.startpos[read_name] = startpos
        if len(endpos) > 0:
            self.endpos[read_name] = endpos
        if len(readpos) > 0:
            self.readpos[read_name] = readpos
        # save edges support
        has_cycled = False
        for pos in range(len(path)-1):
            if path[pos] in cyc_nodes:
                has_cycled = True
            ename = path[pos] + '_' + path[pos+1]
            # add read-positions for this edge
            if ename not in self.edge_to_readpos:
                self.edge_to_readpos[ename] = ReadPosList()
            self.edge_to_readpos[ename].addReadPos(read_name, pos, has_cycled)
            # add node successor information
            if path[pos] not in self.nsuc:
                self.nsuc[path[pos]] = {}
            if path[pos+1] not in self.nsuc[path[pos]]:
                self.nsuc[path[pos]][path[pos+1]] = 0
            self.nsuc[path[pos]][path[pos+1]] += 1

    def getPath(self, read_name):
        if read_name not in self.path:
            return ([])
        return (self.path[read_name])

    def nReads(self):
        return (len(self.path))
    
    def hasRead(self, read_name):
        return (read_name in self.path)

    def getStartPos(self, readn, path_pos):
        if readn in self.startpos and len(self.startpos[readn]) > 0:
            return (self.startpos[readn][path_pos])
        else:
            return ("NA")

    def getEndPos(self, readn, path_pos):
        if readn in self.endpos and len(self.endpos[readn]) > 0:
            return (self.endpos[readn][path_pos])
        else:
            return ("NA")

    def getReadPos(self, readn, path_pos):
        if readn in self.readpos and len(self.readpos[readn]) > 0:
            return (self.readpos[readn][path_pos])
        else:
            return ("NA")

    def getReadPosList(self, snode, enode):
        ename = '{}_{}'.format(snode, enode)
        return (self.edge_to_readpos[ename])

    def getAllReadPosList(self, snode):
        rp = ReadPosList()
        if len(self.nsuc[snode]) == 0:
            return (rp)
        elif len(self.nsuc[snode]) == 1:
            enode = list(self.nsuc[snode])[0]
            ename = '{}_{}'.format(snode, enode)
            return (self.edge_to_readpos[ename])
        else:
            for enode in self.nsuc[snode]:
                ename = '{}_{}'.format(snode, enode)
                rp.addReadPosList(self.edge_to_readpos[ename])
            return (rp)

    def print(self):
        nedge = 0
        nreads = {}
        nnodes = {}
        for ename in self.edge_to_readpos:
            if len(self.edge_to_readpos) > 0:
                nedge += 1
                ereads = self.edge_to_readpos[ename].getReads()
                for readn in ereads:
                    nreads[readn] = True
        for nod in self.nsuc:
            nnodes[nod] = True
            for snod in self.nsuc[nod]:
                nnodes[snod] = True
        tp = "{} reads, {} unique nodes, {} unique edges"
        print(tp.format(len(nreads), len(nnodes), nedge))

    def listSuccessors(self, node_name, min_read_support=3):
        sucs_l = []
        if node_name in self.nsuc:
            for snod in self.nsuc[node_name]:
                if self.nsuc[node_name][snod] >= min_read_support:
                    sucs_l.append(snod)
        return (sucs_l)

    def overlapWithPath(self, path, max_node_gap=10, all_reads=False,
                        in_alns={}, path_pos_s=0):
        # only consider reads overlapping the last node of the path
        last_edge_rp = self.getReadPosList(path[-2], path[-1])
        alns = {}
        if len(in_alns) > 0 and path_pos_s != 0:
            alns = in_alns
        # find reads starting at each possible edge spelled by the current path
        for nii in range(path_pos_s, len(path)-1):
            rpii = self.getReadPosList(path[nii], path[nii+1])
            # check overlap for each of those reads
            for readn in rpii.read_to_pos:
                # skip if not read of interest
                if not all_reads and not last_edge_rp.hasRead(readn):
                    continue
                # rpath = self.getPath(readn)
                if readn not in alns:
                    alns[readn] = []
                # and for each "starting" position
                for read_pos in rpii.read_to_pos[readn]:
                    if len(alns[readn]) == 0:
                        # if we haven't started alignments, start one
                        aln = Alignment(read_pos, nii)
                        alns[readn].append(aln)
                    else:
                        # otherwise try to extend current alignments
                        aln_extended = False
                        for aln in alns[readn]:
                            if aln.extendAl(read_pos + 1,
                                            nii + 1,
                                            max_node_gap=max_node_gap):
                                aln_extended = True
                        # if none were extended, start a new one
                        if not aln_extended:
                            aln = Alignment(read_pos, nii)
                            alns[readn].append(aln)
        return (alns)

    def overlapWithWalks(self, walks, min_prop=.9,
                         min_len=10, all_reads=False):
        # save alignments in a dict: read_name -> array of Alignment objects
        alns = {}
        # loop through the array of walks
        for wi, walk in enumerate(walks):
            # skip that walk is too short/empty
            if len(walk.path) < 2:
                continue
            # overlap reads with that walk
            walns = self.overlapWithPath(walk.path, all_reads=all_reads)
            min_len_w = min(len(walk.path), min_len)
            for readn in walns:
                # make sure this read is in the dict to return
                if readn not in alns:
                    alns[readn] = []
                for aln in walns[readn]:
                    # only save good alignements
                    if aln.propMatch1() >= min_prop and \
                       aln.propMatch2() >= min_prop and\
                       aln.len1() >= min_len_w and \
                       aln.len2() >= min_len_w and \
                       aln.lenLeftSoftclip1() < min_len:
                        # remember which walk this aln is on
                        aln.setWalk(wi)
                        # save to dict to return
                        alns[readn].append(aln)
        return (alns)


class Subreads:
    def __init__(self):
        # saves the cycling "edge"
        self.cyc_edge = ['', '']
        # saves default successor nodes
        self.nsucs = {}
        # saves the reference flanking nodes
        self.flankn = ['', '']
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
        # saves the site markers where subreads differ to use for EM clustering
        self.site_markers = {}
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
                    elif nodes[nod]['class'] == 'ref' and \
                         nodes[nod]['rpos_min'] > rpos_flankl_end and \
                         nodes[nod]['rpos_max'] < rpos_flankr_start:
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

    def keepTopMarkers(self, n_top=10):
        if self.nbMarkers() > n_top:
            mark_ord = sorted(list(self.markers.keys()),
                              key=lambda nn: -self.markers[nn])
            mark_ord = mark_ord[:n_top]
            mark_torm = []
            for nod in self.markers:
                if nod not in mark_ord:
                    mark_torm.append(nod)
            for nod in mark_torm:
                del self.markers[nod]

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
            rscores_cl = pkpath.pathReadGraphAlign(cl_path, self.path)
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

    def clusterEM(self, nprofiles, nalleles=2, verbose=False):
        # prepare markers to use
        nmarkers = 0
        site_markers = {}
        cur_pos = 0
        for smark in self.site_markers:
            mnodes = smark['next']
            if len(mnodes) == nalleles:
                nmarkers += 1
                for val in range(nalleles):
                    en = '{}_{}'.format(smark['prev'], mnodes[val])
                    site_markers[en] = [cur_pos, val]
                cur_pos += 1
        # prepare a map (sread name -> Observation) marker status
        read_prof = {}
        for sreadn in self.path:
            pos_v = []
            val_v = []
            for rpos in range(len(self.path[sreadn])-1):
                en = '{}_{}'.format(self.path[sreadn][rpos],
                                    self.path[sreadn][rpos+1])
                if en in site_markers:
                    pos_v.append(site_markers[en][0])
                    val_v.append(site_markers[en][1])
            obs = Observation(pos_v, val_v)
            # obs.imputeGaps()
            read_prof[sreadn] = obs
        # if verbose:
        #     print('\t\t\t\t{} markers, {} subreads'.format(nmarkers,
        #                                                    len(read_prof)))
        # to record the best set of N profiles
        best_em = None
        # init to a high score, like all reads having all pos wrong
        best_em_ndiff = len(read_prof) * nmarkers
        # try a couple of attempts and keep the best
        for attempt in range(50):
            # init the EM profiles
            em = EM(nmarkers, nprofiles, values=list(range(nalleles)))
            # iterate
            for iter in range(100):
                em.updateWithObs(read_prof)
            em.consensusWithObs(read_prof)
            ndiff = em.consensusWithObs(read_prof)
            if ndiff < best_em_ndiff:
                best_em_ndiff = ndiff
                best_em = em
        # assign each subread to the most similar profile
        # best_em.print()
        for sreadn in read_prof:
            cl = read_prof[sreadn].findBestEMProfile(best_em)
            # print('{} {}'.format(sreadn, cl))
            self.sparts[sreadn] = cl

    def enumAllelePair(self, cluster_list, min_read_support=3,
                       verbose=False):
        # print the cluster sizes
        if verbose:
            cl_size = [len(scl.path) for scl in cluster_list]
            cl_size.sort()
            cl_size = [str(cls) for cls in cl_size]
            print('\t\tEnumerating haplotypes from {} '
                  'clusters (size {})...'.format(len(cluster_list),
                                                 ', '.join(cl_size)))
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
        all_paths = []
        while len(cur_paths) > 0:
            cpath = cur_paths.pop(0)
            path = cpath['path']
            if path[-1] == 'flankr':
                all_paths.append(cpath)
                continue
            if path[-1] in block_e:
                for cl in block_e[path[-1]]:
                    if block_e[path[-1]][cl] >= min_read_support and \
                       path.count(cl) == 0:
                        npath = {}
                        npath['path'] = path + [cl]
                        npath['support'] = cpath['support'] + [block_e[path[-1]][cl]]
                        cur_paths.insert(0, npath)
            if verbose and len(cur_paths) % 5000 == 0 and len(cur_paths) > 0:
                print('\t\tWatchdog: {} haplotype candidate in '
                      'progress ({} finished).'.format(len(cur_paths),
                                                       len(all_paths)))
        # enumerate all compatible pairs
        all_pairs = []
        for p1 in all_paths:
            for p2 in all_paths:
                if p1 == p2:
                    continue
                if len(p1['path']) + len(p2['path']) != 4 + len(cl_paths):
                    continue
                compatible = True
                for cl in p2['path']:
                    if cl == 'flankr' or cl == 'flankl':
                        continue
                    if cl in p1['path']:
                        compatible = False
                        break
                if compatible:
                    all_pairs.append([p1, p2])
        if len(all_pairs) == 0:
            print('Warning: could not find a supported diplotype from {} '
                  'subread clusters.'.format(len(cluster_list)))
            return ({})
        # pick best pair
        best_supp = -1
        best_pair = 0
        for pp in all_pairs:
            # compute the normalized support of the pair
            # supp = sum(pp[0]['support']) / len(pp[0]['path'])
            # supp += sum(pp[1]['support']) / len(pp[1]['path'])
            supp = sum(pp[0]['support'])
            supp += sum(pp[1]['support'])
            # save best one
            if supp > best_supp:
                best_supp = supp
                best_pair = pp
        # enumerate (node) paths
        return (self.clustersToPaths(best_pair, cl_paths, verbose=verbose))

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
                    path_exp += self.flankn[0]
                else:
                    path_exp += cl_paths[cl]
            paths_n['_'.join([str(p) for p in path])] = path_exp
        return (paths_n)

    def print(self):
        print('{} subreads, {} markers, '
              '{} partitioned subreads'.format(len(self.path),
                                               len(self.markers),
                                               len(self.sparts)))


class Observation:
    def __init__(self, positions, values):
        self.pos = positions
        self.val = values

    def print(self):
        out_h = ''
        out_v = ''
        prev_size = 1
        for ii, pos in enumerate(self.pos):
            out_h += ' ' * (4 - prev_size) + str(pos)
            if self.val[ii]:
                out_v += ' ' * (4 - 1) + '+'
            else:
                out_v += ' ' * (4 - 1) + '-'
            prev_size = len(str(pos))
        print(out_h)
        print(out_v)

    def findBestEMProfile(self, em):
        best_ss = 0
        best_count = 0
        for ss in range(em.nstates):
            match_count = 0
            for ii, pp in enumerate(self.pos):
                if em.prof[pp][ss] == self.val[ii]:
                    match_count += 1
            if match_count > best_count:
                best_count = match_count
                best_ss = ss
        return (best_ss)

    def imputeGaps(self):
        new_pos = []
        new_val = []
        # check if gap in positions
        cpos = None
        for ii in range(len(self.pos)):
            if cpos is not None:
                while self.pos[ii] != cpos + 1:
                    # if self.pos[ii] > cpos + 1:
                    #     print('gap at {}'.format(ii))
                    #     print(self.pos)
                    #     print(self.val)
                    #     break
                    # if self.pos[ii] < cpos + 1:
                    #     print('what!?')
                    #     print(self.pos)
                    #     print(self.val)
                    #     break
                    # missing some positions
                    new_pos.append(cpos + 1)
                    new_val.append(self.val[ii])
                    cpos += 1
            new_pos.append(self.pos[ii])
            new_val.append(self.val[ii])
            cpos = self.pos[ii]
        # update data
        self.pos = new_pos
        self.val = new_val


class EM:
    def __init__(self, npos, nstates_per_pos, values=['c1', 'c2']):
        self.prof = []
        self.npos = npos
        self.nstates = nstates_per_pos
        self.values = values
        for pp in range(self.npos):
            vals = []
            for ss in range(self.nstates):
                # init with a random value
                vals.append(random.sample(values, 1)[0])
            self.prof.append(vals)

    def print(self, digits=4):
        for pp in range(self.npos):
            tp = []
            for val in self.prof[pp]:
                tp.append('{}'.format(val))
            print('S{}:\t'.format(pp) + '  '.join(tp))

    def updateWithObs(self, obs):
        obs_names = list(obs.keys())
        random.shuffle(obs_names)
        # loop over observations and integrate them one by one
        for obsn in obs_names:
            ob = obs[obsn]
            # find the best profile
            bprof = ob.findBestEMProfile(self)
            # overwrite profile
            for obs_idx, pp in enumerate(ob.pos):
                self.prof[pp][bprof] = ob.val[obs_idx]

    def consensusWithObs(self, obs):
        # loop over observations, assign, and tally
        counts = {}
        for obn in obs:
            ob = obs[obn]
            # find the best profile
            bprof = ob.findBestEMProfile(self)
            # overwrite profile
            for obs_idx, pp in enumerate(ob.pos):
                oname = '{}_{}_{}'.format(pp, bprof, ob.val[obs_idx])
                if oname not in counts:
                    counts[oname] = 0
                counts[oname] += 1
        # pick the most supported value for each profile position
        ndiffs = 0
        for pp in range(self.npos):
            for ss in range(self.nstates):
                top_val = self.values[0]
                top_val_c = 0
                for val in self.values:
                    oname = '{}_{}_{}'.format(pp, ss, val)
                    if oname in counts and counts[oname] > top_val_c:
                        top_val_c = counts[oname]
                        top_val = val
                self.prof[pp][ss] = top_val
                # how many observations were different at this position
                ndiffs_c = 0
                for val in self.values:
                    oname = '{}_{}_{}'.format(pp, ss, val)
                    if oname in counts and val != top_val:
                        ndiffs_c += counts[oname]
                ndiffs += ndiffs_c
        return (ndiffs)
