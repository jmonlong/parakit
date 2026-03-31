class ReadPosList:
    def __init__(self):
        # record reads and position(s)
        self.read_to_pos = {}
        # list of reads to exclude from future updates
        # (because they were removed earlier)
        self.exc_list = {}

    def addReadPos(self, read_name, read_pos):
        if read_name not in self.read_to_pos:
            self.read_to_pos[read_name] = []
        self.read_to_pos[read_name].append(read_pos)

    def addReadPosList(self, rpl):
        # add read positions
        for readn in rpl.read_to_pos:
            if readn not in self.read_to_pos:
                self.read_to_pos[readn] = []
            for read_pos in rpl.read_to_pos[readn]:
                self.read_to_pos[readn].append(read_pos)

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
    """Set of reads aligned to the pangenome graph

    Saves reads' alignment as path and where the alignment start/end
    in each node and in the read.

    Attributes:
        edge_to_readpos : dict edge name (n_m) -> ReadPosList object
        nsuc : dict counting the number of edges (e.g. nsuc[n][m])
        path : dict read name -> list of node IDs
        startpos : dict read name -> list of starting positions (on nodes)
        endpos : dict read name -> list of ending positions (on nodes)
        readpos : dict read name -> list of starting positions (in the read)

    Methods:
        addRead : add a new read alignment
        getPath :
        nReads :
        hasRead :
        getStartPos :
        getEndPos :
        getReadPos :
        getReadPosList :
        getAllReadPosList :
        predictLocalCopy :
        print :
        listSuccessors :
        overlapWithPath :
        overlapWithWalks :

    """
    def __init__(self):
        self.edge_to_readpos = {}
        self.nsuc = {}
        self.path = {}
        self.startpos = {}
        self.endpos = {}
        self.readpos = {}

    def addRead(self, read_name, path, startpos=[], endpos=[], readpos=[]):
        """Add a new read alignment

        If provided startpos/endpos/readpos should have the same
        length as path. They represent where the alignment
        started/ended in each node of the alignment path. For readpos,
        it records the (starting) position in the read for each node
        in the alignment path.

        Args:
            read_name : the name of the read
            path : list of nodes traversed by the read
            startpos : where do the alignment start in each node
            endpos : where do the alignment start in each node
            readpos : where do the alignment start in the read

        """
        # save path and position in the sequenced read
        self.path[read_name] = path
        if len(startpos) > 0:
            self.startpos[read_name] = startpos
        if len(endpos) > 0:
            self.endpos[read_name] = endpos
        if len(readpos) > 0:
            self.readpos[read_name] = readpos
        # save edges support
        for pos in range(len(path)-1):
            ename = path[pos] + '_' + path[pos+1]
            # add read-positions for this edge
            if ename not in self.edge_to_readpos:
                self.edge_to_readpos[ename] = ReadPosList()
            self.edge_to_readpos[ename].addReadPos(read_name, pos)
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

    def predictLocalCopy(self, read_name, var_pos, nodes, nmarkers):
        read = self.path[read_name]
        # check X markers downstream
        c2_marks_d = 0
        c1_marks_d = 0
        nmarks = 0
        pos = var_pos + 1
        while nmarks < nmarkers / 2 and pos < len(read):
            if nodes[read[pos]]['class'] == 'c1':
                c1_marks_d += 1
                nmarks += 1
            elif nodes[read[pos]]['class'] == 'c2':
                c2_marks_d += 1
                nmarks += 1
            pos += 1
        # check X markers upstream
        c2_marks_u = 0
        c1_marks_u = 0
        nmarks = 0
        pos = var_pos - 1
        while nmarks < nmarkers / 2 and pos >= 0:
            if nodes[read[pos]]['class'] == 'c1':
                c1_marks_u += 1
                nmarks += 1
            elif nodes[read[pos]]['class'] == 'c2':
                c2_marks_u += 1
                nmarks += 1
            pos += -1
        c1_marks = c1_marks_d + c1_marks_u
        c2_marks = c2_marks_d + c2_marks_u
        if c2_marks > 3 * c1_marks:
            return ('c2')
        elif c2_marks_d > 3 * c1_marks_d:
            return ('c2')
        elif c2_marks_u > 3 * c1_marks_u:
            return ('c2')
        elif c1_marks > 3 * c2_marks:
            return ('c1')
        else:
            return (None)

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
