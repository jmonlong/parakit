class Variant:
    def __init__(self, ref_trav=None, alt_trav=None,
                 ref_seq=None, alt_seq=None, varid=None):
        # list of reference traversals
        self.ref_trav = ref_trav
        self.alt_trav = alt_trav
        self.ref_seq = ref_seq
        self.alt_seq = alt_seq
        self.reads_ref = set()
        self.reads_alt = set()
        self.pos = 0
        self.end = None
        self.pos_error = None
        self.copy = None
        self.varid = varid
        self.clinvar = None

    def addRefRead(self, readn):
        self.reads_ref.add(readn)

    def addAltRead(self, readn):
        self.reads_alt.add(readn)

    def nAltReads(self):
        return (len(self.reads_alt))

    def toTsv(self, include_headers=False):
        res = []
        # potentially include the header
        headers = ['variant', 'pos', 'end', 'node', 'ref_trav',
                   'alt_trav', 'alt_seq', 'copy']
        if include_headers:
            res.append('\t'.join(headers))
        # prepare row for this variant
        fmt = '\t'.join(['{}'] * len(headers))
        pos = self.pos if self.pos != 0 else 'NA'
        end = self.end if self.end is not None else self.pos
        ref_trav = ['_'.join(trav) for trav in self.ref_trav]
        ref_trav = '-'.join(ref_trav)
        res_r = fmt.format(self.getVariantID(), pos, end, self.alt_trav[0],
                           ref_trav, '_'.join(self.alt_trav),
                           self.alt_seq, self.copy)
        res.append(res_r)
        return (res)

    def toReadTsv(self, include_headers=False):
        res = []
        # potentially include the header
        headers = ['read', 'variant', 'pos', 'end', 'node', 'allele', 'sig',
                   'copy', 'node_end', 'pos_ci']
        if include_headers:
            res.append('\t'.join(headers))
        # prepare one row per read
        fmt = '\t'.join(['{}'] * len(headers))
        n_start = self.alt_trav[0]
        pos = self.pos if self.pos != 0 else 'NA'
        end = self.end if self.end is not None else self.pos
        node_end = self.alt_trav[1] if len(self.alt_trav) == 2 else 'NA'
        pos_ci = self.pos_error if self.pos_error is not None else 'NA'
        for readn in self.reads_ref:
            res_r = fmt.format(readn, self.getVariantID(), pos, end, n_start,
                               'ref', self.clinvar, self.copy,
                               node_end, pos_ci)
            res.append(res_r)
        for readn in self.reads_alt:
            res_r = fmt.format(readn, self.getVariantID(), pos, end, n_start,
                               'alt', self.clinvar, self.copy,
                               node_end, pos_ci)
            res.append(res_r)
        return (res)

    def assignCopy(self, nodes):
        # is it a variant path relative to c1 or c2?
        c1_balance = 0
        for node in self.alt_trav:
            if nodes[node]['class'] == 'c1':
                c1_balance += 1
            if nodes[node]['class'] == 'c2':
                c1_balance += -1
        # assert c1_balance != 0, "No copy assignment for variant " + self.getVariantID()
        if c1_balance > 0:
            self.copy = 'c2'
        else:
            self.copy = 'c1'

    def findRefTraversal(self, nodes):
        ref_paths = [[self.alt_trav[0]]]
        self.ref_trav = []
        while len(ref_paths) > 0 and len(ref_paths) < 1000:
            ref_path = ref_paths.pop()
            if len(ref_path) > 20:
                continue
            for n_next in nodes[ref_path[-1]]['sucs']:
                if n_next == self.alt_trav[-1]:
                    self.ref_trav.append(ref_path + [n_next])
                elif nodes[n_next]['class'] == 'none' or (nodes[n_next]['class'] == 'c2' and self.copy == 'c2') or (nodes[n_next]['class'] == 'c1' and self.copy == 'c1'):
                    ref_paths.append(ref_path + [n_next])

    def getVariantID(self):
        if self.varid is None:
            return ('_'.join(self.alt_trav))
        else:
            return (self.varid)

    def findPosition(self, nodes):
        if self.copy is None:
            self.assignCopy(nodes)
        if self.copy == 'c1':
            self.pos = nodes[self.alt_trav[0]]['rpos_min']
        else:
            self.pos = nodes[self.alt_trav[0]]['rpos_max']
        self.pos += nodes[self.alt_trav[0]]['size']

    def findRefSequence(self, nodes):
        if self.ref_seq is not None:
            print("Warning: Overwriting reference sequence for variant " +
                  self.getVariantID())
        if self.ref_trav is not None:
            self.ref_seq = ''
            if len(self.ref_trav) > 0:
                for node in self.ref_trav[0][1:-1]:
                    self.ref_seq += nodes[node]['seq']

    def findAltSequence(self, nodes):
        if self.alt_seq is not None:
            print("Warning: Overwriting alternate sequence for variant " +
                  self.getVariantID())
        self.alt_seq = ''
        for node in self.alt_trav[1:-1]:
            self.alt_seq += nodes[node]['seq']

    def fillInfo(self, nodes):
        self.assignCopy(nodes)
        self.findRefTraversal(nodes)
        self.findPosition(nodes)
        self.findAltSequence(nodes)
        self.findRefSequence(nodes)

    def fillFusionInfo(self, nodes):
        # find the start and end position of a fusion
        # based on the alt traversal and copy
        assert self.alt_trav is not None and self.copy is not None, \
            'ALT traversal and copy information missing for fusion'
        # extract position on module 1 and 2 for both boundary nodes
        u1 = nodes[self.alt_trav[0]]['rpos_min']
        d1 = nodes[self.alt_trav[1]]['rpos_min']
        u2 = nodes[self.alt_trav[0]]['rpos_max']
        d2 = nodes[self.alt_trav[1]]['rpos_max']
        # assign pos as the starting position, so it depends on configuration
        if self.copy == 'c1':
            # fusion c1-c2
            self.pos = round((u1 + d1) / 2)
            self.end = round((u2 + d2) / 2)
        elif self.copy == 'c2':
            # fusion c2-c1
            self.pos = round((u2 + d2) / 2)
            self.end = round((u1 + d1) / 2)
        # error in position is always the same no matter the configuration
        self.pos_error = round((abs(u1 - d1) + abs(u2 - d2)) / 4)


class Fusions:
    def __init__(self, nmarkers=10):
        # potential fusion variants (varid -> Variant object)
        self.variants = {}
        # number of markers used to decide if a variant is (likely) converted
        self.nmarkers = nmarkers

    def importReads(self, reads, nodes, support_only=False):
        # find cycling nodes to split into subreads
        cyc_nodes = set()
        for node in nodes:
            if 'cyc' in nodes[node]['class']:
                cyc_nodes.add(node)
        # import each read
        for readn in reads.path:
            path = reads.path[readn]
            # extract marker sequence for each subread
            subread_to_markers = {}
            cur_markers = []
            cur_nodes = []
            for node in path:
                if node in cyc_nodes:
                    if len(cur_nodes) > 2 * self.nmarkers:
                        # save current subread
                        sreadn = readn + '_' + str(len(subread_to_markers))
                        subread_to_markers[sreadn] = {
                            'markers': cur_markers,
                            'nodes': cur_nodes
                        }
                    cur_markers = []
                    cur_nodes = []
                elif nodes[node]['class'] in ['c1', 'c2']:
                    # this is a marker to record
                    cur_markers.append(nodes[node]['class'])
                    cur_nodes.append(node)
            # import each subread
            for sreadn in subread_to_markers:
                self.importSubread(subread_to_markers[sreadn]['markers'],
                                   subread_to_markers[sreadn]['nodes'],
                                   sreadn, support_only=support_only)

    def importSubread(self, markers, nodes, sreadn, support_only=False):
        u_c1 = 0
        u_c2 = 0
        d_c1 = markers.count('c1')
        d_c2 = markers.count('c2')
        # test for a module switch across the markers by subread
        for idx, marker in enumerate(markers):
            # update the counts
            if marker == 'c1':
                u_c1 += 1
                d_c1 += -1
            if marker == 'c2':
                u_c2 += 1
                d_c2 += -1
            # don't test switch points too close to the ends of the subread
            if idx < self.nmarkers or idx > len(markers) - self.nmarkers:
                continue
            # skip if no switch between this marker and the next
            # except if we are looking for supporting reads
            if marker == markers[idx + 1] and not support_only:
                continue
            # check proportion of c1/c2 upstream vs downstream
            fusion_start_copy = None
            fusion_end_copy = None
            if u_c1 > 4 * u_c2:
                fusion_start_copy = 'c1'
            elif u_c2 > 4 * u_c1:
                fusion_start_copy = 'c2'
            if d_c1 > 4 * d_c2:
                fusion_end_copy = 'c1'
            elif d_c2 > 4 * d_c1:
                fusion_end_copy = 'c2'
            fusid = 'fus_{}_{}'.format(fusion_start_copy, nodes[idx])
            if marker != fusion_start_copy and not support_only:
                continue
            fusion_signal = fusion_start_copy is not None and \
                fusion_end_copy is not None and \
                fusion_start_copy != fusion_end_copy
            if support_only:
                # annotate support only
                if fusid in self.variants:
                    if fusion_signal:
                        self.variants[fusid].addAltRead(sreadn)
                    elif self.variants[fusid].copy == fusion_start_copy:
                        self.variants[fusid].addRefRead(sreadn)
                # also annotate reference reads on the other copy?
                other_c = 'c1' if fusion_start_copy == 'c2' else 'c2'
                fusid2 = 'fus_{}_{}'.format(other_c, nodes[idx])
                if fusid2 in self.variants and not fusion_signal:
                    self.variants[fusid2].addRefRead(sreadn)
            if not support_only and fusion_signal:
                # looking for new fusion, potentially create a variant
                if fusid not in self.variants:
                    self.variants[fusid] = Variant(alt_trav=[nodes[idx],
                                                             nodes[idx+1]],
                                                   varid=fusid)
                    self.variants[fusid].copy = fusion_start_copy
                self.variants[fusid].addAltRead(sreadn)

    def filterVariants(self):
        selected_vars = {}
        # sort by the most supporting reads
        vars_s = sorted(list(self.variants.keys()),
                        key=lambda vv: -self.variants[vv].nAltReads())
        used_reads = set()
        # assign some (sub)reads
        for varid in vars_s:
            supp_reads = set()
            for readn in self.variants[varid].reads_alt:
                if readn not in used_reads:
                    supp_reads.add(readn)
            if len(supp_reads) >= 3:
                selected_vars[varid] = self.variants[varid]
                selected_vars[varid].reads_alt = supp_reads
                for readn in supp_reads:
                    used_reads.add(readn)
        self.variants = selected_vars

    def fillInfo(self, nodes):
        for varid in self.variants:
            self.variants[varid].fillFusionInfo(nodes)


class ConvertedVariants:
    def __init__(self, nmarkers=10):
        # potential converted variants (varid -> Variant object)
        self.variants = {}
        # number of markers used to decide if a variant is (likely) converted
        self.nmarkers = nmarkers
        # index the position/nodes of the current variant list
        self.pos_to_varids = {}
        self.node_to_varids = {}

    def addVariant(self, var):
        varid = var.getVariantID()
        self.variants[varid] = var
        # update the pos/node index
        if var.pos not in self.pos_to_varids:
            self.pos_to_varids[var.pos] = []
        self.pos_to_varids[var.pos].append(varid)
        if var.alt_trav[0] not in self.node_to_varids:
            self.node_to_varids[var.alt_trav[0]] = []
        self.node_to_varids[var.alt_trav[0]].append(varid)

    def decomposePangenome(self, nodes):
        # check for variant starting at any node
        for n_start in nodes:
            # well not any node, only if it's on the reference path
            if nodes[n_start]['class'] != 'none':
                continue
            # try to traverse from that node, looking for variant paths
            # as long as the assigned reference node (rnode) is the start node
            cand_paths = [[n_start]]
            var_paths = []
            while len(cand_paths) > 0:
                cur_path = cand_paths.pop()
                # try to extend the path
                for n_next in nodes[cur_path[-1]]['sucs']:
                    if nodes[n_next]['rnode'] == n_start:
                        cand_paths.append(cur_path + [n_next])
                        continue
                    if nodes[n_next]['class'] != 'none' or len(cur_path) == 1:
                        continue
                    # it reached the reference path again and
                    # traversed something before
                    # check if it's close enough from the starting position
                    start_dist = min(abs(nodes[n_next]['rpos_min'] -
                                         nodes[n_start]['rpos_min'] -
                                         nodes[n_start]['size']),
                                     abs(nodes[n_next]['rpos_max'] -
                                         nodes[n_start]['rpos_max'] -
                                         nodes[n_start]['size']))
                    if start_dist < 50:
                        var_paths.append(cur_path + [n_next])
            # create a Variant object from the variant paths
            for var_path in var_paths:
                var = Variant(alt_trav=var_path)
                var.fillInfo(nodes)
                if var.ref_trav is not None:
                    self.addVariant(var)

    def matchAnnotation(self, filen, offset):
        # read the input tsv file
        inf = open(filen, 'rt')
        heads = next(inf).rstrip().split('\t')
        n_matched = 0
        for line in inf:
            # parse the line and get some info
            line = line.rstrip().split('\t')
            pos = int(line[heads.index('start')]) - offset
            ref = line[heads.index('ref')]
            alt = line[heads.index('alt')]
            cvid = '{}_{}_{}'.format(line[heads.index('id')],
                                     line[heads.index('nuc.change')],
                                     line[heads.index('prot.change')])
            # if padding, remove it and update position
            if ref[0] == alt[0]:
                ref = '' if len(ref) == 1 else ref[1:]
                alt = '' if len(alt) == 1 else alt[1:]
                pos += 1
            # look for a variant starting at that position and
            # with matching alt sequences
            if pos not in self.pos_to_varids:
                # print('W: {} not matched (no variant at {})'.format(cvid, pos))
                continue
            for varid in self.pos_to_varids[pos]:
                if self.variants[varid].alt_seq == alt:
                    self.variants[varid].clinvar = cvid
                    n_matched += 1
            continue
            # print('W: {} not matched (variants at {} but different alt'
            #       ' sequence)'.format(cvid, pos))
        inf.close()
        print('{} annotated variants matched.'.format(n_matched))

    def importReads(self, reads, nodes):
        for readn in reads.path:
            path = reads.path[readn]
            # save which variants edges are taken by this read
            for pos, node in enumerate(path):
                if pos == len(path) - 1:
                    # skip if there is no next node
                    # (the edge is used to match the allele)
                    continue
                next_node = path[pos + 1]
                if node not in self.node_to_varids:
                    # no variant starts at this node, skip
                    continue
                # check if this region could be converted
                reg_copy = reads.predictLocalCopy(readn, pos + 1,
                                                  nodes,
                                                  self.nmarkers)
                if reg_copy is None:
                    # we somehow can't tell if we're in c1 or c2
                    continue
                # check if any variant matches
                for varid in self.node_to_varids[node]:
                    var = self.variants[varid]
                    # check that it's the appropriate conversion
                    if reg_copy != var.copy:
                        continue
                    # assign to alt allele if first edge matches
                    if var.alt_trav[1] == next_node:
                        var.addAltRead(readn)
                    # assign to ref allele if first edge matches
                    for ref_trav in var.ref_trav:
                        if ref_trav[1] == next_node:
                            var.addRefRead(readn)

    def filterVariants(self):
        selected_vars = {}
        for varid in self.variants:
            if self.variants[varid].nAltReads() >= 3:
                selected_vars[varid] = self.variants[varid]
        self.variants = selected_vars


def findVariants(nodes, annot_fn, reads, nmarkers=10, pos_offset=0,
                 output_tsv='calls.tsv'):

    decomposePangenome(nodes, pos_offset)
    
    # look for read support for variant edges in vedges
    vars = ConvertedVariants(nmarkers)
    vars.decomposePangenome(nodes)
    vars.matchAnnotation(annot_fn, pos_offset)
    vars.importReads(reads, nodes)
    vars.filterVariants()

    # look for deletions/fusions
    print('Looking for deletion/fusion variants...')
    fusions = Fusions()
    fusions.importReads(reads, nodes)
    fusions.importReads(reads, nodes, support_only=True)
    fusions.filterVariants()
    fusions.fillInfo(nodes)

    # merge variants and print
    all_vars = vars.variants
    for varid in fusions.variants:
        all_vars[varid] = fusions.variants[varid]

    # sort them by position
    var_ids = sorted(list(all_vars), key=lambda vid: all_vars[vid].pos)
    inc_headers = True
    for_tsv = []
    for vid in var_ids:
        # print or prepare the tsv output for this variant
        for_tsv += all_vars[vid].toReadTsv(inc_headers)
        inc_headers = False

    # write TSV output
    print('Writing summary in ' + output_tsv + ' TSV...')
    with open(output_tsv, 'wt') as outf:
        outf.write('\n'.join(for_tsv) + '\n')


def decomposePangenome(nodes, pos_offset=0, output_tsv='variants.tsv'):
    # look for read support for variant edges in vedges
    vars = ConvertedVariants()
    vars.decomposePangenome(nodes)
    vars = vars.variants

    # sort them by position
    var_ids = sorted(list(vars), key=lambda vid: vars[vid].pos)
    inc_headers = True
    for_tsv = []
    for vid in var_ids:
        # print or prepare the tsv output for this variant
        for_tsv += vars[vid].toTsv(inc_headers)
        inc_headers = False

    # write TSV output
    print('Writing summary in ' + output_tsv + ' TSV...')
    with open(output_tsv, 'wt') as outf:
        outf.write('\n'.join(for_tsv) + '\n')


def estimateCopyNumber(nodes, reads, window_size=20):
    fl = 0
    cyc = 0
    # get flanking reference nodes
    fl_nodes = []
    for noden in nodes:
        # if nodes[noden]['rpos_min'] == nodes[noden]['rpos_max']:
        #     fl_nodes.append(noden)
        if nodes[noden]['class'] == 'ref':
            fl_nodes.append(noden)
    # count the reads taking the flank or cycling edges
    for readn in reads.path:
        path = reads.path[readn]
        for pos, noden in enumerate(path):
            # increment counts for each informative edge
            if nodes[noden]['class'] in ['cyc_l', 'cyc_r']:
                # check up to 'window_size' nodes upstream and downstream
                up_pos = max(0, pos - window_size)
                dw_pos = min(len(path), pos + window_size)
                flank_found = False
                for fln in fl_nodes:
                    if fln in path[up_pos:dw_pos]:
                        flank_found = True
                        break
                if flank_found:
                    # flank -> module start
                    fl += 1
                else:
                    # cycle -> module start
                    cyc += 1
    print('flank_total\t{}'.format(fl))
    fl = fl / 2
    print('flank\t{}'.format(fl))
    print('cycle_total\t{}'.format(cyc))
    cyc = cyc / 2
    print('cycle\t{}'.format(cyc))
    cn = (fl + cyc) / fl
    # assume diploid genome
    cn *= 2
    print('cn\t{}'.format(round(cn, 4)))
