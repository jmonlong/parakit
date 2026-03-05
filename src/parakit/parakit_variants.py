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

    def getEnd(self):
        if self.end is not None:
            return (self.end)
        else:
            return (self.pos)

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
        end = self.getEnd()
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
        headers = ['variant', 'pos', 'end', 'node', 'sig', 'copy',
                   'fusion_start_node', 'pos_ci', 'supp_reads',
                   'reads_alt', 'reads_ref']
        if include_headers:
            res.append('\t'.join(headers))
        # prepare one row per read
        fmt = '\t'.join(['{}'] * len(headers))
        pos = self.pos if self.pos != 0 else 'NA'
        end = self.getEnd()
        node = self.alt_trav[1] if len(self.alt_trav) > 1 else self.alt_trav[0]
        node_fstart = self.alt_trav[0] if self.pos_error is not None else 'NA'
        pos_ci = self.pos_error if self.pos_error is not None else 'NA'
        reads_ref = ','.join(list(self.reads_ref))
        if reads_ref == '':
            reads_ref = 'NA'
        reads_alt = ','.join(list(self.reads_alt))
        if reads_alt == '':
            reads_alt = 'NA'
        res_r = fmt.format(self.getVariantID(), pos, end, node,
                           self.clinvar, self.copy, node_fstart, pos_ci,
                           len(self.reads_alt), reads_alt, reads_ref)
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

    def offsetPosition(self, offset):
        self.pos += offset
        if self.end is not None:
            self.end += offset


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
                        sreadn = readn + '_sr' + str(len(subread_to_markers))
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

    def offsetPositions(self, offset):
        for varid in self.variants:
            self.variants[varid].offsetPosition(offset)


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
                    if reg_copy != self.variants[varid].copy:
                        continue
                    self.checkReadSupport(varid, path[pos:], readn)

    def checkReadSupport(self, varid, read, readn):
        var = self.variants[varid]
        # assign to alt allele if first edge matches
        alt_support = True
        node_ii = 0
        while node_ii < len(read) and node_ii < len(var.alt_trav):
            if var.alt_trav[node_ii] != read[node_ii]:
                alt_support = False
                break
            node_ii += 1
        # assign to ref allele if first edge matches
        ref_support = False
        for ref_trav in var.ref_trav:
            ref_support_c = True
            node_ii = 0
            while node_ii < len(read) and node_ii < len(ref_trav):
                if ref_trav[node_ii] != read[node_ii]:
                    ref_support_c = False
                    break
                node_ii += 1
            if ref_support_c:
                ref_support = True
                break
        # assign if unambiguous
        if alt_support and not ref_support:
            var.addAltRead(readn)
        elif ref_support and not alt_support:
            var.addRefRead(readn)
    
    def filterVariants(self):
        selected_vars = {}
        for varid in self.variants:
            if self.variants[varid].nAltReads() >= 3:
                selected_vars[varid] = self.variants[varid]
        self.variants = selected_vars

    def offsetPositions(self, offset):
        for varid in self.variants:
            self.variants[varid].offsetPosition(offset)


def findVariants(nodes, annot_fn, reads, nmarkers=10, pos_offset=0,
                 output_tsv='calls.tsv'):
    # decomposePangenome(nodes, pos_offset)
    # look for read support for variant edges in vedges
    vars = ConvertedVariants(nmarkers)
    vars.decomposePangenome(nodes)
    vars.matchAnnotation(annot_fn, pos_offset)
    vars.importReads(reads, nodes)
    vars.filterVariants()
    vars.offsetPositions(pos_offset)

    # look for deletions/fusions
    print('Looking for deletion/fusion variants...')
    fusions = Fusions()
    fusions.importReads(reads, nodes)
    fusions.importReads(reads, nodes, support_only=True)
    fusions.filterVariants()
    fusions.fillInfo(nodes)
    fusions.offsetPositions(pos_offset)

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


def readVariantCalls(filen):
    variants = {}
    inf = open(filen, 'rt')
    heads = next(inf).rstrip().split('\t')
    for line in inf:
        line = line.rstrip().split('\t')
        varid = line[heads.index('variant')]
        if varid not in variants:
            var = Variant(varid=varid)
            var.pos = int(line[heads.index('pos')])
            var.end = int(line[heads.index('end')])
            pos_ci = line[heads.index('pos_ci')]
            if pos_ci != 'NA':
                var.pos_error = int(pos_ci)
            sig = line[heads.index('sig')]
            if sig != 'None':
                var.clinvar = sig
            var.alt_trav = []
            node_fstart = line[heads.index('fusion_start_node')]
            if node_fstart != 'NA':
                var.alt_trav.append(node_fstart)
            var.alt_trav.append(line[heads.index('node')])
            var.copy = line[heads.index('copy')]
            variants[varid] = var
        # update read support
        reads_ref = line[heads.index('reads_ref')]
        if reads_ref != 'NA':
            for readn in reads_ref.split(','):
                variants[varid].addRefRead(readn)
        reads_alt = line[heads.index('reads_alt')]
        if reads_alt != 'NA':
            for readn in reads_alt.split(','):
                variants[varid].addAltRead(readn)
    inf.close()
    print('Read {} variants.'.format(len(variants)))
    return (variants)


def filterVariants(nodes, calls_fn, module=None, annotated_only=False,
                   fusion_only=False,
                   start_pos=None, end_pos=None,
                   output_tsv='calls.filtered.tsv'):
    variants = readVariantCalls(calls_fn)
    sel_ids = []
    for varid in variants:
        var = variants[varid]
        # filter by module
        if module is not None and var.copy != 'c' + module:
            continue
        # filter everything except annotated or fusions
        if annotated_only and fusion_only:
            if var.clinvar is None and var.pos_error is None:
                continue
        else:
            # filter by annotated
            if annotated_only and var.clinvar is None:
                continue
            # filter by fusion
            if fusion_only and var.pos_error is None:
                continue
        # filter by position
        if start_pos is not None and end_pos is not None:
            if var.pos > end_pos or var.getEnd() < start_pos:
                continue
        # keep this variant
        sel_ids.append(varid)
    # sort them by position
    var_ids = sorted(sel_ids, key=lambda vid: variants[vid].pos)
    inc_headers = True
    for_tsv = []
    for vid in var_ids:
        # print or prepare the tsv output for this variant
        for_tsv += variants[vid].toReadTsv(inc_headers)
        inc_headers = False

    # write TSV output
    with open(output_tsv, 'wt') as outf:
        outf.write('\n'.join(for_tsv) + '\n')


def decomposePangenome(nodes, pos_offset=0, output_tsv='variants.tsv',
                       annot_fn=None):
    # look for read support for variant edges in vedges
    vars = ConvertedVariants()
    vars.decomposePangenome(nodes)
    if annot_fn is not None:
        vars.matchAnnotation(annot_fn, pos_offset)
    vars.offsetPositions(pos_offset)
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


class NodeCoverage:
    def __init__(self, nodeid, size):
        self.nodeid = nodeid
        self.values = [0]
        self.lengths = [size]
        self.bins_start = []
        self.bins_end = []
        self.bins_cov = []

    def addRead(self, start, end):
        n_values = []
        n_lengths = []
        # find the first RLE block overlapping the start position
        ii = 0
        cpos = 0
        while ii < len(self.lengths):
            if start < cpos + self.lengths[ii]:
                # input range starts within this block
                if start != cpos:
                    # it doesn't starts exactly at the beginning of the block
                    # cut off the first part, without changing its value
                    n_values.append(self.values[ii])
                    n_lengths.append(start - cpos)
                    assert cpos < start
                    self.lengths[ii] += cpos - start
                break
            else:
                # we're before the range, add this block unchanged
                cpos += self.lengths[ii]
                n_values.append(self.values[ii])
                n_lengths.append(self.lengths[ii])
                ii += 1
        # move to next block while decrementing the range
        range_size = end - start
        while range_size > 0:
            if self.lengths[ii] <= range_size:
                # this block is completely within the range, just increment
                n_values.append(self.values[ii] + 1)
                n_lengths.append(self.lengths[ii])
                # decrease the range size
                range_size += -self.lengths[ii]
                ii += 1
            else:
                # the range ends within this block,
                # split in two and increment first one
                n_values.append(self.values[ii] + 1)
                n_lengths.append(range_size)
                self.lengths[ii] += -range_size
                range_size = 0
        # add last blocks, unchanged
        while ii < len(self.lengths):
            n_values.append(self.values[ii])
            n_lengths.append(self.lengths[ii])
            ii += 1
        self.values = n_values
        self.lengths = n_lengths

    def addRead2(self, start, end):
        cpos = 0
        n_values = []
        n_lengths = []
        ii = 0
        start_done = False
        end_done = False
        while ii < len(self.lengths):
            if start < cpos + self.lengths[ii] and not start_done:
                # overlaps the current block
                # cut into a first block like the current one
                n_values.append(self.values[ii])
                n_lengths.append(start - cpos)
                if start - cpos < 0:
                    print('start - cpos')
                    print('\t', self.lengths[ii], cpos, start, end)
                # potentially create the second block with incremented coverage
                if end > cpos + self.lengths[ii]:
                    n_values.append(self.values[ii] + 1)
                    n_lengths.append(self.lengths[ii] + cpos - start)
                    if self.lengths[ii] + cpos - start < 0:
                        print('self.lengths[ii] + cpos - start')
                        print('\t', self.lengths[ii], cpos, start, end)
                else:
                    # or the it's completely included within the current block
                    n_values.append(self.values[ii] + 1)
                    n_lengths.append(end - start)
                    if end - start < 0:
                        print('end - start')
                        print('\t', self.lengths[ii], cpos, start, end)
                    n_values.append(self.values[ii])
                    n_lengths.append(self.lengths[ii] + cpos - end)
                    if self.lengths[ii] + cpos - end < 0:
                        print('self.lengths[ii] + cpos - end')
                        print('\t', self.lengths[ii], cpos, start, end)
                start_done = True
            elif end < cpos + self.lengths[ii] and not end_done:
                # split current block in 2
                n_values.append(self.values[ii] + 1)
                n_lengths.append(end - cpos)
                if end - cpos < 0:
                    print('end - cpos')
                    print('\t', self.lengths[ii], cpos, start, end)
                n_values.append(self.values[ii])
                n_lengths.append(self.lengths[ii] + cpos - end)
                if self.lengths[ii] + cpos - end < 0:
                    print('self.lengths[ii] + cpos - end')
                    print('\t', self.lengths[ii], cpos, start, end)
                end_done = True
            else:
                # nothing to do
                if start_done and not end_done:
                    n_values.append(self.values[ii] + 1)
                else:
                    n_values.append(self.values[ii])
                n_lengths.append(self.lengths[ii])
                if self.lengths[ii] < 0:
                    print('self.lengths[ii]')
                    print('\t', self.lengths[ii], cpos, start, end)
            # update pointers
            cpos += self.lengths[ii]
            ii += 1
        self.values = n_values
        self.lengths = n_lengths

    def binCoverage(self, bin_size=100):
        self.bins_start = []
        self.bins_end = []
        self.bins_cov = []
        ii = 0
        cur_sum = 0
        cur_bin = bin_size
        offset = 0
        bin_start = 0
        while ii < len(self.lengths):
            if self.lengths[ii] - offset < cur_bin:
                cur_bin += -(self.lengths[ii] - offset)
                cur_sum += (self.lengths[ii] - offset) * self.values[ii]
                offset = 0
                ii += 1
            else:
                # save that bin
                cur_sum += cur_bin * self.values[ii]
                self.bins_cov.append(cur_sum / bin_size)
                self.bins_start.append(bin_start)
                self.bins_end.append(bin_start + bin_size)
                bin_start += bin_size
                # start a new one
                offset += cur_bin
                cur_bin = bin_size
                cur_sum = 0
        # last bin
        if cur_bin != bin_size:
            self.bins_cov.append(cur_sum / (bin_size - cur_bin))
            self.bins_start.append(bin_start)
            self.bins_end.append(bin_start + bin_size - cur_bin)

    def print(self):
        top = ''
        for ii in range(len(self.lengths)):
            top += str(self.values[ii]) * self.lengths[ii]
        print(top)

    def printCoverage(self):
        print('node\tstart\tend\tcoverage')
        for ii in range(len(self.bins_cov)):
            print('{}\t{}\t{}\t{}'.format(self.nodeid, self.bins_start[ii],
                                          self.bins_end[ii],
                                          self.bins_cov[ii]))
    

def estimateCopyNumberFromFlanks(nodes, reads, window_size=20):
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
    fl = fl / 2
    print('flank\t{}'.format(fl))
    cyc = cyc / 2
    print('cycle\t{}'.format(cyc))
    cn = (fl + cyc) / fl
    # assume diploid genome
    cn *= 2
    print('cn\t{}'.format(round(cn, 4)))


def mode(x, min_value=0, bin_size_pct=10):
    x.sort()
    x_mean = statistics.mean(x)
    # compute the number of element within 1% of the mean for each position
    diff_th = x_mean * bin_size_pct / 100
    bin_max_ii = None
    bin_max_n = 0
    for ii in range(len(x)):
        if x[ii] < min_value:
            continue
        jj = 0
        while ii+jj < len(x) and abs(x[ii] - x[ii+jj]) < diff_th:
            jj += 1
        if jj > bin_max_n:
            bin_max_n = jj
            bin_max_ii = ii
    # return the median of the best 'bin'
    return (statistics.median(x[bin_max_ii:(bin_max_ii+bin_max_n)]))


def estimateCopyNumber(nodes, reads, out_tsv_prefix='parakit.cn', window_size=20,
                       bin_size=50):
    # estimate from the flanks
    estimateCopyNumberFromFlanks(nodes, reads, window_size=window_size)
    # compute coverage on reference nodes
    node_cov = {}
    # also compute the edge coverage for reference nodes
    edge_cov = {}
    for noden in nodes:
        if nodes[noden]['class'] in ['ref', 'none']:
            node_cov[noden] = NodeCoverage(noden, nodes[noden]['size'])
            edge_cov[noden] = {'in': {}, 'out': {}}
    # process all reads
    for readn in reads.path:
        path = reads.path[readn]
        for nii, noden in enumerate(path):
            if noden in node_cov:
                # update node coverage
                rstart = reads.getStartPos(readn, nii)
                rend = reads.getEndPos(readn, nii)
                # print(readn)
                node_cov[noden].addRead(min(rstart, rend), max(rstart, rend))
                # update edge coverage
                if nii > 0:
                    enode = path[nii-1]
                    if enode not in edge_cov[noden]['in']:
                        edge_cov[noden]['in'][enode] = 0
                    edge_cov[noden]['in'][enode] += 1
                if nii + 1 < len(path):
                    enode = path[nii+1]
                    if enode not in edge_cov[noden]['out']:
                        edge_cov[noden]['out'][enode] = 0
                    edge_cov[noden]['out'][enode] += 1
    # debug
    # write output files
    fmt = '{}\t{}\t{}\t{}\n'
    outf = open(out_tsv_prefix + '.node.tsv', 'wt')
    outf.write(fmt.format('node', 'start', 'end', 'coverage'))
    # also list bin coverage in ref and none nodes
    flank_cov = []
    mod_cov = []
    for noden in node_cov:
        # bin node
        node_cov[noden].binCoverage(bin_size)
        for ii in range(len(node_cov[noden].bins_cov)):
            outf.write(fmt.format(noden, node_cov[noden].bins_start[ii],
                                  node_cov[noden].bins_end[ii],
                                  node_cov[noden].bins_cov[ii]))
        if nodes[noden]['class'] == 'ref':
            flank_cov += node_cov[noden].bins_cov
        elif nodes[noden]['class'] == 'none':
            mod_cov += node_cov[noden].bins_cov
    outf.close()
    outf = open(out_tsv_prefix + '.edge.tsv', 'wt')
    outf.write(fmt.format('node', 'enode', 'side', 'coverage'))
    for noden in edge_cov:
        for enode in edge_cov[noden]['in']:
            outf.write(fmt.format(noden, enode, 'in',
                                  edge_cov[noden]['in'][enode]))
        for enode in edge_cov[noden]['out']:
            outf.write(fmt.format(noden, enode, 'out',
                                  edge_cov[noden]['out'][enode]))
    outf.close()
    # compute the mode across all bins in 'ref' nodes
    flank_cov = mode(flank_cov, min_value=3)
    # normalize the average bin coverage in 'none' nodes within the module
    mod_cov = mode(mod_cov, min_value=3)
    # compute module copy estimate
    print('flank\t{}'.format(round(flank_cov, 4)))
    print('module\t{}'.format(round(mod_cov, 4)))
    cn = 2 * mod_cov / flank_cov
    print('cn\t{}'.format(round(cn, 4)))
