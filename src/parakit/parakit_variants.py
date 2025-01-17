def findVariants(nodes, vedges, reads, nmarkers=10, pos_offset=0,
                 output_tsv='calls.tsv'):
    # look for read support for variant edges in vedges
    var_reads = {}
    var_reads_ref = {}
    for readn in reads.path:
        path = reads.path[readn]
        # save which variants edges are taken by this read
        for pos, nod in enumerate(path):
            if nod in vedges and pos + 1 < len(path) and \
               path[pos + 1] in vedges[nod]:
                # found a variant edge
                for ve in vedges[nod][path[pos + 1]]:
                    cvid = ve['cvid']
                    if ve['ref']:
                        # "non-pathogenic" allele at an interesting site
                        if cvid not in var_reads_ref:
                            var_reads_ref[cvid] = []
                        var_reads_ref[cvid].append(readn)
                    else:
                        # either the "pathogenic" allele
                        if cvid not in var_reads:
                            var_reads[cvid] = []
                        var_reads[cvid].append({'read': readn, 'pos': pos + 1})
    # Looking for isolate variants, e.g. from a short gene conversion event
    # Test if any covered node is surrounded by marker from module 2
    print('Looking for isolated known variants...')
    var_reads_f = {}
    # for each read, does it traverse each supported variant, yes/no
    read_sum = {}
    # check every variant
    for varid in var_reads:
        supp_reads = []
        # for each (read, position) check X markers around
        for rp in var_reads[varid]:
            c2_marks = 0
            c1_marks = 0
            read = reads.path[rp['read']]
            # check X markers downstream
            nmarks = 0
            pos = rp['pos'] + 1
            while nmarks < nmarkers / 2 and pos < len(read):
                if nodes[read[pos]]['class'] == 'c1':
                    c1_marks += 1
                    nmarks += 1
                elif nodes[read[pos]]['class'] == 'c2':
                    c2_marks += 1
                    nmarks += 1
                pos += 1
            # check X markers upstream
            nmarks = 0
            pos = rp['pos'] - 1
            while nmarks < nmarkers / 2 and pos >= 0:
                if nodes[read[pos]]['class'] == 'c1':
                    c1_marks += 1
                    nmarks += 1
                elif nodes[read[pos]]['class'] == 'c2':
                    c2_marks += 1
                    nmarks += 1
                pos += -1
            if c2_marks > 3 * c1_marks:
                supp_reads.append(rp['read'])
        # we keep cases with at least 3 supporting reads
        if len(supp_reads) >= 3:
            var_reads_f[varid] = supp_reads
            for read in supp_reads:
                if read not in read_sum:
                    read_sum[read] = {}
                read_sum[read][varid] = True

    # to record the two nodes associated with the alt edge of a variant
    # (for reporting later)
    var_node = {}
    # for fusions we'll only look for breakpoint after the ealiest clinvar
    # variant (minus 1kb say)
    min_pos_var = []
    for nod1 in vedges:
        for nod2 in vedges[nod1]:
            for ve in vedges[nod1][nod2]:
                if not ve['ref']:
                    var_node[ve['cvid']] = [nod1, nod2]
                min_pos_var.append(nodes[nod1]['rpos_min'])
    min_pos_var = min(min_pos_var) - 1000

    # Looking for deletion/fusion events
    # Test if any node in module/copy 2 is flanked by markers from
    # copy 1 upstream, and from copy 2 downstream
    print('Looking for deletion/fusion variants...')
    # for later: list candidate reads supporting a functional allele
    # because containing at least one module 2 "window"
    cand_func_reads = {}
    # we'll keep cases with at least 3 supporting reads
    # we don't want to "discover" the cycle going from module 1 to 2
    # hence we don't want fusion around the position of the cycle
    cyc_pos = 0
    for nod in nodes:
        if nodes[nod]['class'] == 'cyc_r':
            cyc_pos = nodes[nod]['rpos_min']
            break
    # to record fusion-supporting reads by node start
    # fus_reads['nod1'] = ['read1', 'read3']
    fus_reads = {}
    for readn in reads.path:
        # prepare a marker profile vector
        marks_n = []
        marks_c = []
        marks_p1 = []
        marks_p2 = []
        for ii, nod in enumerate(reads.path[readn]):
            if nodes[nod]['class'] == 'c1':
                marks_n.append(nod)
                marks_c.append('c1')
                # save the position in module 1 and 2
                if nodes[nod]['rpos_min'] != nodes[nod]['rpos_min']:
                    # use pos min and max if different
                    marks_p1.append(nodes[nod]['rpos_min'])
                    marks_p2.append(nodes[nod]['rpos_max'])
                else:
                    # if they're the same, we're missing the max
                    marks_p1.append(nodes[nod]['rpos_min'])
                    # look for the next node with a informative rpos_max
                    tpos = -1
                    jj = ii + 1
                    while jj < len(reads.path[readn]) and tpos == -1:
                        nodj = reads.path[readn][jj]
                        if nodes[nodj]['rpos_min'] != nodes[nodj]['rpos_max']:
                            tpos = nodes[nodj]['rpos_max']
                        jj += 1
                    marks_p2.append(tpos)
            elif nodes[nod]['class'] == 'c2':
                marks_n.append(nod)
                marks_c.append('c2')
                # save the position in module 1 and 2
                if nodes[nod]['rpos_min'] != nodes[nod]['rpos_min']:
                    # use pos min and max if different
                    marks_p1.append(nodes[nod]['rpos_min'])
                    marks_p2.append(nodes[nod]['rpos_max'])
                else:
                    # if they're the same, we're missing the min
                    marks_p2.append(nodes[nod]['rpos_max'])
                    # look for the next node with a informative rpos_max
                    tpos = -1
                    jj = ii - 1
                    while jj > 0 and tpos == -1:
                        nodj = reads.path[readn][jj]
                        if nodes[nodj]['rpos_min'] != nodes[nodj]['rpos_max']:
                            tpos = nodes[nodj]['rpos_min']
                        jj += -1
                    marks_p1.append(tpos)
        # skip if not enough markers for the test
        if len(marks_n) < nmarkers * 2:
            continue
        # compute a score to measure how much the flank are c1-c2 specific
        for ii in range(nmarkers, len(marks_n) - nmarkers):
            c1_n = marks_c[(ii-nmarkers+1):(ii+1)].count('c1')
            c2_n = marks_c[(ii+1):(ii+nmarkers+1)].count('c2')
            score = min(c2_n / nmarkers, c1_n / nmarkers)
            # save any position over the minimum score threshold
            if score > .8:
                # save info about the breakpoint
                fus_info = {'node_l': marks_n[ii],
                            'pos_l_1': marks_p1[ii],
                            'pos_l_2': marks_p2[ii],
                            'node_u': marks_n[ii+1],
                            'pos_u_1': marks_p1[ii+1],
                            'pos_u_2': marks_p2[ii+1],
                            'readn': readn}
                if marks_n[ii] not in fus_reads:
                    fus_reads[marks_n[ii]] = []
                fus_reads[marks_n[ii]].append(fus_info)
            if c2_n / nmarkers > .5 and c1_n / nmarkers < .5 and \
               nodes[marks_n[ii]]['rpos_min'] > min_pos_var:
                cand_func_reads[readn] = True
    # filter candidate fusions
    # keep nodes in range and one per read
    # (starting with the ones supported by the most reads)
    fus_reads_f = {}
    # sort nodes with some read support by number of reads
    fusnodes = sorted(list(fus_reads.keys()), key=lambda nod: -len(fus_reads[nod]))
    used_reads = {}
    for fusnod in fusnodes:
        # skip if node out of range
        if (nodes[fusnod]['rpos_min'] < min_pos_var) or \
           (nodes[fusnod]['rpos_min'] > cyc_pos - 100):
            continue
        # check that enough unused supporting reads
        fus_sup_reads = 0
        for fus_info in fus_reads[fusnod]:
            if fus_info['readn'] not in used_reads:
                fus_sup_reads += 1
        if fus_sup_reads >= 3:
            # save fusion
            # fusion name
            fus_n = 'FUS_{}'.format(fusnod)
            # lower and upper bounds of the confidence interval
            pos_l_1 = []
            pos_l_2 = []
            pos_u_1 = []
            pos_u_2 = []
            node_u = []
            for fus_info in fus_reads[fusnod]:
                if fus_info['pos_l_1'] != -1:
                    pos_l_1.append(fus_info['pos_l_1'])
                if fus_info['pos_l_2'] != -1:
                    pos_l_2.append(fus_info['pos_l_2'])
                if fus_info['pos_u_1'] != -1:
                    pos_u_1.append(fus_info['pos_u_1'])
                if fus_info['pos_u_2'] != -1:
                    pos_u_2.append(fus_info['pos_u_2'])
                node_u.append(fus_info['node_u'])
            pos_l_1 = min(pos_l_1)
            pos_l_2 = min(pos_l_2)
            pos_u_1 = max(pos_u_1)
            pos_u_2 = max(pos_u_2)
            node_u = max(node_u)
            # used to report fusions and other variants
            var_node[fus_n] = [fusnod, fusnod]
            # save info
            fus_reads_f[fus_n] = {'pos_l_1': pos_l_1,
                                  'pos_l_2': pos_l_2,
                                  'pos_u_1': pos_u_1,
                                  'pos_u_2': pos_u_2,
                                  'node_u': node_u}
            for fus_info in fus_reads[fusnod]:
                readn = fus_info['readn']
                # mark reads as "used"
                used_reads[readn] = True
                # remember that this read support this fusion
                if readn not in read_sum:
                    read_sum[readn] = {}
                read_sum[readn][fus_n] = True

    # looking for reads supporting a normal/functional allele
    # reads mostly from module 2 and with none of the
    # candidate variants identified above
    print('Looking for reads supporting a functional gene allele...')
    # list variants identified above (both clinvar and fusions)
    varids = {}
    for readn in read_sum:
        for varid in read_sum[readn]:
            varids[varid] = True
    # mark reads supporting the ref allele
    for varid in varids:
        if varid in var_reads_ref:
            for readn in var_reads_ref[varid]:
                # and were not supporting a pathogenic variant
                if readn not in read_sum:
                    cand_func_reads[readn] = True

    # print read support summary
    print("Covered variants and their supporting reads:\n")
    # order variants by position
    var_names = sorted(list(var_reads_f) + list(fus_reads_f),
                       key=lambda k: nodes[var_node[k][0]]['rpos_max'])
    # prepare TSV output
    for_tsv = ['\t'.join(['read', 'variant', 'node', 'allele', 'node_u',
                          'pos_l_1', 'pos_l_2', 'pos_u_1', 'pos_u_2'])]
    for_tsv_nas = '\t'.join(['NA']*5)
    printed_reads = {}
    for read in list(read_sum.keys()) + list(cand_func_reads.keys()):
        # don't write multiple time the same read
        if read in printed_reads:
            continue
        to_print = read + '\t'
        for varid in var_names:
            for_tsv_r = '\t'.join([read, varid, var_node[varid][1]])
            if read not in read_sum or varid not in read_sum[read]:
                # check if the read support a reference node for this variant
                if varid in var_reads_ref and read in var_reads_ref[varid]:
                    # if so, print ----
                    to_print += '-' * len(varid) + '\t'
                    for_tsv.append(for_tsv_r + '\tref\t' + for_tsv_nas)
                else:
                    # if not, print a gap
                    to_print += ' ' * len(varid) + '\t'
                    for_tsv.append(for_tsv_r + '\tNA\t' + for_tsv_nas)
            else:
                # the read goes through the variant, print the variant name
                to_print += varid + '\t'
                if varid in fus_reads_f:
                    # it's a fusion, let's also print breakpoint information
                    fus_info = fus_reads_f[varid]
                    fus_bkpt = [fus_info['node_u'],
                                fus_info['pos_l_1'] + pos_offset,
                                fus_info['pos_l_2'] + pos_offset,
                                fus_info['pos_u_1'] + pos_offset,
                                fus_info['pos_u_2'] + pos_offset]
                    fus_bkpt = '\t'.join([str(ii) for ii in fus_bkpt])
                    for_tsv.append(for_tsv_r + '\talt\t' + fus_bkpt)
                else:
                    for_tsv.append(for_tsv_r + '\talt\t' + for_tsv_nas)
        print(to_print)
        printed_reads[read] = True
    print()

    # write TSV output
    print('Writing summary in ' + output_tsv + ' TSV...')
    with open(output_tsv, 'wt') as outf:
        outf.write('\n'.join(for_tsv) + '\n')


def estimateCopyNumber(nodes, reads, window_size=20):
    fl_l = 0
    fl_r = 0
    cyc_l = 0
    cyc_r = 0
    # get flanking reference nodes
    fl_nodes = []
    for noden in nodes:
        if nodes[noden]['rpos_min'] == nodes[noden]['rpos_max']:
            fl_nodes.append(noden)
    # count the reads taking the flank or cycling edges
    for readn in reads.path:
        path = reads.path[readn]
        for pos, noden in enumerate(path):
            # increment counts for each informative edge
            if nodes[noden]['class'] == 'cyc_l':
                # check up to 'window_size' nodes upstream
                up_pos = max(0, pos - window_size)
                flank_found = False
                for fln in fl_nodes:
                    if fln in path[up_pos:pos]:
                        flank_found = True
                        break
                if flank_found:
                    # flank -> module start
                    fl_l += 1
                else:
                    # cycle -> module start
                    cyc_l += 1
            elif nodes[noden]['class'] == 'cyc_r':
                # check up to 'window_size' nodes downstream
                dw_pos = min(len(path), pos + window_size)
                flank_found = False
                for fln in fl_nodes:
                    if fln in path[pos:dw_pos]:
                        flank_found = True
                        break
                if flank_found:
                    # module end -> flank
                    fl_r += 1
                else:
                    # module end -> cycle
                    cyc_r += 1
    print('flank_left\t{}'.format(fl_l))
    print('flank_right\t{}'.format(fl_r))
    print('cycle_left\t{}'.format(cyc_l))
    print('cycle_right\t{}'.format(cyc_r))
    fl = (fl_l + fl_r) / 2
    print('flank_mean\t{}'.format(fl))
    cyc = (cyc_l + cyc_r) / 2
    print('cycle_mean\t{}'.format(cyc))
    cn = (fl + cyc) / fl
    # assume diploid genome
    cn *= 2
    print('cn\t{}'.format(round(cn, 4)))
