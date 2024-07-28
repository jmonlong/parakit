def findVariants(nodes, vedges, reads, nmarkers=10, pos_offset=0):
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

    # to record the node associate with a variant (for reporting later)
    var_node = {}
    # for fusions we'll only look for breakpoint after the ealiest clinvar
    # variant (minus 1kb say)
    min_pos_var = []
    for nod1 in vedges:
        for nod2 in vedges[nod1]:
            for ve in vedges[nod1][nod2]:
                var_node[ve['cvid']] = nod2
                min_pos_var.append(nodes[nod2]['rpos_max'])
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
            cyc_pos = nodes[nod]['rpos_max']
            break
    # to record fusion-supporting reads by node start
    # fus_reads['nod1'] = ['read1', 'read3']
    fus_reads = {}
    for readn in reads.path:
        # prepare a marker profile vector
        marks = []
        marks_c = []
        for nod in reads.path[readn]:
            if nodes[nod]['class'] == 'c1':
                marks.append(nod)
                marks_c.append('c1')
            elif nodes[nod]['class'] == 'c2':
                marks.append(nod)
                marks_c.append('c2')
        # skip if not enough markers for the test
        if len(marks) < nmarkers * 2:
            continue
        # compute a score to measure how much the flank are c1-c2 specific
        fus_score_max = 0
        fus_score_pos = ''
        for ii in range(nmarkers, len(marks) - nmarkers):
            c1_n = marks_c[(ii-nmarkers):ii].count('c1')
            c2_n = marks_c[(ii+1):(ii+nmarkers+1)].count('c2')
            score = min(c2_n / nmarkers, c1_n / nmarkers)
            if score > fus_score_max:
                fus_score_max = score
                fus_score_pos = ii
            if c2_n / nmarkers > .5 and c1_n / nmarkers < .5 and \
               nodes[marks[ii]]['rpos_max'] > min_pos_var:
                cand_func_reads[readn] = True
        # if best score is high, save fusion
        if fus_score_max > .5:
            if marks[fus_score_pos] not in fus_reads:
                fus_reads[marks[fus_score_pos]] = []
            fus_reads[marks[fus_score_pos]].append(readn)
    # filter candidate fusions
    fus_reads_f = {}
    for fusnod in fus_reads:
        if (len(fus_reads[fusnod]) >= 3) and \
           (nodes[fusnod]['rpos_max'] > min_pos_var) and \
           (nodes[fusnod]['rpos_max'] < cyc_pos - 100):
            fus_n = '{}_FUS_{}'.format(fusnod,
                                       pos_offset + nodes[fusnod]['rpos_max'])
            # anno['node'][fus_n] = fus
            var_node[fus_n] = fusnod
            fus_reads_f[fus_n] = fus_reads[fusnod]
            for readn in fus_reads[fusnod]:
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

    var_calls = {}
    var_calls['var_reads_f'] = var_reads_f
    var_calls['var_reads_ref'] = var_reads_ref
    var_calls['fus_reads_f'] = fus_reads_f
    var_calls['read_sum'] = read_sum
    var_calls['cand_func_reads'] = cand_func_reads
    var_calls['var_node'] = var_node
    return (var_calls)
