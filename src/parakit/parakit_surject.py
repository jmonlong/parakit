import parakit.parakit_path as pkpath


def selectSubreads(reads, nodes, module='c2'):
    """Select subreads from a module

    Split a set of reads and return the subreads that seem to
    originate from a specific module. The module assignment is based
    on the informative markers.

    Args:
        reads : a Reads object
        module : which module to select

    Returns: a dict with the subreads informations
        subread name -> - read name
                        - start position in read
                        - end position in read
    """
    # split reads into subreads
    sreads = splitReads(reads, nodes)
    # select those that are mostly in the module and the interval
    sreads_sel = {}
    for sreadn in sreads:
        # count how many informative nodes are from the target module
        inf_nodes = 0
        mod_nodes = 0
        for node in sreads[sreadn]['path']:
            if nodes[node]['class'] == 'c1' or nodes[node]['class'] == 'c2':
                inf_nodes += 1
            if nodes[node]['class'] == module:
                mod_nodes += 1
        if inf_nodes > 5 and float(mod_nodes) / inf_nodes > .5:
            # consider this subreads as potentially part of the target module
            sreads_sel[sreadn] = {}
            sreads_sel[sreadn]['read'] = sreads[sreadn]['read']
            sreads_sel[sreadn]['start'] = min(sreads[sreadn]['pos'])
            sreads_sel[sreadn]['end'] = max(sreads[sreadn]['pos'])
    return sreads_sel


def selectSubreadsOnDip(reads, nodes, dip_paths):
    """Select subreads on modules relative to a diplotype

    Align the reads to a diplotype, then split them into subreads,
    annotated by which haplotype/module they mapped to.

    Args:
        reads : a Reads object
        nodes : dict with node information
        dip_paths : dict haplotype name -> path (node list)

    Returns: a dict with the subreads informations
        subread name -> - 'read' name
                        - 'start' position in read
                        - 'end' position in read
                        - 'haplotype' name
                        - 'module' name
    """
    # map reads to haplotypes
    hns = list(dip_paths.keys())
    hap_alns = {}
    for hapn in hns:
        hap_alns[hapn] = pkpath.mapReadsToHap(dip_paths[hapn], reads.path)
    # keep the best alignment
    hap_b_alns = {}
    for rn in reads.path:
        if hap_alns[hns[0]][rn]['score'] > hap_alns[hns[1]][rn]['score']:
            hap_b_alns[rn] = hap_alns[hns[0]][rn]
            hap_b_alns[rn]['hap'] = hns[0]
        elif hap_alns[hns[0]][rn]['score'] < hap_alns[hns[1]][rn]['score']:
            hap_b_alns[rn] = hap_alns[hns[1]][rn]
            hap_b_alns[rn]['hap'] = hns[1]
        else:
            # same score, assign semi-randomly to one haplotype
            if hash(rn) % 2 == 0:
                hap_b_alns[rn] = hap_alns[hns[0]][rn]
                hap_b_alns[rn]['hap'] = hns[0]
            else:
                hap_b_alns[rn] = hap_alns[hns[1]][rn]
                hap_b_alns[rn]['hap'] = hns[1]
    # split haplotypes
    haps_annot = {}
    haps_annot[hns[0]] = annotateHaplotype(dip_paths[hapn], nodes)
    haps_annot[hns[1]] = annotateHaplotype(dip_paths[hapn], nodes)
    # split reads into subreads
    sreads = splitReads(reads, nodes)
    # annotate with the  those that are mostly in the module and the interval
    sreads_sel = {}
    for sreadn in sreads:
        readn = sreads[sreadn]['read']
        # find assigned haplotype
        hapn = hap_b_alns[readn]['hap']
        # find most common module
        mod_d = {}
        for node_idx in sreads[sreadn]['node_idx']:
            if node_idx in hap_b_alns[readn]['read_hap_pos']:
                hap_idx = hap_b_alns[readn]['read_hap_pos'][node_idx]
                mod = haps_annot[hapn][hap_idx]
                if mod not in mod_d:
                    mod_d[mod] = 0
                mod_d[mod] += 1
        mod_s = sorted(list(mod_d.keys()), key=lambda k: mod_d[k],
                       reverse=True)
        mod_s = mod_s[0]
        # prepare the subread info
        sreads_sel[sreadn] = {}
        sreads_sel[sreadn]['read'] = sreads[sreadn]['read']
        sreads_sel[sreadn]['start'] = min(sreads[sreadn]['pos'])
        sreads_sel[sreadn]['end'] = max(sreads[sreadn]['pos'])
        sreads_sel[sreadn]['haplotype'] = hapn
        sreads_sel[sreadn]['module'] = '{}_{}'.format(hapn, mod_s)
    return sreads_sel


def splitReads(reads, nodes):
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
    res = {}
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
        # save the read positions
        spos = [[]]
        # save the node index range
        node_idx = [[]]
        add_subread = False
        for node_ii in range(len(reads.path[readn])):
            nod = reads.path[readn][node_ii]
            if nod == cyc_l:
                if len(subreads[-1]) > 0:
                    add_subread = True
            if nod == cyc_r:
                subreads[-1].append(nod)
                spos[-1].append(reads.readpos[readn][node_ii])
                node_idx[-1].append(node_ii)
                add_subread = True
            else:
                if add_subread:
                    subreads.append([])
                    spos.append([])
                    node_idx.append([])
                    add_subread = False
                subreads[-1].append(nod)
                spos[-1].append(reads.readpos[readn][node_ii])
                node_idx[-1].append(node_ii)
        subreads_t = []
        for spath in subreads:
            group = None
            for nod in spath:
                if nod in node_group:
                    group = node_group[nod]
                    break
            if group is None:
                print('Warning: problem splitting read')
                print(subreads)
            subreads_t.append(group)
        # skip if just one subread and not in the region of interest
        if len(subreads_t) == 1 and subreads_t[0] != 'module':
            continue
        # otherwise, save the subreads
        for sbi, sr_path in enumerate(subreads):
            sreadn = '{}_{}'.format(readn, sbi)
            subread = {}
            subread['read'] = readn
            subread['path'] = sr_path
            subread['pos'] = spos[sbi]
            subread['node_idx'] = node_idx[sbi]
            res[sreadn] = subread
    return res


def annotateHaplotype(path, nodes):
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
    # group the nodes into flankl, flankr, buffer and module
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
    # process the path
    annot = []
    cur_mod = 0
    cur_group = 'flankl'
    for nod in path:
        if nod in node_group:
            cur_group = node_group[nod]
        if nod == cyc_r:
            cur_mod += 1
        if cur_group == 'module':
            annot.append(str(cur_mod))
        else:
            annot.append(cur_group)
    return annot
