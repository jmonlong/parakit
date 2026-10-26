import matplotlib.pyplot as plt
import parakit.parakit_path as pkpath


def plotPackedSubreads(nodes, reads, genes, out_fn):
    fig, ax = plt.subplots(2, 1, sharex=True,
                           height_ratios=[10, 1])

    preparePackedSubreadPlot(ax[0], reads, nodes, size=3)

    prepareGenePlot(ax[1], genes)

    plt.subplots_adjust(hspace=.05, top=.99, bottom=.07, right=.99)

    fig.set_size_inches(10, 5)
    plt.savefig(out_fn)


def plotDiplotypeReads(nodes, reads, dip_paths, out_fn):
    fig, ax = plt.subplots(4, 1, sharex=True,
                           height_ratios=[4, 1, 4, 1])

    # haplotype names
    hns = list(dip_paths.keys())

    # map reads on diplotype
    reads_on_dip = mapReadsToDiplotype(reads, dip_paths)

    # prepare the read panel for each haplotype
    prepareReadOnHapPlot(ax[0], reads, reads_on_dip, dip_paths, hns[0], nodes)
    prepareReadOnHapPlot(ax[2], reads, reads_on_dip, dip_paths, hns[1], nodes)

    # prepare the colored node panel for each haplotype
    prepareHapPlot(ax[1], nodes, dip_paths[hns[0]])
    prepareHapPlot(ax[3], nodes, dip_paths[hns[1]])

    ax[0].axis('off')
    ax[2].axis('off')

    ax[0].title.set_text(hns[0])
    ax[2].title.set_text(hns[1])

    fig.set_size_inches(9, 6)
    plt.tight_layout()
    plt.savefig(out_fn)


def plotLongReads(nodes, reads, genes, out_fn):
    fig, ax = plt.subplots(2, 1, sharex=True,
                           height_ratios=[10, 1])

    prepareLongReadPlot(ax[0], reads, nodes, size=6, n_top_reads=30)

    prepareGenePlot(ax[1], genes)

    plt.subplots_adjust(hspace=.05, top=.99, bottom=.07, right=.99)

    fig.set_size_inches(10, 5)
    plt.savefig(out_fn)


def prepareReadOnHapsPlot(ax, reads, dip_paths, nodes):
    # haplotype names
    hns = list(dip_paths.keys())
    # how much to shift the position of the second haplotype
    offset = len(dip_paths[hns[0]]) + 10
    # prepare flank positions
    flank_pos = set()
    hap_offset = 0
    for hapn in hns:
        for idx, node in enumerate(dip_paths[hapn]):
            if nodes[node]['class'] == 'ref':
                flank_pos.add(idx + hap_offset)
        hap_offset = offset
    # map reads to haplotypes
    hap_alns = {}
    for hapn in hns:
        hap_alns[hapn] = pkpath.mapReadsToHap(dip_paths[hapn], reads.path)
    # keep the best alignment
    hap_b_alns = {}
    for rn in reads.path:
        if hap_alns[hns[0]][rn]['score'] > hap_alns[hns[1]][rn]['score']:
            hap_b_alns[rn] = hap_alns[hns[0]][rn]
        elif hap_alns[hns[0]][rn]['score'] < hap_alns[hns[1]][rn]['score']:
            hap_b_alns[rn] = hap_alns[hns[1]][rn]
            # shift position on second haplotype
            aln_len = len(hap_b_alns[rn]['aln_pos'])
            for apos in range(aln_len):
                hap_b_alns[rn]['aln_pos'][apos] += offset
    # prepare the read ranges
    read_r = {}
    for rn in hap_b_alns:
        read_r[rn] = (hap_b_alns[rn]['aln_pos'][0],
                      hap_b_alns[rn]['aln_pos'][-1])
    # sort them by decreasing length
    rn_sorted = sorted(read_r, key=lambda k: read_r[k][1] - read_r[k][0],
                       reverse=True)
    read_y = rangesToYpos(read_r, rn_sorted)

    # add lines
    for rn in rn_sorted:
        x = [read_r[rn][0], read_r[rn][1]]
        y = [read_y[rn], read_y[rn]]
        ax.plot(x, y, 'b-', linewidth=.1)

    # add tiny points
    pt_x = []
    pt_y = []
    pt_col = []
    pt_size = []
    pt_alpha = []
    fl_x = []
    fl_y = []
    for rn in rn_sorted:
        rlen = len(hap_b_alns[rn]['aln_pos'])
        for pp in range(rlen):
            pt_x.append(hap_b_alns[rn]['aln_pos'][pp])
            pt_y.append(read_y[rn])
            if hap_b_alns[rn]['aln_type'][pp] == 'mismatch':
                pt_col.append('black')
                pt_size.append(5)
                pt_alpha.append(.8)
            else:
                pt_col.append('b')
                pt_size.append(1)
                pt_alpha.append(.3)
            if hap_b_alns[rn]['aln_pos'][pp] in flank_pos:
                fl_x.append(hap_b_alns[rn]['aln_pos'][pp])
                fl_y.append(read_y[rn])
    # plot flank squares first
    if len(fl_x) > 0:
        ax.scatter(fl_x, fl_y, color='orange', s=8, alpha=.8,
                   marker='s')
    ax.scatter(pt_x, pt_y, color=pt_col, s=pt_size, alpha=pt_alpha,
               marker='o')


def mapReadsToDiplotype(reads, dip_paths):
    # haplotype names
    hns = list(dip_paths.keys())
    # map reads to haplotypes
    hap_alns = {}
    for hapn in hns:
        hap_alns[hapn] = pkpath.mapReadsToHap(dip_paths[hapn], reads.path)
    # keep the best alignment
    hap_b_alns = {}
    hap_b_alns[hns[0]] = {}
    hap_b_alns[hns[1]] = {}
    for rn in reads.path:
        if hap_alns[hns[0]][rn]['score'] > hap_alns[hns[1]][rn]['score']:
            hap_b_alns[hns[0]][rn] = hap_alns[hns[0]][rn]
        elif hap_alns[hns[0]][rn]['score'] < hap_alns[hns[1]][rn]['score']:
            hap_b_alns[hns[1]][rn] = hap_alns[hns[1]][rn]
    return hap_b_alns


def prepareReadOnHapPlot(ax, reads, reads_on_dip, dip_paths, hap_name, nodes):
    hap_b_alns = reads_on_dip[hap_name]
    # prepare flank positions
    flank_pos = set()
    for idx, node in enumerate(dip_paths[hap_name]):
        if nodes[node]['class'] == 'ref':
            flank_pos.add(idx)
    # prepare the read ranges
    read_r = {}
    for rn in hap_b_alns:
        read_r[rn] = (hap_b_alns[rn]['aln_pos'][0],
                      hap_b_alns[rn]['aln_pos'][-1])
    # sort them by decreasing length
    rn_sorted = sorted(read_r, key=lambda k: read_r[k][1] - read_r[k][0],
                       reverse=True)
    read_y = rangesToYpos(read_r, rn_sorted)

    # add lines
    for rn in rn_sorted:
        x = [read_r[rn][0], read_r[rn][1]]
        y = [read_y[rn], read_y[rn]]
        ax.plot(x, y, 'b-', linewidth=.1)

    # add tiny points
    pt_x = []
    pt_y = []
    pt_col = []
    pt_size = []
    pt_alpha = []
    fl_x = []
    fl_y = []
    for rn in rn_sorted:
        rlen = len(hap_b_alns[rn]['aln_pos'])
        for pp in range(rlen):
            pt_x.append(hap_b_alns[rn]['aln_pos'][pp])
            pt_y.append(read_y[rn])
            if hap_b_alns[rn]['aln_type'][pp] == 'mismatch':
                pt_col.append('black')
                pt_size.append(5)
                pt_alpha.append(.8)
            else:
                pt_col.append('b')
                pt_size.append(1)
                pt_alpha.append(.3)
            if hap_b_alns[rn]['aln_pos'][pp] in flank_pos:
                fl_x.append(hap_b_alns[rn]['aln_pos'][pp])
                fl_y.append(read_y[rn])
    # plot flank squares first
    if len(fl_x) > 0:
        ax.scatter(fl_x, fl_y, color='orange', s=8, alpha=.8,
                   marker='s')
    ax.scatter(pt_x, pt_y, color=pt_col, s=pt_size, alpha=pt_alpha,
               marker='o')


def prepareHapsPlot(ax, nodes, dip_paths):
    xpos = 0
    pt_x = []
    pt_y = []
    pt_colors = []
    pt_size = []
    pt_alpha = []
    sep_pos = None
    mod_sep = []
    for hapn in dip_paths:
        for node in dip_paths[hapn]:
            pt_x.append(xpos)
            if nodes[node]['class'] == 'c1':
                pt_colors.append('r')
                pt_alpha.append(.8)
                pt_size.append(4)
                pt_y.append(-1)
            elif nodes[node]['class'] == 'c2':
                pt_colors.append('g')
                pt_alpha.append(.8)
                pt_size.append(4)
                pt_y.append(1)
            else:
                pt_colors.append('b')
                pt_alpha.append(.01)
                pt_size.append(.2)
                pt_y.append(0)
            if nodes[node]['class'] in ['cyc_l', 'cyc_r']:
                mod_sep.append(xpos)
            xpos += 1
        if sep_pos is None:
            sep_pos = xpos + 5
        xpos += 10

    ax.scatter(pt_x, pt_y, s=pt_size, color=pt_colors, alpha=pt_alpha)
    ax.set_yticks([-1, 0, 1], ['c1', 'ref', 'c2'])
    for pos in mod_sep:
        ax.axvline(x=pos, linestyle='--', alpha=.3, linewidth=1, color='black')
    ax.axvline(x=sep_pos, linestyle='-', linewidth=2, color='black')
    ax.margins(0, .3)


def prepareHapPlot(ax, nodes, hap_paths):
    xpos = 0
    pt_x = []
    pt_y = []
    pt_colors = []
    pt_size = []
    pt_alpha = []
    mod_sep = []
    for node in hap_paths:
        pt_x.append(xpos)
        if nodes[node]['class'] == 'c1':
            pt_colors.append('r')
            pt_alpha.append(.8)
            pt_size.append(4)
            pt_y.append(-1)
        elif nodes[node]['class'] == 'c2':
            pt_colors.append('g')
            pt_alpha.append(.8)
            pt_size.append(4)
            pt_y.append(1)
        else:
            pt_colors.append('b')
            pt_alpha.append(.01)
            pt_size.append(.2)
            pt_y.append(0)
        if nodes[node]['class'] in ['cyc_l', 'cyc_r']:
            mod_sep.append(xpos)
        xpos += 1

    ax.scatter(pt_x, pt_y, s=pt_size, color=pt_colors, alpha=pt_alpha)
    ax.set_yticks([-1, 0, 1], ['c1', 'ref', 'c2'])
    for pos in mod_sep:
        ax.axvline(x=pos, linestyle='--', alpha=.3, linewidth=1, color='black')
    ax.margins(0, .3)


def prepareNodePlot(ax, nodes):
    pos = []
    cols = []
    yvals = []
    alphas = []
    for node in nodes:
        pos.append(int(node))
        if nodes[node]['class'] == 'c1':
            cols.append('r')
            yvals.append(-1)
            alphas.append(.8)
        elif nodes[node]['class'] == 'c2':
            cols.append('g')
            yvals.append(1)
            alphas.append(.8)
        else:
            cols.append('b')
            yvals.append(0)
            alphas.append(.01)

    ax.scatter(pos, yvals, color=cols, alpha=alphas, s=2)
    ax.set_yticks([-1, 0, 1], ['c1', 'ref', 'c2'])
    ax.margins(0, .3)


def prepareGenePlot(ax, genes):
    # prepare a y position for each gene
    gene_y = {}
    for elt in genes:
        if elt['gene_name'] not in gene_y:
            gene_y[elt['gene_name']] = len(gene_y)
    # plot the gene body and exons
    tss_x = []
    tss_y = []
    tss_c = []
    for elt in genes:
        if elt['type'] == 'gene':
            # thin line for the gene body
            x = [int(elt['node_start']), int(elt['node_end'])]
            y = [gene_y[elt['gene_name']]] * 2
            c = 'r-' if elt['mod_start'] == 'c1' else 'g-'
            ax.plot(x, y, c, linewidth=.5)
            if elt['strand'] == '+':
                tss_x.append(int(elt['node_start']))
            else:
                tss_x.append(int(elt['node_end']))
            tss_y.append(gene_y[elt['gene_name']] + .3)
            c = 'red' if elt['mod_start'] == 'c1' else 'green'
            tss_c.append(c)
        elif elt['type'] == 'exon':
            # thicker line for the exon
            x = [int(elt['node_start']), int(elt['node_start']),
                 int(elt['node_end']), int(elt['node_end'])]
            y = [gene_y[elt['gene_name']]] * 4
            y[0] += -.3
            y[1] += .3
            y[2] += .3
            y[3] += -.3
            c = 'red' if elt['mod_start'] == 'c1' else 'green'
            ax.fill(x, y, color=c)
    # add TSS marker
    ax.scatter(tss_x, tss_y, color=tss_c, s=10)
    # rename axis with gene names
    gene_names = list(gene_y.keys())
    gene_ypos = [gene_y[gn] for gn in gene_names]
    ax.set_yticks(gene_ypos, gene_names)
    ax.margins(0, .3)


def prepareLongReadPlot(ax, reads, nodes, size=3, n_top_reads=10,
                        distance_ends=False):
    # score reads with the number of cycling nodes and number of
    # informative nodes
    read_scores = {}
    for readn in reads.path:
        score = 0
        cyc = False
        for node in reads.path[readn]:
            nclass = nodes[node]['class']
            if nclass == 'cyc_r' or nclass == 'cyc_l':
                cyc = True
            elif nclass == 'c1' or nclass == 'c2':
                score += .0000001
        if cyc:
            score += 1
        read_scores[readn] = score
    # select the top N reads
    sorted_reads = sorted(read_scores, reverse=True,
                          key=lambda k: read_scores[k])
    if len(sorted_reads) > n_top_reads:
        sorted_reads = sorted_reads[:n_top_reads]
    # try to avoid reads starting right after the previous one ended
    if distance_ends:
        read_list = [sorted_reads[0]]
        while len(read_list) < len(sorted_reads):
            readn_p = read_list[-1]
            # look for the best next read
            best_read = None
            best_dist = None
            for readn in sorted_reads:
                if readn in read_list:
                    continue
                dist = abs(int(reads.path[readn][0]) -
                           int(reads.path[readn_p][-1]))
                if best_read is None or best_dist < dist:
                    best_dist = dist
                    best_read = readn
            read_list.append(best_read)
    else:
        read_list = sorted_reads
    # prepare the points
    pt_x = []
    pt_y = []
    pt_col = []
    pt_size = []
    pt_alpha = []
    cur_y = 0
    yticks = []
    for readn in read_list:
        yticks.append(cur_y)
        for nodei, node in enumerate(reads.path[readn]):
            pt_x.append(int(node))
            pt_y.append(cur_y)
            cur_y += 1
            if nodes[node]['class'] == 'c1':
                pt_col.append('r')
                pt_size.append(size)
                pt_alpha.append(.8)
            elif nodes[node]['class'] == 'c2':
                pt_col.append('g')
                pt_size.append(size)
                pt_alpha.append(.8)
            elif nodes[node]['class'] == 'ref':
                pt_col.append('black')
                pt_size.append(2*size)
                pt_alpha.append(.8)
            else:
                pt_col.append('b')
                pt_size.append(.05)
                pt_alpha.append(.3)
            # add lines
            if nodei < len(reads.path[readn]) - 1:
                x = [int(node), int(reads.path[readn][nodei + 1])]
                y = [cur_y - 1, cur_y]
                ax.plot(x, y, 'b-', linewidth=.1)
    ax.scatter(pt_x, pt_y, color=pt_col, s=pt_size, alpha=pt_alpha)
    ax.set_yticks(yticks, [''] * len(yticks))


def preparePackedSubreadPlot(ax, reads, nodes, size=3):
    # prepare subreads
    sreads = pkpath.Subreads()
    sreads.splitReads(reads, nodes)

    # prepare the subread ranges
    sr_r = {}
    for sr in sreads.sreads:
        sr_r[sr] = (int(sreads.sreads[sr].path[0]),
                    int(sreads.sreads[sr].path[-1]))
    # sort them by decreasing length
    srn_sorted = sorted(sr_r, key=lambda k: sr_r[k][1] - sr_r[k][0],
                        reverse=True)

    # pack position each subreads
    sr_y = rangesToYpos(sr_r, srn_sorted)

    # add lines
    for sr in srn_sorted:
        x = [sr_r[sr][0], sr_r[sr][1]]
        y = [sr_y[sr], sr_y[sr]]
        ax.plot(x, y, 'b-', linewidth=.1)

    # add tiny points
    pt_x = []
    pt_y = []
    pt_col = []
    pt_size = []
    pt_alpha = []
    for sr in srn_sorted:
        for node in sreads.sreads[sr].path:
            pt_x.append(int(node))
            pt_y.append(sr_y[sr])
            if nodes[node]['class'] == 'c1':
                pt_col.append('r')
                pt_size.append(size)
                pt_alpha.append(.8)
            elif nodes[node]['class'] == 'c2':
                pt_col.append('g')
                pt_size.append(size)
                pt_alpha.append(.8)
            else:
                pt_col.append('b')
                pt_size.append(.05)
                pt_alpha.append(.3)
    ax.scatter(pt_x, pt_y, color=pt_col, s=pt_size, alpha=pt_alpha)


def rangesToYpos(ranges_d, range_names):
    # pack position each subreads
    used_y = {}
    rg_y = {}
    for rg in range_names:
        x1 = ranges_d[rg][0]
        x2 = ranges_d[rg][1]
        ypos = 0
        search_ypos = True
        while search_ypos:
            if ypos in used_y:
                gap = True
                for xx in range(x1, x2 + 1):
                    if xx in used_y[ypos]:
                        gap = False
                        break
                if gap:
                    # we've found a gap, stop searching
                    search_ypos = False
                else:
                    # try the next y position
                    ypos += 1
            else:
                # we've found a gap, stop searching
                used_y[ypos] = set()
                search_ypos = False
        rg_y[rg] = ypos
        # mark as used space
        for xx in range(x1, x2 + 1):
            used_y[ypos].add(xx)
    return rg_y
