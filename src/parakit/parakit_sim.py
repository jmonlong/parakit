import parakit.parakit_path as pkpath
import parakit.parakit_io as pkio
import parakit.parakit_variants as pkvar
import random


def readGfaForSim(gfa_fn, refname='grch38'):
    gfa_info = {}
    # read GFA file
    inf = open(gfa_fn, 'rt')
    gfa_info['seq'] = {}
    gfa_info['is_same_sense'] = {}
    refpath = None
    for line in inf:
        line = line.rstrip().split('\t')
        if line[0] == 'P':
            if refpath is None and line[1] == refname:
                refpath = []
                for pnode in line[2].split(','):
                    node = pnode[:-1]
                    if pnode[-1] == '+':
                        refpath.append([node, '>'])
                    else:
                        refpath.append([node, '>'])
        elif line[0] == 'W':
            if refpath is None and line[3] == refname:
                refpath = pkio.parsePath(line[6])
        elif line[0] == 'S':
            # save node sequence
            gfa_info['seq'][line[1]] = line[2]
        elif line[0] == 'L':
            if line[1] not in gfa_info['is_same_sense']:
                gfa_info['is_same_sense'][line[1]] = {}
            gfa_info['is_same_sense'][line[1]][line[3]] = line[2] == line[4]
    inf.close()
    # check sense in the reference path
    gfa_info['is_rev_in_ref'] = {}
    for node in refpath:
        gfa_info['is_rev_in_ref'][node[0]] = node[1] == '<'
    return gfa_info


def simulateHaplotype(hap_info, nodes, gfa_info, annot_fn=None,
                      pos_offset=None):
    # read clinvar variants (if available) to "protect" important
    # variants from the random conversions later
    protected_c2_nodes = set()
    if annot_fn is not None and pos_offset is not None:
        # look for read support for variant edges in vedges
        vars = pkvar.ConvertedVariants()
        vars.decomposePangenome(nodes)
        if annot_fn is not None:
            vars.matchAnnotation(annot_fn, pos_offset)
        vars.offsetPositions(pos_offset)
        for varid in vars.variants:
            var = vars.variants[varid]
            # check if matches something the annotation
            # protect it if not one of the requested
            if var.clinvar is not None:
                protected_c2_nodes.add(var.alt_trav[0])
    # simulate the path
    path = pkpath.makeFlankConsensus(nodes, flank_type='flankl')
    # prepare the buffer path
    buffer_path = pkpath.makeFlankConsensus(nodes, flank_type='buffer')
    for ii, mod_info in enumerate(hap_info):
        path += simulateModulePath(nodes, mod_info, protected_c2_nodes)
        # add buffer if necessary
        if ii < len(hap_info) - 1 and len(buffer_path) > 0:
            path += buffer_path
    path += pkpath.makeFlankConsensus(nodes, flank_type='flankr')
    # create sequence
    seq = []
    prev_sense_fwd = True
    prev_node = None
    for node in path:
        # should we write the node in forward or reverse complement ?
        cur_rev_fwd = True
        if len(seq) == 0:
            # check the sense on the reference path
            if gfa_info['is_rev_in_ref'][node]:
                cur_rev_fwd = False
        else:
            # check the sense compare to previous node
            if gfa_info['is_same_sense'][prev_node][node]:
                cur_rev_fwd = prev_sense_fwd
            else:
                cur_rev_fwd = not prev_sense_fwd
        # add the sequence
        if cur_rev_fwd:
            seq.append(gfa_info['seq'][node])
        else:
            seq.append(reverseComplement(gfa_info['seq'][node]))
        # save the current state
        prev_node = node
        prev_sense_fwd = cur_rev_fwd
    return ''.join(seq)


def simulateModulePath(nodes, mod_info, protected_c2_nodes):
    """Simulate the path of a module

    Uses the input module info to simulate a walk through the pangenome.

    Args:
        nodes : dict with node information
        mod_info : dict with simulation "instructions"
            start_mod -> start with this module (c1/c2)
            alt_paths -> dict associating a node to an alternate path
            fusions -> set with node where to insert a fusion breakpoint
            mod_noise -> proportion of nodes from the other module to include
        protected_c2_nodes: set of nodes to protect from mod_noise

    Returns: a list of nodes
    """
    mod_noise = .05
    if 'mod_noise' in mod_info:
        mod_noise = mod_info['mod_noise']
    if 'alt_paths' not in mod_info:
        mod_info['alt_paths'] = {}
    if 'fusions' not in mod_info:
        mod_info['fusions'] = {}
    path = []
    cur_module = mod_info['start_mod']
    other_module = 'c1' if cur_module == 'c2' else 'c2'
    # start from the left cycling node
    cyc_l, cyc_r = findCyclingNodes(nodes)
    path.append(cyc_l)
    # add new nodes until we reach the end of the module
    alt_paths_done = set()
    while path[-1] != cyc_r:
        node = path[-1]
        # check if a special variant should be inserted here
        if node in mod_info['alt_paths']:
            path += mod_info['alt_paths'][node].split('_')[1:]
            alt_paths_done.add(node)
            continue
        # check if a fusion should happen here
        if node in mod_info['fusions']:
            cur_module = 'c1' if cur_module == 'c2' else 'c2'
            other_module = 'c1' if cur_module == 'c2' else 'c2'
        # add one of the successor nodes, potentially with noise
        add_noise = random.random() < mod_noise
        next_node = None
        for nnode in nodes[path[-1]]['sucs']:
            if next_node is None:
                next_node = nnode
            if nodes[nnode]['class'] != 'c1' and nodes[nnode]['class'] != 'c2':
                next_node = nnode
            if (nodes[nnode]['class'] == cur_module
                    and (not add_noise or node in protected_c2_nodes)):
                next_node = nnode
                break
            if (nodes[nnode]['class'] == other_module
                    and add_noise):
                # don't create noise that affect a protected node (if in c2)
                if cur_module == 'c1' or node not in protected_c2_nodes:
                    next_node = nnode
                    break
        # if (nodes[next_node]['class'] != cur_module
        #         and not add_noise):
        #     print('{}: {} -> {}'.format(cur_module, node, next_node))
        #     for nnode in nodes[path[-1]]['sucs']:
        #         print('\t{} {}'.format(nnode, nodes[nnode]))
        path.append(next_node)
    for node in mod_info['alt_paths']:
        if node not in alt_paths_done:
            print("Warning: could'n simulate specified variant at"
                  "{}. Maybe try again.".format(node))
    return path


def reverseComplement(seq):
    """Reverse-complement a sequence"""
    rcd = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    res = []
    for base in seq:
        res.append(rcd[base])
    res.reverse()
    return ''.join(res)


def findCyclingNodes(nodes):
    # find cycle's boundaries and start node
    cyc_l = None
    cyc_r = None
    for node in nodes:
        if nodes[node]['class'] == 'cyc_l':
            cyc_l = node
        elif nodes[node]['class'] == 'cyc_r':
            cyc_r = node
    return (cyc_l, cyc_r)


def simulateReads(ref_seq, fastq_out_con,
                  read_len=[10000, 20000], depth=10,
                  substitution_rate=.00001, indel_rate=.00001,
                  substitution_erate=.01, indel_erate=.001):
    # add variants
    ref_seq_mut = []
    in_del = 0
    for base in ref_seq:
        if in_del > 0:
            in_del += -1
        elif random.random() < substitution_rate:
            other_nuc = random.sample(['A', 'T', 'C', 'G'], 1)
            while other_nuc == base:
                other_nuc = random.sample(['A', 'T', 'C', 'G'], 1)
            ref_seq_mut.append(other_nuc[0])
        elif random.random() < indel_rate:
            if random.random() < .5:
                in_del = random.randint(2, 7)
            else:
                other_nuc = random.sample(['A', 'T', 'C', 'G'],
                                          random.randint(2, 7),
                                          counts=[30] * 4)
                ref_seq_mut.append(base)
                ref_seq_mut += other_nuc
        else:
            ref_seq_mut.append(base)
    ref_seq_mut = ''.join(ref_seq_mut)
    # simulate reads
    tot_nuc = 0
    ref_len = len(ref_seq_mut)
    max_tot_nuc = ref_len * depth
    ref_hash = abs(hash(ref_seq_mut))
    read_cpt = 0
    while tot_nuc < max_tot_nuc:
        # draw a read length
        rlen = random.randint(read_len[0], read_len[1])
        # pick a staring position
        spos = random.randint(0, ref_len - rlen - 1)
        # get sequence and add errors
        rseq = []
        in_del = 0
        for rpos in range(spos, spos + rlen):
            if in_del > 0:
                in_del += -1
            elif random.random() < substitution_erate:
                other_nuc = random.sample(['A', 'T', 'C', 'G'], 1)
                while other_nuc == base:
                    other_nuc = random.sample(['A', 'T', 'C', 'G'], 1)
                rseq.append(other_nuc[0])
            elif random.random() < indel_erate:
                if random.random() < .5:
                    in_del = random.randint(2, 7)
                else:
                    other_nuc = random.sample(['A', 'T', 'C', 'G'],
                                              random.randint(2, 7),
                                              counts=[30] * 4)
                    rseq.append(ref_seq_mut[rpos])
                    rseq += other_nuc
            else:
                rseq.append(ref_seq_mut[rpos])
        # reverse complement half of the reads
        if random.random() < .5:
            rseq = reverseComplement(rseq)
        # write to FASTQ
        bq = ''.join(['~'] * len(rseq))
        rseq = ''.join(rseq)
        fastq_out_con.write('@{}_{}\n{}\n+\n{}\n'.format(ref_hash,
                                                         read_cpt,
                                                         rseq,
                                                         bq))
        read_cpt += 1
        tot_nuc += len(rseq)
