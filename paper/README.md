# Data and figures

Data to reproduce the figures is available in the [`data`](data) directory.
Scripts to reproduce some of the figures were placed in the [`figures`](figures) directory.

To test and reproduce an analysis with Parakit, we are sharing the data from one sample in `data/testsample.chr6_31800000_32200000.bam`.
The files for the pangenome and annotation are also provided so the *Pangenome construction* step can be skipped.

If you prefer to use the docker container with all the dependencies installed, start it with something like:

```sh
docker run -it -v `pwd`:/app -w /app -u `id -u $USER` quay.io/jmonlong/parakit:1.0.0
```

If you just want to see what the output files would look like, see the [`testsample_output`](testsample_output) directory.

# Methods

Below, we describe briefly the methods implemented in Parakit as used for the manuscript (*link and info soon*).
We built a pangenome for the RCCX region that can contain a module with the CYP21A1P pseudogene (module *P*), and a module with the CYP21A2 gene (module *G*).

## Pangenome construction

Parakit will work with the sequence of the region of interest and the module *P* and *G*.
Hence, the first step is to specify the coordinates of both modules in GRCh38 (chr6:31980532-32013273 for module *P*, chr6:32013273-32046127 for module *G*) and the size of the flanking region to include (300 Kbp). 
Parakit starts by creating a reference sequence without module *G*, to use as a backbone for the pangenome. 
It then creates a pangenome using Minigraph-Cactus, augmenting the backbone with sequences representing the module *G* extracted from GRCh38, and module *P* and *G* extracted from 65 high-quality assemblies from the HPRC. 
The HPRC assemblies were selected because both the CAT and Ensembl gene annotations consistently identified one module *P* and one module *G*.
Of note, Parakit also includes the sequence of module *P* and *G* from GRCh38 (here the sequence corresponding to chr6:31980532-32046127) as an input sequence for Minigraph-Cactus, to make sure there is an edge from the end of the collapsed module to its beginning. 
This cycling edge is necessary to allow a read or haplotype to traverse multiple modules. 
The full GRCh38 reference path is added to the pangenome using the `augment` subcommands of the vg toolkit

After building the pangenome, Parakit annotates the nodes based on the number of module *P* and *G* traversing them.
A node is marked as specific to module *G* if it is traversed by at least 3 times more module *G* than module *P*, and vice versa.

The command used to construct the RCCX pangenome with Parakit was:

```
parakit construct -j rccx.grch38_hprc.mc.config.json
```

The JSON configuration file contains information about the location of the GRCh38 reference FASTA file, coordinates of each module, size of the flanking regions, location of the HPRC assemblies.

Note: to reproduce this step, see detailed instructions in the [`../data`](../data) directory on how to download and prepare all the necessary files. 
Running the command from this directory **will not** work.
Again, this step is not required to for commands below because we are sharing the pangenome files.

## Alignment of long-reads to the pangenome

Parakit first extracts reads in the region of interest from an indexed BAM. 
This assumes that the reads were first aligned to the reference genome, for example using minimap2.
Although the read alignment around the RCCX module might not be accurate, reads of interest should still be mapped to the region. 
Hence, we should retrieve all informative reads by extracting those that were originally mapped to this region of the linear reference.

Once retrieved, the reads are re-mapped to the pangenome using GraphAligner.
Internally, Parakit uses GraphAligner v1.0.17 with variation graph mode (`-x vg`) and 100 bp for the alignment bandwidth (`-b 100`).

To extract and map long reads to the pangenome with Parakit, we ran:

```
parakit map -j rccx.grch38_hprc.mc.config.json -b data/testsample.chr6_31800000_32200000.bam -o testsample.rccx.grch38_hprc.mc.gaf.gz
```

The output is a GAF file representing the alignment of each read through the pangenome.
Because the paralogous regions are collapsed, each read maps confidently to only one position in the pangenome.
A long read can often traverse the collapsed region of the pangenome multiple times when it spans multiple modules.

## Read-based variant calling

Aligned reads that traverse the pangenome through module-specific nodes (as defined above) can be used for inference. 
In Parakit, the *call* subcommand searches for evidence of gene fusion or gene conversion in each read separately. 
Specifically, a sliding window approach looks for positions - within the RCCX module - where a read switches from aligning to *P*-specific nodes to aligning to *G*-specific nodes. 
We used sliding windows of 20 informative markers on each side of a potential variant site.
For fusions, we select sites where at least 80% of markers in the upstream window are specific to module *P*, and 80% of markers in the downstream window are specific to the module *G*.
For small gene conversion event, we look for module *P* nodes, specifically those corresponding to known ClinVar variants, surrounded by module *G* nodes.
Here, Parakit selects candidates if there are more than three times more module *G* nodes than module *P* in the windows upstream and downstream.
For both fusions and gene conversion variants, Parakit reports sites with at least 3 supporting reads.

To call variants with Parakit, we ran:
```
parakit call -r testsample.rccx.grch38_hprc.mc.gaf.gz -j rccx.grch38_hprc.mc.config.json -o testsample.rccx.grch38_hprc.mc.calls.tsv
```

## Coverage and allele ratio along the pangenome

Parakit uses the read alignment through the pangenome to compute estimates of coverage and a ratio of *G* alleles along the RCCX region.

Read coverage helps estimate the total number of module copies. 
We count the number of reads covering each of the non-module-specific nodes (light blue points in the figures).
These nodes are not specific to a module so they should be traversed by all reads, not matter the allele, and provide an estimate of total copy number.
This coverage can be plotted along the pangenome (see Visualization of the results below).
Parakit also computes two global copy number estimates.
The first is simply the median read coverage across the non-specific nodes, as mentioned above, normalized by the coverage in flanking regions.
We note that this approach can underestimate the copy number when there is mapping bias in some parts of the region (typically lower coverage in the middle of the RCCX module).
The second copy number estimate focuses on the reads entering (or leaving) the collapsed region, and the reads cycling back from the end of the module back to the beginning.
Indeed, each cycle supports the presence of an additional copy. 
Parakit estimates the copy number as `2*(Rf+Rc)/Rf` where `Rf` is the number of reads entering (or exiting) from the flanks and `Rc` is the number of cycles supported by reads.

Parakit also helps look for changes in the allelic balance at informative sites, i.e. the proportion of reads traversing nodes specific to module *G*.
Changes in this proportion is expected around breakpoints of a fusion allele or at the boundaries of a large gene-converted region. 
Each module-specific node is assigned an *anchor* node as the first non-specific reference node  upstream.
For each anchor node, we compute the coverage of *G*-specific nodes divided by the coverage of nodes specific to either *P* or *G*.
In the presence of two bimodular haplotypes, the ratio should be centered around 0.5.
An individual carrying a fusion haplotype and a bimodular haplotype (one module *P* and one module *G*) will show consistent allele *G* ratio around 1/3 up to the fusion breakpoint, and around 2/3 afterward.
This evidence is orthogonal to the read-level evidence extracted by the variant calling or haplotype reconstruction approaches.

Both of those analyses are performed by the visualization command (see below) using the read alignment information.

## Haplotype reconstruction

Parakit also infers the most likely pair of haplotypes (or diplotype) by finding the pair of paths through the pangenome graph that are the most consistent with the aligned reads. 
Haplotype reconstruction is performed in two steps: 1) haplotype candidates are enumerated, 2) the best pair of haplotype candidates is identified.

Haplotype candidates are first enumerated. 
The goal is to produce an extensive list of potential haplotypes guided by the read alignments.
This is also done in two steps: module candidates are constructed and then stitched together into haplotypes.
Read alignments are first split at the boundaries of the collapsed region of the pangenome. 
For example, a read spanning the flanking region and two modules will be split in three subreads: the upstream flank, the first module, and the second module traversed.
All the subreads are then clustered using an iterative approach. 
All subreads start in the same cluster and a consensus path is produced. 
Parakit then looks for positions in this consensus where at least three reads disagree.
If such *variant* markers are found, they are used to split the cluster in two.
To do so, a network is built where reads are connected by edges that are weighted by the number of variant markers they have in common.
The reads are split in two sets using the Kernighan–Lin algorithm, implemented in the networkx Python package.
This process is iterated until no *variant* markers are found in the clusters.
The consensus path for all clusters (and all iterations) are saved and combined in the final step into candidate haplotypes.
Here, Parakit uses information about which subreads are assigned to a cluster and how the subreads were originally split to suggest adjacency between the clusters. 
Haplotypes are created by starting at the upstream flank and adding clusters if at least one read supports this adjacency, until it reaches the downstream flank or the maximum allowed number of module is exceeded (currently set to 5).

The enumeration of the haplotype candidates described above is permissive on purpose to maximize the chances of containing the correct ones. 
In the second step of the haplotype reconstruction, Parakit looks for the pair of haplotypes that best match the read alignments, both in terms of read identity and read coverage.
Read identity is measured as the average identity of the best alignment on the pair of haplotype candidates. 
This identity is computed in pangenome space, i.e. representing the proportion of nodes matching when aligning the reads and haplotype.
Hence, this pangenome alignment identity gives more weight to variants than sequence alignment because (potentially long) stretches of common sequences are merged into single nodes.
Read coverage uniformity is measured as the average deviation of the read coverage in each node compared to the expected coverage.
The expected read coverage is computed based on the number of times the node is present in the tested pair of haplotypes.
Both the identity and read coverage metrics are normalized by their maximum values across all evaluated pairs, summed and used to rank the pairs of haplotypes. 
The haplotype pair with the highest normalized score is selected as the most likely diplotype.

In summary, reads are represented by their traversal of the pangenome, naturally putting more weight on known variants, and used to predict potential haplotypes. 
The pair of haplotypes that results in the highest read identity and most uniform read coverage is then reported by Parakit as the most likely diplotype.
The command to reconstruct haplotypes was:

```
parakit diplotype -j rccx.grch38_hprc.mc.config.json -r testsample.rccx.grch38_hprc.mc.gaf.gz -o testsample.rccx.grch38_hprc.mc
```

## Visualization of the results

Parakit uses R and the ggplot2 package to display the different layers of evidence described above on the collapsed pangenome. 
The different panels are horizontally aligned so that sites of variation (e.g. fusion, SNVs) can be compared easily between analyses.
For example, to help check that the fusion predicted by reads is concordant with the switch in module *G* ratio and the predicted fusion in the reconstructed diplotype.
It includes a gene annotation track with the exons and introns for the genes in the region. 
The command used to generate a PDF file with the graph was:

```
parakit viz -j rccx.grch38_hprc.mc.config.json -r testsample.rccx.grch38_hprc.mc.gaf.gz -c testsample.rccx.grch38_hprc.mc.calls.tsv -d testsample.rccx.grch38_hprc.mc.paths-stats.tsv -p testsample.rccx.grch38_hprc.mc.paths-info.tsv -o testsample.rccx.grch38_hprc.mc.pdf
```

The PDF shows evidence of a trimodular allele carrying a pathogenic SNV, and a fusion allele (see [pdf](testsample_output/testsample.rccx.grch38_hprc.mc.pdf)).
