# Parakit

<img src="parakit.logo.svg" width="200">

Parakit is a tool to analyze the [RCCX module](https://en.wikipedia.org/wiki/RCCX), which contain the CYP21A2 gene, using long sequencing reads. 
A first version was used to analyze ONT R10 data and made use of 65 high-quality assemblies from the [Human Pangenome Reference Consortium](https://humanpangenome.org/). 
This version is described and benchmarked in [Monlong et al. medRxiv 2025](https://www.medrxiv.org/content/10.1101/2025.02.07.25321404v1), and more information can be found in the [`paper`](paper) directory.
It is still in active development, especially to extend this approach to other loci and keep improving the resolution (new features, larger pangenome, other sequencing technologies).

- [Installation](#installation)
- [GRCh38+HPRC RCCX pangenome](#GRCh38-HPRC-RCCX-pangenome)
- [Commands](#commands)
- [Output](#output)
- [Citation](#citation)
- [FAQ](#FAQ)

Starting from an indexed BAM file, reads in the RCCX region are extracted and realigned to a local pangenome were both modules are collapsed.
Parakit then looks for pathogenic variants supported by multiple types of signal: read coverage, allele support, read support, diplotype reconstruction. 

![](docs/imgs/overview.png)

## Installation

Clone the repo and install locally with

```sh
git clone https://github.com/jmonlong/parakit.git
cd parakit
python3 -m pip install -e .
```

It might be a good idea to use a virtual environment

```sh
## to create the environment
python3 -m venv parakit_venv
## to activate it
source parakit_venv/bin/activate
## intall
## run from the directory where the repo was cloned
pip3 install -e .
## deactivate env when you're done
deactivate
```

Then, you can use Parakit anytime you activate this environment.

### Dependencies

Parakit will use the following external tools

- [vg](https://github.com/vgteam/vg)
- [samtools](https://samtools.github.io/)
- [docker](https://docs.docker.com/engine/install/)
    - Used to run [GraphAligner](https://github.com/maickrau/GraphAligner), [cactus-pangenome](https://github.com/ComparativeGenomicsToolkit/cactus), [pggb](https://github.com/pangenome/pggb) (if they are available, docker won't be used)
- R with the following packages
    - dplyr
    - ggplot2
    - tidyr
    - RColorBrewer
    - GenomicRanges
    - cowplot
    - rjson
    - rmarkdown
- gunzip

A Docker image is also available with all the dependencies: `quay.io/jmonlong/parakit:1.0.0`. 
The [`Dockerfile`](Dockerfile) can also give hints how to install all the dependencies.

## GRCh38+HPRC RCCX pangenome

Ready-to-use files for the RCCX pangenome, available in the [`data` folder](data)

- `rccx.grch38_hprc.mc.config.json` the configuration file for this pangenome (contains coordinates, flank size, etc used to build the pangenome)
- `rccx.grch38_hprc.mc.pg.gfa` the pangenome in GFA format
- `rccx.grch38_hprc.mc.node_info.tsv` metadata about the nodes in the pangenome, e.g. which one is specific to module 1/2.
- Annotation files
    - `CYP21A2.pathogenic.variant_summary.2024_09_03.tsv` reformatted subset of ClinVar including CYP21A2 pathogenic variants
    - `CYP21A2.gencodev43.nearby_genes.tsv` reformatted subset of GENCODE containing gene annotation in the region.

See [data/rccx.summary.md](data/rccx.summary.md) for some descriptive metrics on this pangenome.

The files mentioned above are the output of the pangenome construction.
They can be used to analyze a new long-read sequencing samples using the [commands](#commands) below.
For info, the steps to construct the pangenome are described in the [`data` folder](data).

## Commands

For example, to analyze one sample with an indexed BAM file (aligned to GRCh38).

To extract relevant reads and map them to the pangenome

```bash
parakit map -j rccx.grch38_hprc.mc.config.json -b input.bam -o reads.gaf.gz
```

This creates the `reads.gaf.gz` GAF file.

Then, to look for variant-supporting reads

```bash
parakit call -j rccx.grch38_hprc.mc.config.json -r reads.gaf.gz -o calls.tsv
```

The reads/calls are saved in `calls.tsv`.

To list and evaluate candidate diplotype:

```bash
parakit diplotype -j rccx.grch38_hprc.mc.config.json -r reads.gaf.gz -o diplotype
```

This command creates two files

- `diplotype.paths-stats.tsv` with the diplotypes ranked by score (based on read alignment and coverage).
- `diplotype.paths-info.tsv` with the path taken by each haplotype through the pangenome.

Of note, the diplotype inference is the module most sensitive to the input reads. 
It might not always work, especially if reads are shorter.
If you notice inconsistencies, believe the calls from the raw reads and aggregated coverage/allele support.

Finally, the visualization command makes a figure. 
The *all* mode, will make a multi-panel figure summarizing all analysis.

```bash
parakit viz -v all -j rccx.grch38_hprc.mc.config.json -r reads.gaf.gz -c calls.tsv -d diplotype.paths-stats.tsv -p diplotype.paths-info.tsv -o parakit.out.pdf
```

Other modes include

- *all_small* (a slightly more compact version of the *all* mode)
- *calls* just the results of variant calling and the reads supporting them
- *allele_support* just the coverage and allele ratio results
- *diplotype* just the diplotype prediction results

## Output

Example of a summary figure generated by Parakit.

![](example.summary.graph.jpg)

This sample shows strong evidence of a fusion and pathogenic SNV in a compound heterozygous configuration.

Points are positioned based on their position in the pangenome (x-axis). 
Because it's based on node position, some large regions are compressed.
This means that the x-axis is not exactly to scale with the genome sequence.

- Reads supporting pathogenic variants (black circle and triangle).
    - Here showing only one supporting read per variant (the longest read).
    - The variant position are highlighted by the black circle (SNV) and triangle (fusion).
    - The variants location are also marked by the vertical dotted lines in all panels, to help compare the different analyses.
- Copy number estimate from read coverage.
    - The dots shows the read coverage on each node across the region
    - The dashed line highlights the global copy number estimate based on this read coverage.
    - The dotted line is another global copy number estimate, often more robust to mapping bias. It is based on the reads around the boundaries and of the module, basically comparing looking at reads supporting extra copies (cycling back) and reads entering/leaving the module.
- Allele balance as the ratio of module 2 support. 
    - Expected around 0.5 if carrying two bimodular alleles (one module 1, one module 2). 
    - Deviation suggests fusions (or large gene conversion regions).
    - Faint blue lines highlight 1/3 and 2/3 marks, the expected levels in the presence of one fusion allele (and one bimodular allele).
- Diplotype candidate. Two haplotypes that match the reads best in term of alignment and coverage.
- Gene annotation

For the reads and diplotype panels, points are colored to highlight informative nodes (specific to module 1 in red or 2 in green).
Informative nodes/points are slightly shifted to help distinguish them.
The reads/haplotypes are split in parts when they loop back in the pangenome.

## Citation

For now, please cite [our preprint on medRxiv](https://www.medrxiv.org/content/10.1101/2025.02.07.25321404v1):

> Long-read sequencing resolves the clinically relevant CYP21A2 locus, supporting a new clinical test for Congenital Adrenal Hyperplasia. Jean Monlong, Xiao Chen, Hayk Barseghyan, William J Rowell, Shloka Negi, Natalie Nokoff, Lauren Mohnach, Josephine Hirsch, Courtney Finlayson, Catherine E. Keegan, Miguel Almalvez, Seth I. Berger, Ivan de Dios, Brandy McNulty, Alex Robertson, Karen H. Miga, Phyllis W. Speiser, Benedict Paten, Eric Vilain, Emmanuèle C. Délot. medRxiv 2025.02.07.25321404; doi: https://doi.org/10.1101/2025.02.07.25321404 

# FAQ

### I tried running Parakit and found a bug/issue. What should I do?

If it looks like it's an issue with external tools used by Parakit (GraphAligner, samtools, Minigraph-Cactus), make sure they are properly installed or, if using Docker, that docker is properly installed.
Another test could be to try running the `parakit` within the docker container that we've prepared: `quay.io/jmonlong/parakit:1.0.0`.
If that doesn't work, please post an [Issue](https://github.com/jmonlong/parakit/issues).

If the error looks like a Python error, you can try to rerun the command with a more verbose output.
Some commands can take a `-t` arguments to toggle a "debug trace" mode, which might give you more information about what Parakit is doing. 
Run `parakit COMMAND -h` to list options and check if there is a "debug trace mode" option.
Then, please post an [Issue](https://github.com/jmonlong/parakit/issues), pasting the output with this `-t` option.

### The diplotype output is inconsistent with the other predictions. Who's right?

When in doubt, believe the calls made directly from the reads (`call` subcommand), and the change in allele ratio (*allele_support* mode for the `viz` subcommand). 
The current diplotyping approach seems to work for samples in our study (~30x WGS of 30kb N50 reads with R10), but if reads are shorter or noisier the predicted diplotypes will most likely be incorrect. 
If your reads tend to be shorter than 30kb, take the diplotyping predictions with a grain of salt and focus on the other analysis to understand the composition of the RCCX module in the sample.

A more stable approach is under development using an EM algorithm that should convergence to good predictions across a wider range of input data. 

### How can I be sure that there is no functional CYP21A2 copy?

Two main recommendations:

First, the coverage and allele support graph can suggest the number of copies for each module.
Combined with the calls found in the reads, one can have a first idea if all the CYP21A2 copies are accounted for by the pathogenic variants called.
For example, if there is evidence for at least one CYP21A2 copy and no variants were identified in the reads, there is likely a functional copy.
If, like [in the example above](#output), the coverage and allele support suggest two CYP21A2 copies, and two mutually exclusive variants were found in the reads, there is no functional CYP21A2 copy.

Second, look at the graph made by `viz` with the *calls* mode. 
It will include reads that might support a functional CYP21A2 copy. 
These are reads that cover the variants identified but support the reference allele.
To make that figure, run something like:

```sh
parakit viz -v calls -j rccx.grch38_hprc.mc.config.json -r reads.gaf.gz -c calls.tsv -o parakit.calls.pdf
```

Note: to get quick stats about the reads aligning to the pangenome, run

```sh
parakit gafstats -r reads.gaf.gz -j rccx.grch38_hprc.mc.config.json
```

## Next

What we plan on the near future.

- [ ] Implement more robust diplotyping and add some information about our confidence in them
- [ ] Take long sequences, e.g. assembled contigs, as input.
- [ ] Test on Pacbio reads.
- [ ] Automate pangenome construction on other regions. 
