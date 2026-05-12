# Simulated samples

We simulate a couple of samples from the pangenome using `parakit sim` and test the pipeline.
The following samples are simulated:

### tsamp1

5 modules, 2 disrupted *gene* modules, one functional.

- *tsamp1_1* haplotype: *pseudogene*, *gene*, *gene* with 719T>A_Met240Lys and 1306C>T_Arg436Cys.
- *tsamp1_2* haplotype: *pseudogene*, fusion *pseudogene*-*gene*
- reads: 20-30kb long, ~20x
- [Simulation config file](sim.config.tsamp1.json)

### tsamp2

3 modules, 2 disrupted *gene* modules.

- *tsamp1_1* haplotpe: *pseudogene*, *gene* with 518T>A_Ile173Asn.
- *tsamp1_2* haplotype: fusion *pseudogene*-*gene*
- reads: 20-30kb long, ~20x
- [Simulation config file](sim.config.tsamp2.json)

### tsamp3

3 modules, 2 disrupted *gene* modules.

- *tsamp3_1* haplotpe: *pseudogene*, *gene* with 719T>A_Met240Lys and 1306C>T_Arg436Cys.
- *tsamp3_2* haplotype: fusion *pseudogene*-*gene* ~2kb upstream of *CYP21A2* and with 518T>A_Ile173Asn.
- reads: 20-30kb long, ~20x
- [Simulation config file](sim.config.tsamp3.json)

## Parakit pangenome

*Soon: how to setup the pangenome HPRC v2. Most likely copy from the ../data folder once it's up-to-date.*

## Snakemake workflow

```sh
## assuming the virtual environment where Parakit was installed is in ../env_parakit
source ../env_parakit/bin/activate

## run one sample (tsamp1)
snakemake --config sample=tsamp1,tsamp2,tsamp3 --cores 8 -p

## run all samples
snakemake --config sample=tsamp1,tsamp2,tsamp3 --cores 8 -p
```

## Commands for one sample

For information, if you want to test Parakit without using the Snakemake pipeline, here are the commands for one sample.

```sh
## list the variant in the pangenome
parakit deconstruct -j rccx.grch38_hprc2.mcc.config.json -o rccx.grch38_hprc2.mcc.decon.tsv

## create a directory where to place the different output files
mkdir -p tsamp1

## simulate haplotypes and reads
parakit sim -j rccx.grch38_hprc2.mcc.config.json -c sim.config.tsamp1.json -o tsamp1/tsamp1.sim.fa

## annotate the simulated diplotype
parakit annotate  -j rccx.grch38_hprc2.mcc.config.json -f tsamp1/tsamp1.sim.fa -o tsamp1/tsamp1.sim.annotate.pdf

## map reads and check some basic stats
parakit map -j rccx.grch38_hprc2.mcc.config.json -b tsamp1/tsamp1.sim.fq -o tsamp1/tsamp1.sim.gaf.gz
parakit gafstats -j rccx.grch38_hprc2.mcc.config.json -r tsamp1/tsamp1.sim.gaf.gz

## call variants from the aligned reads
parakit call -j rccx.grch38_hprc2.mcc.config.json -r tsamp1/tsamp1.sim.gaf.gz -s 3 -o tsamp1/tsamp1.sim.calls.tsv
parakit filtercalls -j rccx.grch38_hprc2.mcc.config.json -c tsamp1/tsamp1.sim.calls.tsv -af -o tsamp1/tsamp1.sim.filtered_calls.tsv

## estimate the module copy number
parakit copy -j rccx.grch38_hprc2.mcc.config.json -r tsamp1/tsamp1.sim.gaf.gz -o tsamp1/tsamp1.sim.copy.tsv

## reconstruct a diplotype
parakit diplotype -j rccx.grch38_hprc2.mcc.config.json -r tsamp1/tsamp1.sim.gaf.gz -o tsamp1/tsamp1.sim.diplotype

## visualize results
parakit viz -v all_small -j rccx.grch38_hprc2.mcc.config.json -r tsamp1/tsamp1.sim.gaf.gz -c tsamp1/tsamp1.sim.filtered_calls.tsv -m 1 -d tsamp1/tsamp1.sim.diplotype.paths-stats.tsv -p tsamp1/tsamp1.sim.diplotype.paths-info.tsv -o tsamp1/parakit.tsamp1.sim.pdf

## surject *gene* reads
parakit surject -j rccx.grch38_hprc2.mcc.config.json -r tsamp1.sim.gaf.gz -f tsamp1.sim.fq -o tsamp1.sim.parakit_surjected.bam

## surject reads and tag them with the predicted diplotype/module
parakit surject -j rccx.grch38_hprc2.mcc.config.json -r tsamp1.sim.gaf.gz -f tsamp1.sim.fq -d tsamp1.sim.diplotype.paths-stats.tsv -p tsamp1.sim.diplotype.paths-info.tsv -o tsamp1.sim.parakit_surjected_diplotype.bam
```

## Tradional mapping

In case we want to compare the surjected reads with the reads from a traditional mapping to a linear reference:

```sh
minimap2 -x map-ont -a -k 17 references/hg38.fa tsamp1/tsamp1.sim.fq -t 8 | samtools sort -O BAM > tsamp1/tsamp1.minimap2.bam
samtools index tsamp1/tsamp1.minimap2.bam
```
