## Provided files

- `rccx.grch38_hprc.mc.config.json` the configuration file for this pangenome (contains coordinates, flank size, etc used to build the pangenome)
- `rccx.grch38_hprc.mc.pg.gfa` the pangenome in GFA format
- `rccx.grch38_hprc.mc.node_info.tsv` metadata about the nodes in the pangenome, e.g. which one is specific to module 1/2.
- Annotation files:
    - `CYP21A2.pathogenic.variant_summary.2024_09_03.tsv` reformatted subset of ClinVar including CYP21A2 pathogenic variants
    - `CYP21A2.gencodev43.nearby_genes.tsv` reformatted subset of GENCODE containing gene annotation in the region.

See [rccx.summary.md](rccx.summary.md) for some descriptive metrics on this pangenome.

## Pangenome construction

### Download HPRC sequence and install AGC

AGC (https://github.com/refresh-bio/agc) will be used to extract the sequence of interest during pangenome construction. 
It needs to be installed (easiest is to download a prebuilt binary from the [release page](https://github.com/refresh-bio/agc/releases)).

The AGC file with the HPRC assemblies will also be needed for pangenome construction.
Download it with:

```
wget -O HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1
```

### Prepare/download HPRC annotation files

The paths to the Ensembl annotations was compiled in [hprc.ensembl.gff3.paths.tsv](hprc.ensembl.gff3.paths.tsv) (with [list.ensembl.paths.R](list.ensembl.paths.R)).

For the CAT and segmental duplications, download the index with:

```sh
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_CAT_genes.index
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_Seg_Dups.index
```

Then, run the Snakemake pipeline to download each file and subset records about `CYP21A*`, `TNXB`, `C4A`/`C4B`
For example:

```sh
snakemake --cores 8
```

### Method 1 - Specifying the position of each module in all input assemblies

The first method developed.

- Pros
    - Complete control on what is considered module 1/2
    - Pangenome is simpler (flat with one loop)
- Cons
    - Works only for tandem duplications (no buffer regions handled yet)
    - Requires trustworthy annotation
    - More "manual" work
    - Can't integrate as much sequence as not all sequence are clearly

#### Find position of each module in the HPRC assemblies

Using the downloaded HPRC annotations, the coordinates for the RCCX modules for each assembly were extracted by identifying haplotypes where one gene and one pseudogene were confidently annotated.
See [extract-rccx-coords.md report](extract-rccx-coords.md) for details.
It creates the [`hprc.cyp21a.coords.tsv`](hprc.cyp21a.coords.tsv) file.
It's a TSV with no header containing the coordinate of each module/copy for the HPRC haplotype to extract.
The three column in the file are: sample name, contig name, label (*c1_\** for a module/copy 1 allele, *c2_\** for a module/copy 2 allele).

#### Config file

Create a config file, see [`rccx.grch38_hprc.mc.config.json`](rccx.grch38_hprc.mc.config.json) file:

```json
{
    "ref_fa": "PATH/TO/hg38.fa"
    "c1": "chr6:31980532-32013273",
    "c2": "chr6:32013273-32046127",
    "flank_size": 300000,
    "method": "mc",
    "hprc_agc": "HPRC-yr1.agc"
    "hprc_coords": "hprc.cyp21a.coords.tsv"
}
```

### Method 2 - Specifying the region of interest in all input assemblies

Newer method that lets the pangenome builder do the collapse. 
We guess the modules after building the pangenome

- Pros
    - More automated process.
    - No/less need of an annotation
    - Works in the presence of buffer region between the modules
- Cons
    - Pangenome might be more complex
    - Not tested as much

#### Find position of each module in the HPRC assemblies

Using the downloaded HPRC annotations, we want to find the contig and coordinates of the region of interest in each HPRC assembly.
See [extract-rccx-region-coords.md report](extract-rccx-region-coords.md) for details.
It creates the [hprc.cyp21a.haps.coords.tsv](hprc.cyp21a.haps.coords.tsv) specified in the config below.
It's a TSV with no header containing the coordinate of the region of interest to extract.
The two columns in the file are: sample name and coordinates.

#### Config file

Create a config file, see [`rccx.grch38_hprc.mcc.config.json`](rccx.grch38_hprc.mcc.config.json) file.
Short vesion:

```json
{
    "ref_fa": "PATH/TO/hg38.fa"
    "c1": "chr6:31980532-32013273",
    "c2": "chr6:32013273-32046127",
    "flank_size": 300000,
    "method": "mc_collapse",
    "hprc_agc": "HPRC-yr1.agc"
    "hprc_coords": "hprc.cyp21a.haps.coords.tsv"
}
```


### Run

```sh
parakit construct -j config.json
```

Will create two output files:

1. A GFA file with the pangenome (`.gfa`)
2. A TSV with node information (`.node_info.tsv`)

### QC - Visualize the pangenome

- Open the GFA in [Bandage](https://rrwick.github.io/Bandage/)
- Compute stats about the nodes classifications using something like [rccx.summary.md](rccx.summary.md) (source [rccx.summary.Rmd](rccx.summary.Rmd))
- Look at the module classification (c1/c2) for all input sequences like in the exploratory analysis in [guess.rccx.modules.pdf](guess.rccx.modules.pdf) (source [guess.rccx.modules.Rmd](guess.rccx.modules.Rmd)).


### Other annotations 

Other annotations, which can be included in the config file or specified later, are used for some analysis or visualization.
For now they include:

1. Pathogenic variants in ClinVar in the gene of interest.
2. Gene/exon annotation from GENCODE.

For example, below are some commands to prepare:

- `CYP21A2.pathogenic.variant_summary.2024_09_03.tsv`
- `CYP21A2.gencodev43.nearby_genes.tsv`

First download ClinVar variants:

```sh
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
date +%Y_%m_%d > variant_summary.txt.date
```

We also keep track on the date it was download in `variant_summary.txt.date` because this file is updated continuously by ClinVar.

Then run the `prepare-annotations.R` R script, for example with:

```sh
Rscript prepare-annotations.R
```

This script will read the ClinVar variant file and extract information (position, protein change, ...) for the ones of interest, i.e. affecting the gene of interest and potentially pathogenic.
It will also download the GENCODE annotation (if needed) and extract the position of the genes, exons, UTRs, coding sequences, in the region of interest. 
Adapt the script if you're interested in a different region/gene.


### Building a pangenome for another region

Currently under development. 
In theory, one would need to:

1. Update the [Snakefile](Snakefile) that downloads HPRC annotation. Change the gene names of interest.
2. Edit the config JSON with the coordinates of the two new regions of interest (*c1* and *c2* fields).
3. Make a new TSV file with the coordinates of the region of interest in the HPRC assemblies (*hprc_coords* in the config JSON)
4. Prepare the other annotations (ClinVar/Gencode) by changing the names of the genes of interest in the [prepare-annotations.R](prepare-annotations.R) script (replace *CYP21A2*).
