## Provided files

- `rccx.grch38_hprc.mc.config.json` the configuration file for this pangenome (contains coordinates, flank size, etc used to build the pangenome)
- `rccx.grch38_hprc.mc.pg.gfa` the pangenome in GFA format
- `rccx.grch38_hprc.mc.node_info.tsv` metadata about the nodes in the pangenome, e.g. which one is specific to module 1/2.
- Annotation files:
    - `CYP21A2.pathogenic.variant_summary.20231127.txt` reformatted subset of ClinVar including CYP21A2 pathogenic variants
    - `CYP21A2.gencodev43.nearby_genes.tsv` reformatted subset of GENCODE containing gene annotation in the region.

See [rccx.mc.summary.md](rccx.mc.summary.md) for some descriptive metrics on this pangenome.

## Pangenome construction

Semi-automated approach. 

### Additional dependencies

- agc https://github.com/refresh-bio/agc


### Find position of each module in the HPRC assemblies

First, prepare/download files with the path to the annotation of each assembly.

```sh
Rscript list.ensembl.paths.R  # creates hprc.ensembl.gff3.paths.tsv
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_CAT_genes.index
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_Seg_Dups.index
```

Run the Snakemake pipeline that download each file and subset them to the record about `CYP21A*`, `TNXB`, `C4A`/`C4B`

```sh
snakemake --cores 8
```

Then the coordinates for the RCCX modules for each assembly were extracted by identifying haplotypes where one gene and one pseudogene were confidently annotated.
See [extract-rccx-coords.md report](extract-rccx-coords.md) for details.

### Config file

Create a `pg.config.mc.cyp21a2.json` file:

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

The `HPRC-yr1.agc` file contains the sequence from the HPRC dataset.
Download it with:

```
wget -O HPRC-yr1.agc https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1
```

The `hprc.cyp21a.coords.tsv` is a TSV with no header containing the coordinate of each module/copy for the HPRC haplotype to extract.
The three column in the file are: sample name, contig name, label (*c1_\** for a module/copy 1 allele, *c2_\** for a module/copy 2 allele).
Making this file is not automated yet.

### Run

```sh
parakit construct -j rccx.grch38_hprc.mc.config.json -o rccx.grch38_hprc.mc
```

Will create two output files:

1. `rccx.grch38_hprc.mc.pg.gfa`
2. `rccx.grch38_hprc.mc.node_info.tsv`

### Annotations 

To prepare:

- `CYP21A2.pathogenic.variant_summary.20231127.txt`
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


