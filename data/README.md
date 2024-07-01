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
parakit construct -j pg.config.mc.cyp21a2.json -o rccx.grch38_hprc.mc
```

Will create two output files:

1. `rccx.grch38_hprc.mc.pg.gfa`
2. `rccx.grch38_hprc.mc.node_info.tsv`
