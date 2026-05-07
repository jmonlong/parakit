# `parakit deconstruct`

*Soon.*

# `parakit sim`

```sh
usage: parakit sim [-h] -c C [-j J] [-g G] [-n N] [-a A] [-o O] [-t]

options:
  -h, --help  show this help message and exit
  -c C        input simulation configuration
  -j J        config JSON file
  -g G        input GFA pangenome (optional)
  -n N        node information (optional)
  -a A        annotation file (e.g. from ClinVar) (optional)
  -o O        output fasta file
  -t          debug trace mode
```

- `-c` JSON file with information on what/how to simulate (details below)
- `-j` JSON configuration for this locus/pangenome, same as used in most commands. 
    - Helps find/guess many helper files.
- `-o` output fasta file with the two haplotype sequences. 
    - A FASTQ might also be written if reads are simulated
- Optional files that should be guessed using config file (`-j`).
    - `-g` GFA file of the pangenome. 
    - `-j` TSV file with node information. 
    - `-a` TSV file with variant annotation. 
    
Minimal usage example:

```sh
parakit sim -j rccx.grch38_hprc2.mcc.config.json -c sim.config.json -o sim.diplotype.fa
```

The simulation configuration (`-c`) is a JSON file with information about the two haplotypes to simulate, and, eventually how to simulate reads.
It has a *haplotypes* Object with two elements (the haplotypes). 
Each haplotype is an array with the modules to simulate, in order. 
A module is an Object with:

- `start_mod`: the type of module to start on, either *c1* or *c2*.
- `mod_noise`: proportion of nodes specific to the other module to include. In real haplotype, even from healthy individuals, a module is never 100% made out of nodes specific to that module type.
- `alt_paths`: if we want to inject specific variants, define an Object associating node start with the alternate path. Extract this from the deconstructed output, for example.
- `fusions`: if we want to simulate a fusion, which node should the module type change. It's recommended to pick a conserved node, i.e. that most/all haplotypes have.

Finally, in addition to the `haplotypes` object, we can define a `read_length` (array of two elements defining a range) and `haplotype_depth` (read depth for each haplotype).
If this information is present, reads will be simulated and a FASTQ will be generated.
Warning: the read simulation is quite simplistic with just random substitution and indel errors, and is used mostly to test the installation or give preliminary insighths on the effect of some parameters (read length, depth).

For example, for a trimodular haplotype and a bimodular haplotype, with two specific short variants and a fusion, and simulating 20-30kb reads for a total coverage of 20x:

```json
{"haplotypes": 
 {
   "tsamp1_1": [
     {"start_mod": "c1", "mod_noise": 0.05},
     {"start_mod": "c2", "mod_noise": 0.01},
     {"start_mod": "c2", "mod_noise": 0.05,
      "alt_paths": {"8547": "8547_8548_8550",
                    "8664": "8664_8666_8667"}}
   ],
   "tsamp1_2": [
     {"start_mod": "c1", "mod_noise": 0.05},
     {"start_mod": "c1", "mod_noise": 0.05, "fusions": ["8333"]}
   ]
 },
 "read_length": [20000, 30000],
 "haplotype_depth": 10
}
```

### How to find which nodes to use with `alt_paths` or `fusions`?

To pick the specific variants to inject, we can look for ClinVar variants in the deconstructed pangenome (i.e. the list of variants in the pangenome, see [*deconstruct* command](#parakit_deconstruct)).

```sh
$ grep -e Met240Lys -e Arg436Cys rccx.grch38_hprc2.mcc.decon.tsv
8547_8548_8550	32039816	32039817	8547	8547_8549_8550	T	8547_8548_8550	A	c2	40032_719T>A_Met240Lys
8664_8666_8667	32040952	32040953	8664	8664_8665_8667	C	8664_8666_8667	T	c2	190758_1306C>T_Arg436Cys
```

For the fusion, let's say we want it to happen after the variant at position 32,038,419 in the genome.
We find the corresponding node in the variant list, then double-check that the node is indeed conserved enough:

```sh
$ grep 32038419 rccx.grch38_hprc2.mcc.decon.tsv
8330_8331_8333	32038419	32038420	8330	8330_8332_8333	C	8330_8331_8333	T	c2	None
$ grep -C 3 -w 8333 rccx.grch38_hprc2.mcc.node_info.tsv
8330	8	GGGCGTCT	2	449	479	325144	357878	8330	none
8331	1	T	1	414	13	325152	357886	8330	c1
8332	1	C	1	35	466	325152	357886	8330	c2
8333	17	GCCATGCTGC+	2	449	479	325153	357887	8333	none
8334	1	C	2	449	478	325170	357904	8334	none
8335	3	CTG	2	449	385	325171	357905	8335	none
8336	14	CTGCTGCTGC+	2	449	478	325174	357908	8336	none
```

Node `8333` is good because it has class *none* (not specific to a module type), and a high number of reference modules c1/c2 (449/479) that go through there. 
That means, it's very likely the simulated haplotype will also pass through there which will make it possible the fusion to be included.

> Note: To simulate a large gene conversion, one could specify two `fusions` breakpoints.

Find examples of configurations and commands in the [`test` directory of the Parakit repo](https://github.com/jmonlong/parakit/test).
