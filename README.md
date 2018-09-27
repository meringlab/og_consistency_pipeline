# Consistency pipeline for hierarchies of orthologous groups

This repository contains the python implementation for the methodology described in:

> Heller, D., Szklarczyk, D. and von Mering, C.: Tree reconciliation combined with subsampling improves large scale inference of orthologous group hierarchies (2018) manuscript in preparation

A preprint of the article can be found on bioRxiv at [https://doi.org/10.1101/417840](https://doi.org/10.1101/417840)

---

# Content

The current version of the pipeline (v0.3) is a [Snakemake](https://snakemake.readthedocs.io/en/stable) workflow written in python3, which relies on the python tree library [etetoolkit](http://etetoolkit.org) and the progress bar [tqdm](https://github.com/tqdm/tqdm), as well as on the following software to compute and reconcile gene trees with species trees:

- [MAFFT](https://mafft.cbrc.jp/alignment/software/linuxportable.html) for multiple sequence alignment
- [FastTree](http://www.microbesonline.org/fasttree/#Install) for tree prediction
- [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/) for tree reconciliation

For convenience the binaries of the three tools have been included in the archive `bin.tar.gz`. We also added an example dataset from the [eggNOG database](http://eggnog.embl.de) in `data.tar.gz`. Currently the package structure is tied to how the eggNOG data is organized, but we plan to update the software in the coming weeks for easier application to other datasets.

The software has only been tested on Linux (Ubuntu 12/16/18.04). Other Unix systems might be suitable as well but binaries have to be adapted accordingly.

# Example execution

The configuration file `config.yaml` contains the parameters to run the pipeline on a small example generated from the included eggNOG dataset. The `data.tar.gz` archive contains information regarding the Primates level of eggNOG and its two sublevels, Hominidae and Cercopithecoidea:

```
                             /-314294[prNOG-1][superfamily:Cercopithecoidea]
-9443[prNOG][order:Primates]--
                             \-9604[homNOG][family:Hominidae]
```

For the 15 member species of the Primates level (see `data/9443.primates.species.tsv`), the data directory includes FASTA sequences (in `data/fastafiles`) and orthologous group mappings (in `data/pickles`) as well as clades (in `data/clades`). 

- To run the Snakemake workflow, simply type `snakemake`. This will generate a small inconsistent OG definition in `test_data`. After workflow completion the consistent results can be found in `test_output/consistent_ogs`.

# Contact

Feedback is always welcome. Feel free to write to davide.heller@imls.uzh.ch
