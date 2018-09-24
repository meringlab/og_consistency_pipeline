# Consistency pipeline for hierarchies of orthologous groups

This repository contains the python implementation for the methodology described in:

> Heller, D., Szklarczyk, D. and von Mering, C.: Tree reconciliation combined with subsampling improves large scale inference of orthologous group hierarchies (2018) manuscript in preparation

A preprint of the article can be found on bioRxiv at [https://doi.org/10.1101/417840](https://doi.org/10.1101/417840)

---

# Content

The current version of the pipeline (v0.1) is a first stand-alone version written in python2, which relies on the python tree library [etetoolkit](http://etetoolkit.org), as well as on the following software to compute and reconcile gene trees with species trees:

- [MAFFT](https://mafft.cbrc.jp/alignment/software/linuxportable.html) for multiple sequence alignment
- [FastTree](http://www.microbesonline.org/fasttree/#Install) for tree prediction
- [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/) for tree reconciliation

For convenience the binaries of the three tools have been included in the archive `bin.tar.gz`. We also added an example dataset from the [eggNOG database](http://eggnog.embl.de) in `data.tar.gz`. Currently the package structure is tied to how the eggNOG data is organized, but we plan to update the software in the coming weeks to a [Snakemake](https://snakemake.readthedocs.io/en/stable/) based version, for easier application to other datasets.

# Example execution

By running the `setup.py` script, the two archives are extracted and the small eggNOG related example dataset is generated. The data directory contains information regarding the Primates level of eggNOG and its two sublevels, Hominidae and Cercopithecoidea:

```
                             /-314294[prNOG-1][superfamily:Cercopithecoidea]
-9443[prNOG][order:Primates]--
                             \-9604[homNOG][family:Hominidae]
```

For the 15 member species of the Primates level (see `data/9443.primates.species.tsv`), the data directory includes FASTA sequences (in `data/fastafiles`) and orthologous group mappings (in `data/pickles`) as well as clades (in `data/clades`). 

- To test the consistency pipeline on the small example data ([test_data] directory), execute `test.sh`, which will make the three levels hierarchically consistent. The command in the script is limited to a small number of trees, to run the full example (~ 8000 Trees), an execution with several cores is recommended (`-c` option, see final comment `test.sh`)

# Contact

Feedback is always welcome. Feel free to write to davide.heller@imls.uzh.ch
