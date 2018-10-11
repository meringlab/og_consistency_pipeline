# Consistency pipeline for hierarchies of orthologous groups

This repository contains the python implementation for the methodology described in:

> Heller, D., Szklarczyk, D. and von Mering, C.: Tree reconciliation combined with subsampling improves large scale inference of orthologous group hierarchies (2018) manuscript in preparation

A preprint of the article can be found on bioRxiv at [https://doi.org/10.1101/417840](https://doi.org/10.1101/417840)

---

# Content

The current version of the pipeline (v0.4) is a [Snakemake](https://snakemake.readthedocs.io/en/stable) workflow written in python3, which relies on the python tree library [etetoolkit](http://etetoolkit.org) and the progress bar [tqdm](https://github.com/tqdm/tqdm). By default the following software is used to compute and reconcile gene trees with species trees:

- [MAFFT](https://mafft.cbrc.jp/alignment/software/linuxportable.html) for multiple sequence alignment
- [FastTree](http://www.microbesonline.org/fasttree/#Install) for tree prediction
- [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/) for tree reconciliation

The binaries of the three tools are downloaded automatically using the snakemake rules specified in `rules/tools.smk`. 

Input files are specified through the configuration file `config.yaml`, with parameters explained therein. As a small example we included a dataset from the [eggNOG database](http://eggnog.embl.de) in `data.tar.gz`.

The software has been developed and tested on Linux (Ubuntu 12/16/18.04). Other Unix systems might be suitable as well but binaries will have to be adapted accordingly.

# Installation

The easiest way to use the pipeline is to create a python3 environment with the [Anaconda/Miniconda](https://www.anaconda.com) distribution (installation instructions [here](https://conda.io/docs/user-guide/install/index.html)). Assuming that the distrution has been installed, the following commands create a new environment and install all the required dependencies:

```
# create a new environment named "smk"
conda create -n smk python=3.6
# activate the environment
source activate smk
# install the dependencies (snakemake, ete3, tqdm)
conda install -c bioconda -c conda-forge snakemake
conda install -c etetoolkit ete3 ete_toolchain 
conda install -c conda-forge tqdm
```

Alternatively the dependencies can also be installed natively using pip or compiled from source by following the respective guides in their documentation.

# Example execution

The configuration file `config.yaml` is predefined with the input parameters for the small example included in `data.tar.gz`. The archive contains information regarding the Primates level of eggNOG and its two sublevels, Hominidae and Cercopithecoidea:

```
                             /-314294[prNOG-1][superfamily:Cercopithecoidea]
-9443[prNOG][order:Primates]--
                             \-9604[homNOG][family:Hominidae]
```

For the 15 member species of the Primates level (see `data/9443.primates.species.tsv`), the data directory includes FASTA sequences (in `data/fastafiles`) and orthologous group mappings (in `data/orthologous_groups`) as well as the clades (in `data/clades`). 

To run the Snakemake workflow: 

1. expand the example dataset with `tar -xzf data.tar.gz` 
2. (opt) list the outstanding tasks with `snakemake -n` or `snakemake --dag | dot -Tsvg > dag.svg` to visualize them as SVG graph
3. execute the tasks with `snakemake`

The software will read the test dataset with 100 OGs from `data/orthologous_groups` and resolve the hierarchical inconsistencies. After workflow completion (~2 min on a single core) the consistent OG definition can be found in `test_output/consistent_ogs`. To run a larger example with the complete clustering of the 15 species, change the input parameter in the `config.yaml` file to point at `data/orthologous_groups_full`. Be aware that this will require much more time and multi-core execution is strongly reccomended (~1h using 10 cores, i.e. `snakemake --cores 10`).

# Contact

Feedback is always welcome. Feel free to write to davide.heller@imls.uzh.ch
