This repository contains the scripts used to construct synthetic long reads and the data to recreate the plots from Stapleton et al. 2015, "Haplotype-phased synthetic long reads from short-read sequencing."

Each of the samples described in the paper has its own directory,
Chicken/
Gelsemium2/
HepG2/
MG1655/
multiplex/
Env/
Gelsemium/
HCT116/
Potato/

... which contains a makefile that will construct synthetic reads from the raw data (which can be downloaded from the Sequence Read Archive). In practice, most of these assemblies are best done on a cluster. 

The python scripts called by the makefiles are in Scripts/

Figures.ipynb is an ipython notebook to recreate most of the figures in the paper.

HPCC_Hasher.sub and HPCC_Spades.sub are qsub files for running the hashing and assembly steps on a cluster.

This directory also contains JAStrim.fa and pairedTrim.txt, which are adapter files for trimmomatic.
