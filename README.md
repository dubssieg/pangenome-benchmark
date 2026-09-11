# pangenome-benchmark

# Requirements

+ [MSpangenome](forge.inrae.fr/pangepop/MSpangepop)
+ PGGB
+ Minigraph-Cactus
+ samtools
+ vg toolkit
+ [pancat](github.com/dubssieg/pancat)

# Steps

This bash pipeline:

1. simulates from a single fasta sequence a population using MSpangenome
2. builds MC and PGGB graphs from the simulated pangenome
3. converts files in similar formats (haplotype names, encoding)
4. performs comparisons between simulated and reconstructed variation graphs

Results can be visualized using the joint notebooks (you may have to adapt the code to your setup).

# Data availability

MSpangenome parameters and genomes are available on [Zenodo](https://doi.org/10.5281/zenodo.22685598).