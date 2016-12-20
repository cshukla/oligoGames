# oligoGames
A R package to analyze data from different massively parallel reporter assays

## Motivation

Since the early days of molecular biology, scientists have used reporter assays for different applications including studying the regulation of gene expression. While conventionally, one regulatory element is assayed at a given time, recent breakthroughs in microarray based nucleotide synthesis and next generation sequencing have allowed scientists to test several thousand elements in a single experiment. Such massively parallel reporter assays were first demonstrated in [Patwardhan R.P. *et al.* (2009)](http://www.nature.com/nbt/journal/v27/n12/abs/nbt.1589.html) and [Melnikov A. *et al.* (2012)](http://www.nature.com/nbt/journal/v30/n3/full/nbt.2137.html). Recently, there has been a flurry of papers which have used MPRA's to test potential enhancer elements, eQTL's, SNP's as well as repressive elements [[Ulirsch J.C. *et al.* (2016)](http://www.cell.com/cell/fulltext/S0092-8674(16)30493-7), [Tewhey R. *et al* (2016)](http://www.cell.com/cell/fulltext/S0092-8674(16)30421-4), [Ernst J. *et al.* (2016)](http://www.nature.com/nbt/journal/v34/n11/full/nbt.3678.html)]. Although many "games" are being played by clever scientists to deduce regulatory elements using MPRAs, currently there is no software for analyzing MPRA data. This is especially true for experiments which tile a long region of interest with short oligos to identify regulatory elements contained within.

Here, we provide the user functions for all aspects of MPRA experiments - from design of the oligo pool to learning kmers enriched in the infered regulatory elements.

## What is in the package

The following functions are provided in the package:

1. generateOligos - Input a FASTA file of the regions you want to tile and get an output oligo arrary
2. mapToBarcodes - Map the FASTQ files from your experiment to the barcodes and generate a table of counts for each oligo
3. normCounts - Normalize oligo counts for library size
4. modelNucCounts - Model nucleotide level counts from normalized oligo level counts

Please browse the vigenttes or read the documentation to understand the learn all the available options.
