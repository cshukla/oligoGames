# oligoGames
A R package to analyze data from different massively parallel reporter assays

#### Why do I need another R package? :unamused: (aka "Motivation")

Since the early days of molecular biology, scientists have used reporter assays for different applications including studying the regulation of gene expression. While conventionally, one regulatory element is assayed at a given time, recent breakthroughs in microarray based nucleotide synthesis and next generation sequencing have allowed scientists to test several thousand elements in a single massively parallel reporter assay (MPRA). MPRA's were first demonstrated in [Patwardhan R.P. *et al.* (2009)](http://www.nature.com/nbt/journal/v27/n12/abs/nbt.1589.html) and [Melnikov A. *et al.* (2012)](http://www.nature.com/nbt/journal/v30/n3/full/nbt.2137.html) and the underlying technology has only become stronger and more robust in the following years. Recently, there has been a flurry of papers which have used MPRA's to test potential enhancer elements, eQTL's, SNP's as well as repressive elements - [Ulirsch J.C. *et al.* (2016)](http://www.cell.com/cell/fulltext/S0092-8674(16)30493-7), [Tewhey R. *et al* (2016)](http://www.cell.com/cell/fulltext/S0092-8674(16)30421-4), [Ernst J. *et al.* (2016)](http://www.nature.com/nbt/journal/v34/n11/full/nbt.3678.html). Although many "games" are being played by clever scientists to deduce regulatory elements using MPRAs, currently there is no dedicated software for analyzing MPRAs - which often have constraints a bit different from other *-Seq technologies. This is especially true for games which tile a long region of interest with short oligos and aim to identify regulatory elements contained within. Here, we provide the user functions for all aspects of MPRA experiments - from design of the oligo pool to learning kmers enriched in the infered regulatory elements.

# What is in the package

The following functions are provided in the package:

1. designOligoPool - Input a FASTA file of the regions you want to tile and get an output oligo arrary
2. mapToBarcodes - Map the FASTQ files from your experiment to the barcodes and generate a table of counts for each oligo
3. normCounts - Normalize oligo counts based on library size
4. modelNucCounts - Model nucleotide level counts from normalized oligo level counts
5. DRfinder - Perform inference to detect differential regions in oligo experiments

Please browse the vigenttes or read the documentation to understand the learn all the available options.

# Installation

Currently, oligoGames is under development and we recommend you contact us [here](https://github.com/cshukla/oligoGames) before using it. 

The easiest way to install the development version directly from Github is 

```r
devtools::install_github("oligoGames")
```

After installation, you can load the package into R.
```r
library(oligoGames)
```

Now you are ready to play some games :)

# Let's play a demo game

Here, we describe how to play a demo game. To begin playing, we first need to design an oligo pool.

We do this by using the designOligoPool function. We need a fasta file of sequences which are interested in tiling. Once we have that file, it is very simple to generate our oligo pool.

```r
regionsFile <- system.file('extdata', 'lncRNAs.fa' package='oligoGames')
outDir <- 'demoOligoGame'
designOligoPool(regionsFile)
```

At this stage, we need to order the oligos and do the experiment (likely several times because it never works the first time). Since this is a demo game, we are going to fast forward the 2.5 long years which we spent doing this and assume we have the sequencing data.

Next, we map the data to our oligo barcodes with the mapToBarcodes function.

```r
fastqCases <- c(system.file('extdata', 'fastqFiles', 'case1.fastq.gz', package='oligoGames'), system.file('extdata', 'fastqFiles', 'case2.fastq.gz', package='oligoGames'))

fastqControl <- c(system.file('extdata', 'fastqFiles', 'control1.fastq.gz', package='oligoGames'), system.file('extdata', 'fastqFiles', 'control2.fastq.gz', package='oligoGames'))

oligoMap <- system.file('extdata', 'lncLocOligoPool.fa' package='oligoGames')

mapToBarcodes(fastqCases, fastqControls, oligoMap, demoOligoGame)

```

Our next step is to normalize the counts for library size. We use the normCounts function for this.

```r
countsFile <- "demoOligoGame/allTranscriptsCounts_Raw.tsv"

normalizedCounts <- normCounts(countsFile, normType = "median")
```

We will now go from oligo counts to nucleotide level counts using modelNucCounts. 

```r
metaData <- system.file('extdata', 'oligoMeta.tab' package='oligoGames')

OligoSignal <- modelNucCounts(normalizedCounts, metaData, modelMethod = "median", oligoLen = 110)
```

Now we are in the final round of our game and we will infer the differential regions with the help of DRfinder

```r
diffRegions <- DRfinder(OligoSignal, conditionLabels = c("case", "control"))
```

# Future Work

1. Unit test for each function in the package.
2. Functions to analyze mutation MPRAs (QSAM, Mutual Information)

~~3. Robustify modelNucCounts
4. Add functionality to automatically generate metaData files when using designOligoPool~~

# Bug reports
Report bugs as issues on the [GitHub repository](https://github.com/cshukla/oligoGames)

# Contributors

* [Chinmay Shukla](https://github.com/cshukla)
* [Keegan Korthauer](https://github.com/kdkorthauer)
* [Rafael Irizarry](https://github.com/rafalab)

<a href="http://www.youtube.com/watch?feature=player_embedded&v=47jgoSJiB6w" target="_blank"><img src="http://img.youtube.com/vi/47jgoSJiB6w/0.jpg" alt="Immigrants (We Get The Job Done)" width="240" height="180" border="10" /></a>
