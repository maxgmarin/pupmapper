<img width="500" src="https://github.com/maxgmarin/pupmapper/raw/main/Images/pupmapper.logo.png" alt="pupmapper logo">

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Static Badge](https://img.shields.io/badge/language-Python_3-blue)
<!---[![Build Status]()]()
[![github release version]()]()
[![DOI]()]()
--->

Pupmapper: A **P**ile**up** **Map**pability Calculato**r** 

<!---
> TBD Reference
--->

[TOC]: #
## Table of Contents
- [Motivation](#motivation)
- [Installation](#installation)
  - [Install locally](#install-locally)
  - [`pip`](#pip)
  - [`conda`](#conda)
- [Basic usage](#basic-usage)
  - [test data set](#analyzing-included-test-data-set)
- [Full usage](#full-usage)


## Motivation
The Pileup Mappability metric can be used to quickly identify regions which may be more difficult to perform variant calling with short-read WGS data. `pupmapper` was created to allow users to quickly convert k-mer mappability scores to pileup mappability. 

The first step of the `pupmapper` pipeline is to calculate k-mer uniquness scores using the [Genmap](https://github.com/cpockrandt/genmap) software. Then pupmapper will summarize the pileup mappability of each genomic position using the k-mer uniqueness of all overlapping k-mers.


### How is pileup mappability calculated from individual k-mer uniqueness/mappability scores?

<img width="500" src="https://github.com/maxgmarin/pupmapper/raw/main/Images/Derrien2012.PmapFig.png" alt="PmapFig">
The Pileup mappability of a position is specifically calculated as the mean k-mer mappability of all k-mers overlapping a given position. <br>

A pileup mappability score of 1 indicates that all k-mers overlapping with a position are unique within the genome (using the user defined parameters of uniqueness).

**Simply put, Pileup mappability is useful because it gives a sense of uniquemess of all possible reads (of defined length) that could arise from given genomic position.**


### Useful reading for k-mer mappability and pileup mappability:

> Derrien, T, (2012). Fast Computation and Applications of Genome Mappability. PLOS ONE 7(1): e30377. [https://doi.org/10.1371/journal.pone.0030377](https://doi.org/10.1371/journal.pone.0030377)

> Pockrandt C, (2020) GenMap: ultra-fast computation of genome mappability, Bioinformatics, Volume 36, Issue 12, June 2020, Pages 3687–3692, [https://doi.org/10.1093/bioinformatics/btaa222](https://doi.org/10.1093/bioinformatics/btaa222)

> Lee H, Schatz MC. (2012). Genomic dark matter: the reliability of short read mapping illustrated by the genome mappability score, Bioinformatics, Volume 28, Issue 16, August 2012, Pages 2097–2105, [https://doi.org/10.1093/bioinformatics/bts330](https://doi.org/10.1093/bioinformatics/bts330)

## Installation

You will need to install both the `pupmapper` python package and ensure that the [genmap](https://github.com/cpockrandt/genmap) software is installed (available on your $PATH environmental variable.)
If conda is used for installation software dependencies are installed automatically. If pupmapper is installed via `pip`, then [genmap](https://github.com/cpockrandt/genmap) and [bigtools]([https://github.com/cpockrandt/genmap](https://github.com/jackh726/bigtools)) need to be installed independently available on your $PATH environmental variable.).

### `conda`
Installing pupmapper via conda will install bpoth the pupmapper package and needed software dependencies ([genmap](https://github.com/cpockrandt/genmap) and [bigtools]([https://github.com/cpockrandt/genmap](https://github.com/jackh726/bigtools)) in one command.

To install pupmapper to your current environment:
```
conda install -c bioconda pupmapper
```

To create a dedicated conda environment with pupmapper installed:
```
conda create -n pupmapper_env -c bioconda pupmapper
```

### `pip`
pupmapper can be installed via `pip` python package manager.
```
pip install pupmapper
```

### Install locally
`pupmapper` can be installed by cloning this repository and installing with `pip`.

```
git clone git@github.com:maxgmarin/pupmapper.git

cd pupmapper

pip install . 
```

## Basic usage

#### 1) `pupmapper all` - Run the full pipeline starting with an input genome

```
pupmapper all -i Input.Genome.fasta -o output_directory/ -k 50 -e 1
```
The above command will first use genmap to calculate k-mer mappability scores for the input genome and then calculate pileup mappability scores.

Arguments:

`-i, --in_genome_fa`: Input genome FASTA file. <br>
`-o, --outdir`: Directory for output files. <br>
`-k, --kmer_len`: K-mer length (e.g., 50 bp). <br>
`-e, --errors`: Number of allowed mismatches (hamming distance) in k-mer mapping step. <br>
`-g, --gff`: (Optional) Input genome annotations in GFF format. <br>
`--save-numpy`: (Optional) Save results as compressed numpy arrays. <br>


### Analyzing included test sequence

If you wish to run an `pupmapper` on a small test sequence (15 bp), you can run the following commands:
```
cd tests/data/Genmap_Ex1/

pupmapper all -i Ex1.genome.fasta -o Ex1_OutputDir -k 4 -e 0
```
This command will analyze the pileup mappability of the test sequence with a k-mer size of 4 bp and a max mismatch of 0 (K=4,E=0).


## Full usage
```
pupmapper all --help
usage: pupmapper run_all [-h] -i IN_GENOME_FA -o OUTDIR -k KMER_LEN -e ERRORS [-g GFF] [--save-numpy]

optional arguments:
  -h, --help            show this help message and exit
  -i IN_GENOME_FA, --in_genome_fa IN_GENOME_FA
                        Input genome fasta file (.fasta)
  -o OUTDIR, --outdir OUTDIR
                        Directory for all outputs of k-mer and pileup mappability processing.
  -k KMER_LEN, --kmer_len KMER_LEN
                        k-mer length (bp) used to generate the k-mer mappability values
  -e ERRORS, --errors ERRORS
                        Number of errors (mismatches) allowed in Genmap's k-mer mappability calculation
  -g GFF, --gff GFF     GFF formatted genome annotations for input genome (.gff) (Optional)
  --save-numpy          If enabled, all pileup mappability scores will be output as compressed numpy arrays (.npz).
```






