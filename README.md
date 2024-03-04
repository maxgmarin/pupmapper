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
- [Basic usage](#basic-usage)
  - [test data set](#analyzing-included-test-data-set)
- [Full usage](#full-usage)
- [FAQ](#FAQ)


## Motivation

Pileup mappability is a useful metric for accessing sequence uniqueness within a genome and identifying repetive regions to remove from analysis when using short-read sequencing. It is calculated as the mean k-mer mappability of all k-mers overlapping a given position. This tool is meant to make it as easy as possible to quickly calculate pileup mappability from k-mer mappability values generated by the [Genmap](https://github.com/cpockrandt/genmap) software.


## Installation
### Install locally
`pupmapper` can be installed by cloning this repository and installing with `pip`.

```
git clone git@github.com:maxgmarin/pupmapper.git

cd pgqc

pip install . 
```

### `pip`
```
pip install pupmapper
```

## Basic usage
```
pupmapper -i kmap.K50E0.bedgraph -o pupmap.K50E0.bedgraph -k 50
```
The above command will run `pupmapper` on the input k-mer mappabilities (k= 50 bp, E = 0 mismatches) and output pileup mappability scores. 


### Analyzing included test sequence

If you wish to run an `panqc nrc` on a small test sequence, you can run the following commands:
>🚧 Check back soon 🚧
```
cd tests/data

pupmapper -i ${in_kmap} -o ${out_pupmap} -k 4
```

## Full usage
```
$ pupmapper --help

usage: pupmapper [-h] -i INPUT -o OUTPUT -k KMER_SIZE

Tool for calculating genome wide pileup mappability based on k-mer mappability values

optional arguments:
  -h, --help            show this help message and exit
  -i, --input INPUT
                        Input k-mer mappability values in bedgraph format (.bedgraph). Ideally,
                        generated with genmap software
  -o, --output OUTPUT
                        Output table of pileup mappability (.bedgraph)
  -k --kmer_size KMER_SIZE
                        k-mer size (bp) used to generate the input k-mer mappability values
....
```




## FAQ

### 1) How do I go from my genome of interest to identifying regions with low pileup mappability?

- 1.1) Use `genmap` (with your desired parameters) to calculate k-mer mappability for your genome of interest. (Output to .bedgraph)
- 1.2) Use `pupmapper` to calculate pileup mappability from k-mer mappability (output to .bedgraph)
- 1.3) Use `awk` and `bedtools` to identify regions of the genome which have pileup mappability < 1 (or below your desired threshold).


### 2) How do I generate the k-mer mappability values that pupmapper needs?
To calculate pileup mappability with `pupmapper` you must first generate k-mer mappability values with [genmap](https://github.com/cpockrandt/genmap). Refer to [getting started section](https://github.com/cpockrandt/genmap?tab=readme-ov-file#getting-started) of `genmap`'s README for more details.

You can use `genmap` in two steps:
#### 1) Index your target sequence
```
$ ./genmap index -F /path/to/fasta.fasta -I /path/to/index/folder
```

#### 2) Calculate k-mer mappability with desired parameters (in this example k = 30 bp and E = up to 2 mismatches)
```
$ ./genmap map -K 30 -E 2 -I /path/to/index/folder -O /path/to/output/folder -t -w -bg
```






