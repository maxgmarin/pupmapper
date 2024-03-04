#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Pupmapper: A tool for converting k-mer mappability to pileup mappability


import sys
import argparse
import os 
from ._version import __version__

import pandas as pd


def _pup_cli(args):
    ## 1) Set input parameters and PATHs ####
    input_KMap_BG = args.input

    output_PileupMap_BG = args.output

    kmer_size = args.kmer_size

    ## 2) Parse k-mer mappability bedgraph file as Pandas DF
    # ______

    ## 3) Convert kmap_DF to dict of numpy arrays of k-mer mappability

    ## 4) For each chromosome, calculate Pileup Mappability from k-mer mappability


    ## 5) Output the calculated pileup mappability values as a .bedgraph table

    print(f" Calculated pileup mappability scores (.bedgraph) output to: {output_PileupMap_BG}")
    


def main():
    pup_parser = argparse.ArgumentParser(description = "Tool for calculating genome wide pileup mappability based on k-mer mappability values")        

    pup_parser.add_argument('-i', '--input', type=str, required=True,
                                  help="Input k-mer mappability values in bedgraph format (.bedgraph). Ideally, generated with genmap software")

    pup_parser.add_argument('-o', '--output', type=str, required=True,
                                                    help="Output table of pileup mappability (.bedgraph)")

    pup_parser.add_argument('-k', '--kmer_size',type=int, required=True,
                                  help="k-mer size (bp) used to generate the input k-mer mappability values")
    
    pup_parser.set_defaults(func=_pup_cli)

    args = parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    sys.exit(main())
