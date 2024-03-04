#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Pupmapper: A tool for converting k-mer mappability to pileup mappability


import sys
import argparse
import os 
from ._version import __version__

import pandas as pd
import numpy as np
from tqdm import tqdm

from .utils import kmap_bedgraph_to_DF, kmap_DF_To_ArrayDict
from .utils import convert_kmap_to_pmap_arrays
from .utils import process_nparrays_to_bedgraph_df


def _pup_cli(args):
    ## 1) Set input parameters and PATHs ####
    input_KMap_BG = args.input

    output_PileupMap_BG = args.output

    kmer_size = args.kmer_size

    ## 2) Parse k-mer mappability bedgraph file as Pandas DF
    Kmap_DF = kmap_bedgraph_to_DF(input_KMap_BG)

    ## 3) Convert kmap_DF to dict of numpy arrays of k-mer mappability
    Kmap_Arrays = kmap_DF_To_ArrayDict(Kmap_DF)

    ## 4) For each chromosome, calculate Pileup Mappability from k-mer mappability
    Pmap_Arrays = convert_kmap_to_pmap_arrays(Kmap_Arrays, kmer_size)

    ## 5) Output the calculated pileup mappability values as a .bedgraph table

    Pmap_BEDGRAPH_DF = process_nparrays_to_bedgraph_df(Pmap_Arrays)

    Pmap_BEDGRAPH_DF.to_csv(output_PileupMap_BG,
                          sep = "\t", index = False)

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

    args = pup_parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    sys.exit(main())
