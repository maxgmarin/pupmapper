#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Pupmapper: A tool for converting k-mer mappability to pileup mappability


import sys
import argparse
import os
import shutil

from ._version import __version__

from .utils import kmap_bedgraph_to_DF, kmap_DF_To_ArrayDict
from .utils import convert_kmap_to_pmap_arrays
from .utils import process_nparrays_to_bedgraph_df, save_numpy_array_dict
from .utils import getRegions_BelowThreshold_PupMap
from .utils import calc_pupmap_per_gene

from .genmap import get_fasta_basename
from .genmap import is_genmap_available, run_genmap_index, run_genmap_map




def _genmap_2steps_cli(args):

    # Step 1: Set input parameters and define variables for relevant paths
    input_Genome_FA = args.in_genome_fa 

    kmer_length = args.kmer_len

    errors = args.errors

    index_dir = args.outdir + "/genmap_index"

    kmap_dir = args.outdir + "/genmap_kmap_K" + str(kmer_length) + "_E" + str(errors)


    FA_BaseName = get_fasta_basename(input_Genome_FA)

    genmap_out_prefix = kmap_dir + "/" + FA_BaseName + ".kmermap." + "K" + str(kmer_length) + "_E" + str(errors)

    # Create the target output directory
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(kmap_dir, exist_ok=True)

    # Step 2: Check if Genmap is available in the path
    if not is_genmap_available():
        raise EnvironmentError("genmap is not available on the PATH. Please install it or fix your PATH.")

    # Step 3: Run Genmap Index step for input genome

    run_genmap_index(input_Genome_FA,
                     index_dir)

    # Step 4: Run Genmap MAP step for genome (after indexing) to calculate k-mer mappability

    run_genmap_map(index_dir,
                   genmap_out_prefix,
                   kmer_length,
                   errors)
    
    print(f"Genmap processing completed.\nOutput files saved to: {args.outdir}")
    print(f"Indexed genome files saved to: {index_dir}")
    print(f"K-mer mappability files saved to: {kmap_dir}")
    print("")





def _pup_cli(args):
    ## 1) Set input parameters and PATHs ####
    input_KMap_BG = args.input

    output_PileupMap_BG = args.output

    kmer_len = args.kmer_len

    input_Genome_GFF = args.gff 

    save_npz = args.save_numpy


    ## 2) Parse k-mer mappability bedgraph file as Pandas DF
    Kmap_DF = kmap_bedgraph_to_DF(input_KMap_BG)

    ## 3) Convert kmap_DF to dict of numpy arrays of k-mer mappability
    Kmap_Arrays = kmap_DF_To_ArrayDict(Kmap_DF)

    ## 4) For each chromosome, calculate Pileup Mappability from k-mer mappability
    Pmap_Arrays = convert_kmap_to_pmap_arrays(Kmap_Arrays, kmer_len)

    ## 5) Output the calculated pileup mappability values as a .bedgraph table

    Pmap_BEDGRAPH_DF = process_nparrays_to_bedgraph_df(Pmap_Arrays)

    Pmap_BEDGRAPH_DF.to_csv(output_PileupMap_BG,
                          sep = "\t", index = False)

    print(f" Calculated pileup mappability scores (.bedgraph) output to: {output_PileupMap_BG}")
    
    ## 6) Output all regions with Pileup Mappability below 1 to a .bed file

    o_PileupMap_Below1_BED = args.outdir + f"/{Genome_FA_BaseName}.PileupMap.K{kmer_length}_E{errors}.bed"

    Pmap_Below1Merge_DF = getRegions_BelowThreshold_PupMap(Pmap_BEDGRAPH_DF, threshold = 1)
    
    Pmap_Below1Merge_DF.to_csv(o_PileupMap_Below1_BED, 
                          sep = "\t", index = False, header=None)

    print(f" All regions w/ pileup mappability scores below 1 (.bed) output to: {o_PileupMap_Below1_BED}")
    

    # Step 7: OPTIONAL - Calculate average pileup mappability across all features in genome's GFF file
    if input_Genome_GFF != None:
        
        print(f" Starting calculation of mean pileup mappability per feature in provided GFF file ({input_Genome_GFF})")

        PupMap_PerFeature_DF = calc_pupmap_per_gene(Pmap_Arrays,
                                                    input_Genome_GFF)

        o_Features_Wi_MeanPupmap_TSV = args.outdir + f"/{Genome_FA_BaseName}.PerFeature.MeanPileupMap.K{kmer_length}_E{errors}.tsv"

        PupMap_PerFeature_DF.to_csv(o_Features_Wi_MeanPupmap_TSV, 
                                    sep = "\t", index = False)

    print(f"\n Annotated features scored by mean pileup mappability output to: \n   {o_Features_Wi_MeanPupmap_TSV}\n")

    # Step 8: Optional - output all calculated pileup mappability values as compressed numpy arrays (.npz)

    if save_npz:
        pmap_npz_dir = args.outdir + f"/{Genome_FA_BaseName}.PileupMap.NumpyArrays/"
        save_numpy_array_dict(Pmap_Arrays, pmap_npz_dir)


def _genmap_and_pup_cli(args):
    """
    Run Genmap indexing and mapping steps on input genome fasta file, then calculate pileup mappability from k-mer mappability values. 

    This function combines all of the steps of the 2 sub-pipelines: 'run_genmap' and 'pup' into a single function.
    """
    
    # Step 1: Set input parameters and define variables for relevant paths
    input_Genome_FA = args.in_genome_fa 
    kmer_length = args.kmer_len
    errors = args.errors

    input_Genome_GFF = args.gff 

    save_npz = args.save_numpy

    index_dir = args.outdir + "/genmap_index"
    kmap_dir = args.outdir + "/genmap_kmap_K" + str(kmer_length) + "_E" + str(errors)
    
    Genome_FA_BaseName = get_fasta_basename(input_Genome_FA)

    genmap_out_prefix = kmap_dir + "/" + Genome_FA_BaseName + ".kmermap." + "K" + str(kmer_length) + "_E" + str(errors)

    # Create the target output directory
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(kmap_dir, exist_ok=True)

    # Step 2: Check if Genmap is available in the path
    if not is_genmap_available():
        raise EnvironmentError("genmap is not available on the PATH. Please install it or fix your PATH.")
        
    # Step 3: Run Genmap Index step for input genome

    run_genmap_index(input_Genome_FA,
                     index_dir)

    # Step 4: Run Genmap MAP step for genome (after indexing) to calculate k-mer mappability

    run_genmap_map(index_dir,
                   genmap_out_prefix,
                   kmer_length,
                   errors)
    
    print(f"Genmap processing completed.\nOutput files saved to: {args.outdir}")
    print(f"\nIndexed genome files saved to: {index_dir}")
    print(f"\nK-mer mappability files saved to: {kmap_dir}")
    print("\n")

    # Step 5: Run pileup mappability calculation on k-mer mappability values

    ## 5.1) Define input and output file paths (.bedgraph format) ####

    i_KMap_BG = f'{genmap_out_prefix}.bedgraph'

    o_PileupMap_BG = args.outdir + f"/{Genome_FA_BaseName}.PileupMap.K{kmer_length}_E{errors}.bedgraph"

    kmer_len = args.kmer_len

    ## 5.2) Parse k-mer mappability bedgraph file as Pandas DF
    Kmap_DF = kmap_bedgraph_to_DF(i_KMap_BG)

    ## 3) Convert kmap_DF to dict of numpy arrays of k-mer mappability
    Kmap_Arrays = kmap_DF_To_ArrayDict(Kmap_DF)

    ## 4) For each chromosome, calculate Pileup Mappability from k-mer mappability
    Pmap_Arrays = convert_kmap_to_pmap_arrays(Kmap_Arrays, kmer_len)

    ## 5) Output the calculated pileup mappability values as a .bedgraph table

    Pmap_BEDGRAPH_DF = process_nparrays_to_bedgraph_df(Pmap_Arrays)

    Pmap_BEDGRAPH_DF.to_csv(o_PileupMap_BG,
                          sep = "\t", index = False, header=None)

    print(f"\n Calculated pileup mappability scores (.bedgraph) output to: {o_PileupMap_BG}")
    
    ## 6) Output all regions with Pileup Mappability below 1 to a .bed file

    o_PileupMap_Below1_BED = args.outdir + f"/{Genome_FA_BaseName}.LowPileupMap.Below1.K{kmer_length}_E{errors}.bed"

    Pmap_Below1Merge_DF = getRegions_BelowThreshold_PupMap(Pmap_BEDGRAPH_DF, threshold = 1)
    
    Pmap_Below1Merge_DF.to_csv(o_PileupMap_Below1_BED, 
                          sep = "\t", index = False, header=None)

    print(f"\n All regions w/ pileup mappability scores below 1 (.bed) output to: {o_PileupMap_Below1_BED}")
    

    # N_Pos_PmapBelowThreshold = Pmap_BEDGRAPH_DF.query(f"score < 1 & score >= 0").shape[0]
    # print(f"\n\n    Number of REGIONs with pileup mappability scores below 1: {N_Pos_PmapBelowThreshold} ")


    # Step 7: OPTIONAL - Calculate average pileup mappability across all features in genome's GFF file
    if input_Genome_GFF != None:
        
        print(f"\n Starting calculation of mean pileup mappability per feature in user provided GFF file:  \n   {input_Genome_GFF}")

        PupMap_PerFeature_DF = calc_pupmap_per_gene(Pmap_Arrays,
                                                    input_Genome_GFF)

        o_Features_Wi_MeanPupmap_TSV = args.outdir + f"/{Genome_FA_BaseName}.PerFeature.MeanPileupMap.K{kmer_length}_E{errors}.tsv"

        PupMap_PerFeature_DF.to_csv(o_Features_Wi_MeanPupmap_TSV, 
                                    sep = "\t", index = False)

    print(f"\n Annotated features scored by mean pileup mappability output to: \n   {o_Features_Wi_MeanPupmap_TSV}\n")


    # Step 8: Optional - output all calculated pileup mappability values as compressed numpy arrays (.npz)

    if save_npz:
        pmap_npz_dir = args.outdir + f"/{Genome_FA_BaseName}.PileupMap.NumpyArrays/"
        save_numpy_array_dict(Pmap_Arrays, pmap_npz_dir)



def main():

    ascii_art = r"""
  _____             __  __                             
 |  __ \           |  \/  |                            
 | |__) |   _ _ __ | \  / | __ _ _ __  _ __   ___ _ __ 
 |  ___/ | | | '_ \| |\/| |/ _` | '_ \| '_ \ / _ \ '__|
 | |   | |_| | |_) | |  | | (_| | |_) | |_) |  __/ |   
 |_|    \__,_| .__/|_|  |_|\__,_| .__/| .__/ \___|_|   
             | |                | |   | |              
             |_|                |_|   |_|              
    """

    # Dynamically get terminal width
    terminal_width = shutil.get_terminal_size().columns
    
    pupmapper_general_description = """
Toolkit for calculating Pileup and k-mer mappability values from a genome sequence. 

The first step is a wrapper for running the Genmap software's indexing and k-mer mappability steps.
The second step calculates pileup mappability from k-mer mappability values generated by Genmap.
                                    """

    parser = argparse.ArgumentParser(description = f"{ascii_art} \n {pupmapper_general_description}",
                                     formatter_class = lambda prog: argparse.RawTextHelpFormatter(prog, width=terminal_width) )

    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')



    sub_parser_1 = parser.add_subparsers(required=True, help='Please select one of the sub-pipelines of the Pupmapper.\n')

    # 1) Create parser for running ALL steps (Genmap then Pileup mappability calculation)
    run_genmap_and_pup_parser = sub_parser_1.add_parser("run_all", help="Run Genmap indexing and mapping steps, then calculate pileup mappability from k-mer mappability values.")

    run_genmap_and_pup_parser.add_argument('-i', '--in_genome_fa', type=str, required=True,
                                  help="Input genome fasta file (.fasta)")
    
    run_genmap_and_pup_parser.add_argument('-o', '--outdir', type=str, required=True,
                                  help="Directory for all outputs of k-mer and pileup mappability processing.")

    run_genmap_and_pup_parser.add_argument('-k', '--kmer_len',type=int, required=True,
                                  help="k-mer length (bp) used to generate the k-mer mappability values")
    
    run_genmap_and_pup_parser.add_argument('-e', '--errors',type=int, required=True,
                                  help="Number of errors (mismatches) allowed in Genmap's k-mer mappability calculation")
    
    run_genmap_and_pup_parser.add_argument('-g', '--gff', type=str, required=False,
                                  help="GFF formatted genome annotations for input genome (.gff) (Optional)")

    run_genmap_and_pup_parser.add_argument('--save-numpy', action='store_true', required=False,
                                  help="If enabled, all pileup mappability scores will be output as compressed numpy arrays (.npz).")

    run_genmap_and_pup_parser.set_defaults(func=_genmap_and_pup_cli)


    # 2) Create parser for running Genmap Index + Map steps
    run_genmap_parser = sub_parser_1.add_parser("run_genmap", help="Wrapper command for running Genmap processing on input FASTA (2 Steps, indexing + k-mer mappability calculation. For more control of parameters, run 'genmap index' and 'genmap map' separately)")

    run_genmap_parser.add_argument('-i', '--in_genome_fa', type=str, required=True,
                                  help="Input genome fasta file (.fasta)")
    
    run_genmap_parser.add_argument('-o', '--outdir', type=str, required=True,
                                  help="Directory for all outputs of k-mer mappability processing w/ Genmap (indexing and mapping steps) ")

    run_genmap_parser.add_argument('-k', '--kmer_len',type=int, required=True,
                                  help="k-mer length (bp) used to generate the k-mer mappability values")
    
    run_genmap_parser.add_argument('-e', '--errors',type=int, required=True,
                                  help="Number of errors (mismatches) allowed in Genmap's k-mer mappability calculation")
    
    run_genmap_parser.set_defaults(func=_genmap_2steps_cli)


    # 3) Create k-map to pup parser
    kmap_to_pup_parser = sub_parser_1.add_parser("run_pileup", help="Calculate pileup mappability based on k-mer mappability values (bedgraph format)")

    kmap_to_pup_parser.add_argument('-i', '--input', type=str, required=True,
                                  help="Input k-mer mappability values in bedgraph format (.bedgraph). Ideally, generated with genmap software")

    kmap_to_pup_parser.add_argument('-o', '--output', type=str, required=True,
                                                    help="Output table of pileup mappability (.bedgraph)")

    kmap_to_pup_parser.add_argument('-k', '--kmer_len',type=int, required=True,
                                  help="k-mer length (bp) used to generate the input k-mer mappability values")

    kmap_to_pup_parser.add_argument('-g', '--gff', type=str, required=False,
                                  help="GFF formatted genome annotations for input genome (.gff) (Optional Input)")

    kmap_to_pup_parser.add_argument('--save-numpy', action='store_true', required=False,
                                  help="If enabled, all pileup mappability scores will be output as compressed numpy arrays (.npz).")

    kmap_to_pup_parser.set_defaults(func=_pup_cli)


    ### Check if no arguments were provided, if so exit with help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    # 4) Handle arguments passed by user
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())

