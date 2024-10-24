#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import bioframe as bf
#import pybigtools as pbt
#from pybigtools import BigWigWrite
import pybigtools 

def kmap_bedgraph_to_DF(in_bedgraph_PATH):
    """ """
    i_Kmap_DF = pd.read_csv(in_bedgraph_PATH, sep = "\t", header=None)
    i_Kmap_DF.columns = ("chrom", "start", "end", "kmap")

    return i_Kmap_DF


def kmap_DF_To_ArrayDict(i_Kmap_DF):
    """ """
    # Get the maximum end values for each chromosome
    # (Note this end value is the length of the chromosome - K).
    # This is how Genmap outputs the kmer mappability info

    Dict_ChrToMaxEnd = i_Kmap_DF.groupby("chrom")["end"].max().to_dict()

    # Initialize a dictionary to hold the numpy array for each chromosome
    Kmap_Arrays = {}
    
    # Create a numpy array for each chromosome based on the maximum end value
    for chrom, max_end in Dict_ChrToMaxEnd.items():
        Kmap_Arrays[chrom] = np.zeros(max_end, dtype=float)  # Assuming the 'end' is 1-indexed and inclusive
    
    # Update the numpy array for each chromosome based on the bedgraph values
    for _, row in i_Kmap_DF.iterrows():
        chrom, start, end, value = row
        Kmap_Arrays[chrom][start:end] = value  # Assuming 'start' is 0-indexed and 'end' is exclusive

    return Kmap_Arrays


def kmap_bedgraph_to_DFandArray(kmap_bedgraph_PATH):
    """ """
    Kmap_DF = kmap_bedgraph_to_DF(kmap_bedgraph_PATH)
    Kmap_Arrays = kmap_DF_To_ArrayDict(Kmap_DF)

    return Kmap_DF, Kmap_Arrays


def compute_pmap(kmap_array, k_size):
    """ """
    # This list comprehension will iterate over all windows of the array and calculate the mean kmer mappability within that window
    pmap_Array = np.array([ np.mean(kmap_array[e - (k_size): e ]) for e in np.arange(k_size, len(kmap_array) + 1 ) ])

    # Create empty array of NaN values for padding pmap array
    #Pad_NaN = np.full((k_size - 1,), np.nan) # Length is k -1
    Pad_NaN = np.full((k_size - 1,), -1) # Length is k -1

    pmap_Array_Padded = np.concatenate( (Pad_NaN,
                                         pmap_Array,
                                         Pad_NaN) )

    return pmap_Array_Padded


def convert_kmap_to_pmap_arrays(Kmap_Arrays, k_size):

    Pmap_Arrays = {}

    for seq_id, kmap_array in Kmap_Arrays.items():

        Pmap_Arrays[seq_id] = compute_pmap(kmap_array, k_size)

    return Pmap_Arrays



# BED format specifications: https://useast.ensembl.org/info/website/upload/bed.html

def convert_GenomeNParray_To_BEDGRAPH_DF(input_PmapArray, seq_id):
    """ """
    last_Score = input_PmapArray[0]

    startOfRegion = 0
    listOfBED_Tuples = []

    for RefPos_0based in tqdm(np.arange(1, len(input_PmapArray))):

        Score = input_PmapArray[RefPos_0based]

        if Score != last_Score:

            endOfRegion = RefPos_0based

            BED_EntryTuple = (seq_id, startOfRegion, endOfRegion, last_Score,)
            
            listOfBED_Tuples.append(BED_EntryTuple)

            startOfRegion = RefPos_0based 

            #1 Output the last range
            #2 Store the new score    

        last_Score = Score #2 Store the new score   

    endOfRegion = RefPos_0based + 1

    BED_EntryTuple = (seq_id, startOfRegion, endOfRegion, last_Score)
    listOfBED_Tuples.append(BED_EntryTuple)       

    BEDGRAPH_DF = pd.DataFrame(listOfBED_Tuples)
    
    BEDGRAPH_DF.columns = ["chrom", "start", "end", "score" ]
    
    return BEDGRAPH_DF



def process_nparrays_to_bedgraph_df(nparrays_dict):
    all_dfs = []

    for seq_id, input_PmapArray in nparrays_dict.items():
        BED_DF = convert_GenomeNParray_To_BEDGRAPH_DF(input_PmapArray, seq_id)
        all_dfs.append(BED_DF)

    Pmap_DF = pd.concat(all_dfs, ignore_index=True)
    return Pmap_DF



def validate_PupMap_DF(i_Pmap_DF):
    # Validate that the input is a DataFrame with required columns
    if not isinstance(i_Pmap_DF, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame")
    
    required_columns = ["chrom", "start", "end", "score"]
    if list(i_Pmap_DF.columns) != required_columns:
        raise ValueError(f"Pileup Mappability DataFrame must have columns: {required_columns}")


def infer_ChrLengths(i_Pmap_DF):
    """
    Computes the maximum 'end' position for each chromosome in the DataFrame,
    returning a dictionary where keys are chromosome names and values are their lengths.
    """
    Dict_ChrToMaxEnd = i_Pmap_DF.groupby("chrom")["end"].max().to_dict()
    return Dict_ChrToMaxEnd



def pmap_df_to_bigwig(i_Pmap_DF, out_bigwig_path):
    """
    Converts a DataFrame to a BigWig file using the BigWigWrite class.

    Parameters:
    - i_Pmap_DF: Pandas DataFrame with columns ["chrom", "start", "end", "score"]
    - out_bigwig_path: Path where the output BigWig file will be written
    """
    
    # Step 1: Validate the input DataFrame
    validate_PupMap_DF(i_Pmap_DF)

    # Step 2: Infer chromosome lengths from the DataFrame
    chrom_lengths_dict = infer_ChrLengths(i_Pmap_DF)
    print(chrom_lengths_dict)

    i_Pmap_DF = i_Pmap_DF.query("score >= 1")

    # Step 3: Create the iterable with values (chromosome, start, end, score)

    values_list = ((row["chrom"], row["start"], row["end"], row["score"]) for _, row in i_Pmap_DF.head().iterrows())
    #values_list = [(row["chrom"], row["start"], row["end"], row["score"]) for _, row in i_Pmap_DF.head(1).iterrows()]

    # Step 4: Open the BigWig file and write the values
    out_BigWig = pybigtools.open(out_bigwig_path, "w")

    out_BigWig.write( chrom_lengths_dict, values_list)



def test_bigwig_write(i_Pmap_DF, out_bigwig_path):

    chroms = ["J02459.1"]
    clengths = {"J02459.1": 48502}

    def genintervals():
        import random

        for chrom in chroms:
            clength = clengths[chrom]
            current = random.randint(0, 300)
            start = current
            while True:
                length = random.randint(1, 200)
                end = start + length
                if end > clength:
                    break
                value = round(random.random(), 5)
                yield (chrom, start, end, value)
                start = end + random.randint(20, 50)

    intervals = list(genintervals())
    for i in intervals:
        print(i)
    print()
    b = pybigtools.open(out_bigwig_path, "w")
    b.write(clengths, iter(intervals))







def save_numpy_array_dict(array_dict, output_dir):
    """
    Save numpy arrays of pileup mappability from a dictionary of arrays to a specified directory.
    Each array is saved in a compressed format (.npz) with the key as the filename.

    Parameters:
    - array_dict: Dictionary where keys are filenames and values are numpy arrays.
    - output_dir: Directory where the arrays will be saved.
    """
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Iterate over the dictionary and save each numpy array
    for key, array in array_dict.items():
        output_file = os.path.join(output_dir, f"{key}.npz")
        np.savez_compressed(output_file, array)



#### Functions to merge and define regions w/ mappability below a threshold ####

def getRegions_BelowThreshold_PupMap(Pmap_BEDGRAPH_DF, threshold = 1):
    """
    This function takes a DF of pileup mappability and 
    merges all positions which are below a defined threshold (Default: 1)
    """

    # Step 1: Filter Pileup Mappability DF for positions BELOW the threshold
    Pmap_BelowThreshold_DF = Pmap_BEDGRAPH_DF.query(f"score < {threshold} & score >= 0")

    # Step 2: Use Bioframe to merge all neighboring ranges w/ a "low" pileup mappability
    if Pmap_BelowThreshold_DF.shape[0] > 0:
        Pmap_BelowAndMerged_DF = bf.merge(Pmap_BelowThreshold_DF, min_dist=0)

    else:
        Pmap_BelowAndMerged_DF = pd.DataFrame(columns=["chrom", "start", "end"])

    return Pmap_BelowAndMerged_DF[["chrom", "start", "end"]]

#####################################################################################################





#### Functions to calculate average pileup mappability across regions of interest (ie genes) ####


def parse_gff3_with_pandas(file_path):
    # Define column names based on the GFF3 format
    column_names = [
        "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"
    ]
    
    # Use pandas to read the file, skip lines starting with '#' (comments), and use tab as the delimiter
    df = pd.read_csv(
        file_path, 
        sep='\t', 
        comment='#',  # Ignores lines that start with '#'
        names=column_names,
        header=None  # No header in GFF files
    )

    return df

# Example usage:
# df = parse_gff3_with_pandas('example.gff3')
# print(df.head())



def calc_pupmap_per_gene(Pmap_Arrays, input_GFF_PATH):
    """

    """

    # Step 1: Parse GFF as a Pandas Dataframe 

    Feat_DF = parse_gff3_with_pandas(input_GFF_PATH)
    Feat_DF.rename(columns={'seqname': 'chrom'}, inplace=True)

    # Step 2: 

    Region_PmapScores = []

    # Iterate over DataFrame rows
    for index, row in Feat_DF.iterrows():
        i_chrom = row['chrom']
        i_start_0idx = row['start'] - 1
        i_end = row['end']


        # Subset the numpy array for the given genomic region
        subset_Pmap_Scores = Pmap_Arrays[i_chrom][i_start_0idx:i_end]

        # Calculate the average of the subset values
        mean_Pmap = np.mean(subset_Pmap_Scores)

        Region_PmapScores.append(mean_Pmap)
    
    Feat_DF["Mean_PupMap"] = Region_PmapScores

    return Feat_DF



