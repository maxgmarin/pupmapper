#!/usr/bin/env python3

import pandas as pd
import numpy as np
from tqdm import tqdm


def kmap_bedgraph_to_DF(in_bedgraph_PATH):
    """ """
    i_Kmap_DF = pd.read_csv(in_bedgraph_PATH, sep = "\t", header=None)
    i_Kmap_DF.columns = ("chrom", "start", "end", "kmap")

    return i_Kmap_DF


def kmap_DF_To_ArrayDict(i_Kmap_DF):
    """ """
    # Get the maximum end values for each chromosome
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
    
    BEDGRAPH_DF.columns = ["chrom", "chromStart", "chromEnd", "score" ]
    
    return BEDGRAPH_DF



def process_nparrays_to_bedgraph_df(nparrays_dict):
    all_dfs = []
    for seq_id, input_PmapArray in nparrays_dict.items():
        BED_DF = convert_GenomeNParray_To_BEDGRAPH_DF(input_PmapArray, seq_id)
        all_dfs.append(BED_DF)
    combined_df = pd.concat(all_dfs, ignore_index=True)
    return combined_df







