#!/usr/bin/env python3
import os
import gzip
import numpy as np
import pandas as pd
import bioframe as bf
from tqdm import tqdm


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

    # The window size is k_size, so the first elements will be the mean of the first k_size elements
    pmap_first_positions = np.array([np.mean(kmap_array[: f]) for f in np.arange(1, k_size)])

    # The last k positions will have the same pileup mappability value as the last kmer mappability value
    last_k_mappability_value = kmap_array[-1]
    pmap_last_positions = np.full((k_size - 1,), last_k_mappability_value) # Length is k -1

    pmap_Array_Final = np.concatenate( (pmap_first_positions,
                                         pmap_Array,
                                         pmap_last_positions) )

    return pmap_Array_Final


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




def summarize_pileup_map(i_Pmap_DF, genome_name, k, e):
    
    # Calculate length sums for different score thresholds
    i_Pmap_DF["length"] = i_Pmap_DF["end"] - i_Pmap_DF["start"]
    total_length = i_Pmap_DF["length"].sum()
    pmap_bp_below1 = i_Pmap_DF.query("score >= 0 & score < 1")["length"].sum()
    pmap_bp_below09 = i_Pmap_DF.query("score >= 0 & score < 0.9")["length"].sum()
    pmap_bp_below075 = i_Pmap_DF.query("score >= 0 & score < 0.75")["length"].sum()
    pmap_bp_below05 = i_Pmap_DF.query("score >= 0 & score < 0.5")["length"].sum()
    pmap_bp_below025 = i_Pmap_DF.query("score >= 0 & score < 0.25")["length"].sum()

    # Create a single-row DataFrame with the calculated values
    pmap_summ_df = pd.DataFrame({
        "Genome": [genome_name],
        "K": [k],
        "E": [e],
        "Length": [total_length],
        "Pupmap_Below1": [pmap_bp_below1],
        "Pupmap_Below0.9": [pmap_bp_below09],
        "Pupmap_Below0.7": [pmap_bp_below075],
        "Pupmap_Below0.5": [pmap_bp_below05],
        "Pupmap_Below0.25": [pmap_bp_below025]
    })
    
    return pmap_summ_df



def summarize_pileup_map_per_chromosome(i_Pmap_DF, genome_name, k, e):
    """
    Summarizes the pileup mappability per chromosome.
    """
    # Calculate length sums for different score thresholds per chromosome
    i_Pmap_DF["length"] = i_Pmap_DF["end"] - i_Pmap_DF["start"]
    summary_list = []

    for chrom, group in i_Pmap_DF.groupby("chrom"):
        total_length = group["length"].sum()
        pmap_bp_below1 = group.query("score >= 0 & score < 1")["length"].sum()
        pmap_bp_below09 = group.query("score >= 0 & score < 0.9")["length"].sum()
        pmap_bp_below075 = group.query("score >= 0 & score < 0.75")["length"].sum()
        pmap_bp_below05 = group.query("score >= 0 & score < 0.5")["length"].sum()
        pmap_bp_below025 = group.query("score >= 0 & score < 0.25")["length"].sum()

        summary_list.append({
            "Genome": genome_name,
            "Chromosome": chrom,
            "K": k,
            "E": e,
            "Length": total_length,
            "Pupmap_Below1": pmap_bp_below1,
            "Pupmap_Below0.9": pmap_bp_below09,
            "Pupmap_Below0.7": pmap_bp_below075,
            "Pupmap_Below0.5": pmap_bp_below05,
            "Pupmap_Below0.25": pmap_bp_below025
        })

    summary_df = pd.DataFrame(summary_list)
    return summary_df





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


def calc_pupmap_per_annotated_feature(Pmap_Arrays, input_GFF_PATH):
    """

    """

    # Step 1: Parse GFF as a Pandas Dataframe 

    Feat_DF = parse_gff3_with_pandas(input_GFF_PATH)
    Feat_DF.rename(columns={'seqname': 'chrom'}, inplace=True)


    # Step 2: Filter dataframe of GFF for features on analyzed chromosomes (SequenceIDs)

    Feat_DF = Feat_DF[Feat_DF["chrom"].isin(Pmap_Arrays.keys())]

    # Step 3: Calculate mean Pileup mappability score for each feature

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



def PileupMappability_bedgraph_to_DF(in_bedgraph_PATH):
    """ """
    i_PupMap_DF = pd.read_csv(in_bedgraph_PATH, sep = "\t", header=None)
    i_PupMap_DF.columns = ("chrom", "start", "end", "PupMap")

    return i_PupMap_DF

def parse_3column_BED_to_DF(in_bed_PATH):
    """ """
    BED_DF = pd.read_csv(in_bed_PATH, sep = "\t", header=None)
    BED_DF.columns = ("chrom", "start", "end")

    return BED_DF


def convert_gz_to_txt(gzipped_file, output_file):
    """
    Converts a gzipped file to a plain text file
    
    Parameters:
        gzipped_file (str): Path to the input gzipped file.
        output_file (str): Path to the output plain text file.
    """
    with gzip.open(gzipped_file, 'rt') as gz_file, open(output_file, 'w') as txt_file:
        for line in gz_file:
            txt_file.write(line)