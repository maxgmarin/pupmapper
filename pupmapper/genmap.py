#!/usr/bin/env python3

# This script contains functions relevant to running Genmap (as a wrapper)

import os 
import subprocess
import shutil


def infer_genmap_prefix_from_fasta_path(fasta_path):
    # Get the base name of the file (i.e., with extension)
    base_name = os.path.basename(fasta_path)
    
    # Split the base name into the name and extension
    file_name_without_ext, _ = os.path.splitext(base_name)
    
    genmap_prefix = file_name_without_ext + ".genmap.kmermap"

    return genmap_prefix

def get_fasta_basename(fasta_path):
    # Get the base name of the file (i.e., with extension)
    base_name = os.path.basename(fasta_path)
    
    # Split the base name into the name and extension
    file_name_without_ext, _ = os.path.splitext(base_name)
    
    return file_name_without_ext


def is_genmap_available():
    return shutil.which("genmap") is not None


def run_genmap_index(fasta_path, index_dir):
    cmd = ["genmap", "index", "-F", fasta_path, "-I", index_dir]

    print(f"Running 'genmap index' with command: \n{' '.join(cmd)} \n")

    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running genmap index: {result.stderr}")
        
    else:
        print("genmap index completed successfully.\n")


def run_genmap_map(index_dir, kmap_out_prefix, k_length, errors, threads=1):

    cmd = ["genmap", "map", "-K", str(k_length),
                            "-E", str(errors),
                            "-I", index_dir,
                            "-O", kmap_out_prefix,
                            "--threads", str(threads),
                            "-w",
                            "-bg"]

    print(f"Running 'genmap map' with command: \n{' '.join(cmd)} \n")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running genmap map: {result.stderr}")
    else:
        print("genmap map completed successfully.\n")

