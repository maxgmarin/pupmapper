#!/usr/bin/env python3

# This script contains functions relevant to running bigtools (as a wrapper)

import os 
import subprocess
import shutil


def is_bigtools_available():
    return shutil.which("bigtools") is not None


def run_bigtools_bedgraphtobigwig(bedgraph_file, chromsizes_file, output_file, threads=1):
    # Construct the command
    cmd = ["bigtools", "bedgraphtobigwig", bedgraph_file, chromsizes_file, output_file,
           "-t", str(threads)]

    # Print the command for visibility
    print(f"Running 'bigtools bedgraphtobigwig' with command: \n{' '.join(cmd)} \n")

    # Execute the command
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check for errors
    if result.returncode != 0:
        print(f"Error running bigtools bedgraphtobigwig: {result.stderr}")
    else:
        print("bigtools bedgraphtobigwig completed successfully.\n")

