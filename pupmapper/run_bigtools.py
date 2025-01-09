#!/usr/bin/env python3

# This script contains functions relevant to running bigtools (as a wrapper)

import subprocess
import shutil


def is_bigtools_available():
    return shutil.which("bigtools") is not None

def check_bigtools_version():
    try:
        result = subprocess.run(['bigtools', '--version'], capture_output=True, text=True, check=True)
        print(f"bigtools version:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"Error checking bigtools version:\n{e.stderr}")


def run_bigtools_bedgraphtobigwig(bedgraph_file, chromsizes_file, output_file, threads=1, verbose=False):
    # Construct the command
    cmd = ["bigtools", "bedgraphtobigwig", bedgraph_file, chromsizes_file, output_file,
           "-t", str(threads)]

    if verbose:
        print(f"Running 'bigtools bedgraphtobigwig' with command: \n{' '.join(cmd)} \n")

    # Execute the command
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check for errors
    if result.returncode != 0:
        print(f"Error running bigtools bedgraphtobigwig: {result.stderr}")
    else:
        if verbose:
            print("bigtools bedgraphtobigwig completed successfully.\n")
