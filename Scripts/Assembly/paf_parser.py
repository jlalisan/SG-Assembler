#!/usr/bin/env python

"""
PAF File Generator

This module creates a PAF (Pairwise mApping Format) data set from a FastQ data set,
using the "minimap2" command-line tool.
If the PAF data set already exists, it will not be recreated.

Author: Lisan Eisinga
Date: 26-07-2023

Input data set formats:
- fastq_file: FastQ data set containing biological sequence and quality score information (MinION).
"""

import logging
import os
import subprocess
import pandas as pd


def create_paf(fastq_file: str) -> str:
    """
    Create a PAF data set from a FastQ data set using the "minimap2" command-line tool.
    If the PAF data set already exists, it will not be recreated.

    Parameters:
        fastq_file (str): The path to the FastQ data set.

    Returns:
        str: The path to the PAF data set.

    Raises:
        FileNotFoundError: If the FastQ data set does not exist.
        subprocess.CalledProcessError: If the "minimap2" command fails.
        Exception: If there is an error while creating the PAF data set.
    """
    try:
        # Check if the FastQ file exists
        if not os.path.exists(fastq_file):
            raise FileNotFoundError(f"FastQ data set '{fastq_file}' not found.")

        # Generate the PAF data set name based on the FastQ data set name
        paf_file = f"{os.path.splitext(fastq_file)[0]}_overlaps.paf"

        # Check if the PAF data set already exists
        if os.path.exists(paf_file):
            logging.info("PAF data set '%s' already exists.", paf_file)
            return paf_file

        # Find the minimap2 executable
        minimap_executable = subprocess.run(
            'find / -type f -name "minimap2" -executable -print -quit 2>/dev/null',
            shell=True, capture_output=True, text=True, check=True
        )
        minimap_location = minimap_executable.stdout.strip()

        # Create the PAF data set
        minimap_command = f"{minimap_location} -x ava-ont {fastq_file} {fastq_file} > {paf_file}"
        subprocess.run(minimap_command, shell=True, check=True)
        logging.info("PAF data set '%s' created.", paf_file)

        return paf_file

    except FileNotFoundError as error:
        logging.error("FileNotFoundError: %s", str(error))
        raise

    except subprocess.CalledProcessError as error:
        logging.error("CalledProcessError: %s", str(error))
        raise


def parse_paf(paf_file: str) -> pd.DataFrame:
    """
    Reads a PAF data set and returns a pandas DataFrame with the relevant columns.

    Parameters:
        paf_file (str): The path to the PAF data set.

    Returns:
        pandas.DataFrame: A DataFrame with columns "query_id", "query_length", "query_start",
        "query_end", "strand", "target_id", "target_length", "target_start", and "target_end".

    Raises:
        FileNotFoundError: If the PAF data set does not exist.
        Exception: If there is an error while parsing the PAF data set.
    """
    if not os.path.exists(paf_file):
        raise FileNotFoundError(f"PAF data set '{paf_file}' not found.")

    # Read the PAF data set into a DataFrame
    columns = ["query_id", "query_length", "query_start", "query_end", "strand",
               "target_id", "target_length", "target_start", "target_end",
               "alignment_block_length", "residue_matches", "mapping_quality"]
    paf_df = pd.read_csv(paf_file, sep="\t", header=None, usecols=range(12), names=columns)
    return paf_df
