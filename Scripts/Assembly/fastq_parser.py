#!/usr/bin/env python3

"""
FastQ File Reader

This module reads a FastQ data set and returns a dictionary with the sequence titles as keys
and the sequences as values.

Author: Lisan Eisinga
Date: 26-07-2023

Input data set formats:
- fastq_file: FastQ data set containing biological sequence and quality score information (MinION).
"""

import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def read_fastq_file(fastq_file: str) -> dict[str, str]:
    """
    Reads a FastQ data set and returns a dictionary with the sequence titles as keys
    and the sequences as values.

    Parameters:
        fastq_file (str): The path to the FastQ data set.

    Returns:
        dict: A dictionary with the sequence titles as keys and the sequences as values.

    Raises:
        FileNotFoundError: If the FastQ data set does not exist.
        ValueError: If the FastQ data set is empty or contains invalid records.
    """
    # Check if the FastQ file exists
    if not os.path.exists(fastq_file):
        raise FileNotFoundError(f"FastQ data set '{fastq_file}' not found.")

    # Create an empty dictionary to store the sequences
    sequences = {}

    # Open the FastQ file and read its contents
    with open(fastq_file, "r", encoding="utf-8") as fastq:
        try:
            # Iterate over each record in the FastQ file
            for title, seq, _ in FastqGeneralIterator(fastq):
                # Use the first element of the title as the read ID to avoid possible duplicates
                read_id = title.split()[0].strip()
                if read_id:
                    # Add the sequence to the dictionary with the read ID as the key
                    sequences[read_id] = seq
        except ValueError as v_e:
            # If the FastQ file contains invalid records, raise an error
            raise ValueError("Invalid FastQ format. Unable to read the FastQ data set.") from v_e

    # If no sequences were added to the dictionary, raise an error
    if not sequences:
        raise ValueError("The FastQ data set is empty or does not contain valid records.")

    # Return the dictionary of sequences
    return sequences
