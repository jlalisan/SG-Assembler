#!/usr/bin/env python3

"""
Overlap-Layout Consensus-based Assembler

This script reads a set of DNA sequencing reads from a FASTQ file and their overlap information from a PAF file
to perform overlap-layout consensus assembly. It generates layouts based on overlapping reads and calculates the
consensus sequence for each layout. The resulting assembled consensus sequences are saved in a FASTA file.

Usage:
    python olc_assembler.py

Input:
    - A FASTQ file containing sequencing reads.
    - A PAF file containing overlap information between reads.

Output:
    - A FASTA file (olc_contigs.fasta) containing the assembled consensus sequences.

Author: Lisan Eisinga
"""

import os
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def read_fastq_file(fastq_file):
    """
    Read a FASTQ file and extract sequences.

    Parameters:
        fastq_file (str): Path to the FASTQ file.

    Returns:
        dict: A dictionary with read IDs as keys and sequences as values.
    """
    if not os.path.exists(fastq_file):
        raise FileNotFoundError(f"FastQ data set {fastq_file} not found.")

    sequences = {}
    with open(fastq_file, "r", encoding="utf-8") as fastq:
        for title, seq, _ in FastqGeneralIterator(fastq):
            # Use the first element to avoid possible duplicated read IDs
            sequences[title.split()[0]] = seq

    return sequences


def parse_paf(paf_file):
    """
    Parse a PAF file and extract overlap information.

    Parameters:
        paf_file (str): Path to the PAF file.

    Returns:
        pandas.DataFrame: A DataFrame with overlap information columns.
    """
    if not os.path.exists(paf_file):
        raise FileNotFoundError(f"PAF data set '{paf_file}' not found.")

    # Read the PAF data set into a DataFrame
    columns = ["query_id", "query_length", "query_start", "query_end", "strand",
               "target_id", "target_length", "target_start", "target_end",
               "alignment_block_length", "residue_matches", "mapping_quality"]
    paf_df = pd.read_csv(paf_file, sep="\t", header=None, usecols=range(12), names=columns)
    return paf_df


def calculate_overlap_length(row):
    """
    Calculate the overlap length between two reads based on their alignment positions.

    Parameters:
        row (pandas.Series): A row containing alignment information.

    Returns:
        int: The length of the overlap between the two reads.
    """
    query_start, query_end = row["query_start"], row["query_end"]
    target_start, target_end = row["target_start"], row["target_end"]
    overlap_start = max(query_start, target_start)
    overlap_end = min(query_end, target_end)
    overlap_length = max(0, overlap_end - overlap_start)
    return overlap_length


def extend_layout(layout_id, paf_df, sequences, layout_reads, visited_reads):
    """
    Recursively extend a layout with overlapping reads.

    Parameters:
        layout_id (str): The read ID of the current layout.
        paf_df (pandas.DataFrame): DataFrame containing overlap information.
        sequences (dict): Dictionary with read IDs as keys and sequences as values.
        layout_reads (list): List to store read IDs in the current layout.
        visited_reads (set): Set to keep track of visited reads in the layout extension.
    """
    overlapping_reads = paf_df[paf_df["query_id"] == layout_id]
    for _, row in overlapping_reads.iterrows():
        target_id = row["target_id"]
        overlap_length = calculate_overlap_length(row)
        if overlap_length > 0 and target_id not in visited_reads:
            visited_reads.add(target_id)
            layout_reads.append(target_id)
            extend_layout(target_id, paf_df, sequences, layout_reads, visited_reads)


def layout_generation(paf_df, sequences):
    """
    Generate layouts based on the overlap information.

    Parameters:
        paf_df (pandas.DataFrame): DataFrame containing overlap information.
        sequences (dict): Dictionary with read IDs as keys and sequences as values.

    Returns:
        list: List of lists representing layouts with read sequences.
    """
    layouts = []
    visited_reads = set()
    for layout_id, _ in paf_df.groupby("query_id"):
        if layout_id not in visited_reads:
            layout_reads = [layout_id]
            visited_reads.add(layout_id)
            extend_layout(layout_id, paf_df, sequences, layout_reads, visited_reads)
            layout_sequences = [sequences[read_id] for read_id in layout_reads]
            layouts.append(layout_sequences)
    return layouts


def generate_consensus(layout):
    """
    Generate a consensus sequence from a layout.

    Parameters:
        layout (list): List of lists representing a layout with read sequences.

    Returns:
        str: The consensus sequence generated from the layout.
    """
    layout_length = max(len(read) for read in layout)
    consensus_seq = ""
    for i in range(layout_length):
        base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        valid_reads = [read for read in layout if i < len(read)]
        for read in valid_reads:
            base = read[i]
            base_counts[base] += 1
        consensus_base = max(base_counts, key=base_counts.get)
        consensus_seq += consensus_base
    return consensus_seq


def write_fasta_file(file_path, sequences):
    """
    Write sequences to a FASTA file.

    Parameters:
        file_path (str): Path to the output FASTA file.
        sequences (list): List of sequences to write to the file.
    """
    with open(file_path, "w") as fasta_file:
        for i, seq in enumerate(sequences, 1):
            fasta_file.write(f">Consensus_{i}\n")
            fasta_file.write(seq + "\n")


def main():
    fastq_file = "../../Test_data/foo-reads.fq"
    paf_file = "../../Test_data/foo.paf"
    output_file = "../../Test_data/olc_contigs.fasta"

    # Step 1: Read the FASTQ file and extract the sequences
    sequences = read_fastq_file(fastq_file)

    # Step 2: Parse the PAF file and extract the overlap information
    paf_df = parse_paf(paf_file)

    # Step 3: Generate layouts based on the overlap information
    layouts = layout_generation(paf_df, sequences)

    # Step 4: Generate consensus sequences from the layouts
    consensus_sequences = [generate_consensus(layout) for layout in layouts]

    # Step 5: Output the assembled sequences to a FASTA file
    write_fasta_file(output_file, consensus_sequences)


if __name__ == "__main__":
    main()
