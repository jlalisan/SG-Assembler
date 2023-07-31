#!/usr/bin/env python3

"""
Directed Graph-Based Genome Assembler

This script reads FASTQ and PAF files, constructs a directed graph based on read sequences and alignment information,
and performs genome assembly to generate contigs (Small scale simplified example).

Usage:
    Modify the 'fastq_file_path' and 'paf_file_path' variables in the main function to specify the paths to your
    FASTQ and PAF files. Then, run the script to perform the genome assembly and obtain contigs.

Author: Lisan Eisinga
"""

import os
import networkx as nx
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqRecord import SeqRecord


def read_fastq_file(fastq_file):
    """
    Read a FASTQ file and return the read sequences as a dictionary with read IDs as keys.

    Parameters:
        fastq_file (str): Path to the FASTQ file.

    Returns:
        dict: A dictionary with read IDs as keys and read sequences as values.
    """
    if not os.path.exists(fastq_file):
        raise FileNotFoundError(f"FastQ data set {fastq_file} not found.")

    sequences = { }
    with open(fastq_file, "r", encoding="utf-8") as fastq:
        for title, seq, _ in FastqGeneralIterator(fastq):
            # Use the first element to avoid possible duplicated read_ids
            sequences[title.split()[0]] = seq

    return sequences


def parse_paf(paf_file):
    """
    Parse a PAF file and return the alignment information as a DataFrame.

    Parameters:
        paf_file (str): Path to the PAF file.

    Returns:
        pandas.DataFrame: DataFrame containing alignment information.
    """
    if not os.path.exists(paf_file):
        raise FileNotFoundError(f"PAF data set '{paf_file}' not found.")

    # Read the PAF data set into a DataFrame
    columns = [
        "query_id", "query_length", "query_start", "query_end", "strand",
        "target_id", "target_length", "target_start", "target_end",
        "alignment_block_length", "residue_matches", "mapping_quality"
    ]
    paf_df = pd.read_csv(paf_file, sep="\t", header=None, usecols=range(12), names=columns)
    return paf_df


def construct_directed_graph(sequences, paf_df):
    """
    Construct a directed graph based on read sequences and alignment information.

    Parameters:
        sequences (dict): Dictionary of read sequences with read IDs as keys.
        paf_df (pandas.DataFrame): DataFrame containing alignment information.

    Returns:
        networkx.DiGraph: Directed graph with reads as nodes and alignments as edges.
    """
    # Create an empty directed graph
    graph = nx.DiGraph()

    # Add nodes to the graph (each read sequence will be a node)
    for read_id, sequence in sequences.items():
        graph.add_node(read_id, sequence=sequence)

    # Add directed edges to the graph based on the PAF alignment data
    for _, row in paf_df.iterrows():
        query_id = row["query_id"]
        target_id = row["target_id"]

        # Add an edge from the query read to the target read
        graph.add_edge(query_id, target_id)

    return graph


def genome_assembly(directed_graph):
    """
    Perform genome assembly using the directed graph.

    Parameters:
        directed_graph: Directed graph for traversal.
    Returns:
        list: List of contigs (assembled sequences).
    """
    contigs = []

    def dfs(node, current_contig):
        if directed_graph.out_degree(node) == 0:
            # Reached the end of the contig, add it to the list of contigs
            contigs.append(current_contig + directed_graph.nodes[node]["sequence"])
        else:
            # Continue traversing the graph
            current_contig += directed_graph.nodes[node]["sequence"]
            for neighbor in directed_graph.successors(node):
                dfs(neighbor, current_contig)

    for start_node in directed_graph.nodes():
        if directed_graph.in_degree(start_node) > 0 and directed_graph.out_degree(start_node) == 0:
            # Start a new contig from this start node
            dfs(start_node, "")

    return contigs


def write_to_file(contigs, output_file):
    """
    Write contigs to a FASTA file.

    Parameters:
        contigs (list): List of contigs (assembled sequences).
        output_file (str): Path to the output FASTA file.
    """
    with open(output_file, "w") as fasta_file:
        for i, contig_sequence in enumerate(contigs, start=1):
            seq_record = SeqRecord(Seq(contig_sequence), id=f"Contig_{i}", description="")
            SeqIO.write(seq_record, fasta_file, "fasta")


def main():
    # Example usage:
    fastq_file_path = "../../Test_data/foo-reads.fq"
    paf_file_path = "../../Test_data/foo.paf"
    output_file_path = "../../Test_data/foo_directed_contigs.fq"

    sequences = read_fastq_file(fastq_file_path)
    paf_df = parse_paf(paf_file_path)
    directed_graph = construct_directed_graph(sequences, paf_df)
    assembled_contigs = genome_assembly(directed_graph)

    # Write contigs to a FASTA file
    write_to_file(assembled_contigs, output_file_path)


if __name__ == "__main__":
    main()
