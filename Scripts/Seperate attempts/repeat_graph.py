#!/usr/bin/env python3

"""
Repeat Graph-Based Genome Assembler.

This script implements a simplified repeat graph-based assembler to generate contigs
from a set of reads provided in a FASTQ file (Small scale simplified example).

The assembler works as follows:
1. Read the FASTQ file to obtain a dictionary of read sequences.
2. Find overlaps between reads to construct a directed graph (repeat graph).
3. Identify repeat regions in the repeat graph using strongly connected components.
4. Generate contigs by concatenating reads based on overlaps within each repeat region.
5. Write the contigs to a FASTA file.

Author: Lisan Eisinga
"""

import os
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqRecord import SeqRecord


def read_fastq_file(fastq_file):
    """
    Read a FASTQ file and store the sequences in a dictionary.

    Parameters:
        fastq_file (str): Path to the input FASTQ file.

    Returns:
        dict: A dictionary containing read IDs as keys and sequences as values.
    """
    if not os.path.exists(fastq_file):
        raise FileNotFoundError(f"FastQ data set {fastq_file} not found.")

    sequences = { }
    with open(fastq_file, "r", encoding="utf-8") as fastq:
        for title, seq, _ in FastqGeneralIterator(fastq):
            # Use the first element to avoid possible duplicated read_ids
            sequences[title.split()[0]] = seq

    return sequences


def find_overlaps(sequences, min_overlap=30):
    """
    Find overlaps between reads.

    Parameters:
        sequences (dict): A dictionary containing read IDs as keys and sequences as values.
        min_overlap (int, optional): Minimum required overlap length between reads. Defaults to 30.

    Returns:
        list: A list of overlaps as tuples (read_id1, read_id2, overlap_length).
    """
    overlaps = []

    # Iterate through all pairs of reads
    for read_id1, seq1 in sequences.items():
        for read_id2, seq2 in sequences.items():
            if read_id1 == read_id2:
                continue  # Skip comparing a read to itself

            # Find the minimum length to consider an overlap
            min_length = min(len(seq1), len(seq2)) - min_overlap + 1

            # Only compare reads if their combined length is greater than min_overlap
            if min_length <= 0:
                continue

            # Look for seed matches in both directions (forward and reverse complement)
            for i in range(min_overlap):
                seed = seq1[-min_overlap + i:]
                if seq2.startswith(seed):
                    overlap_length = len(seq1) - i
                    overlaps.append((read_id1, read_id2, overlap_length))
                    break

                seed = seq2[-min_overlap + i:]
                if seq1.startswith(seed):
                    overlap_length = len(seq2) - i
                    overlaps.append((read_id2, read_id1, overlap_length))
                    break

    return overlaps


def construct_repeat_graph(sequences, overlaps):
    """
    Construct a directed graph representing read overlaps.

    Parameters:
        sequences (dict): A dictionary containing read IDs as keys and sequences as values.
        overlaps (list): A list of overlaps as tuples (read_id1, read_id2, overlap_length).

    Returns:
        nx.DiGraph: A directed graph representing the read overlaps.
    """
    # Create an empty directed graph
    repeat_graph = nx.DiGraph()

    # Add nodes to the graph (representing reads)
    repeat_graph.add_nodes_from(sequences.keys())

    # Add edges to the graph (representing overlaps)
    for read_id1, read_id2, overlap_length in overlaps:
        repeat_graph.add_edge(read_id1, read_id2, weight=overlap_length)

    return repeat_graph


def find_repeat_regions(repeat_graph):
    """
    Find repeat regions in the repeat graph using strongly connected components (SCCs).

    Parameters:
        repeat_graph (nx.DiGraph): A directed graph representing the read overlaps.

    Returns:
        list: A list of repeat regions, where each region is a set of read IDs.
    """
    # Find all strongly connected components in the graph
    strongly_connected_components = list(nx.strongly_connected_components(repeat_graph.to_directed()))

    # Filter out components with size > 1 as they indicate repeat regions
    repeat_regions = [component for component in strongly_connected_components if len(component) > 1]

    return repeat_regions


def dfs(graph, start_node, visited, contig, sequences):
    """
    Perform depth-first search (DFS) to generate a contig from a repeat region.

    Parameters:
        graph (nx.DiGraph): A directed graph representing the read overlaps.
        start_node: The starting node for the DFS traversal.
        visited (dict): A dictionary to keep track of visited nodes during DFS.
        contig (list): A list to store the contig as a sequence of reads.
        sequences (dict): A dictionary containing read IDs as keys and sequences as values.

    Returns:
        None: The contig is updated directly in the `contig` list.
    """
    # Mark the current node as visited
    visited[start_node] = True

    # Add the current node's sequence to the contig
    contig.append(sequences[start_node])

    # Explore all neighbors of the current node
    for neighbor in graph.neighbors(start_node):
        if not visited[neighbor]:
            dfs(graph, neighbor, visited, contig, sequences)


def generate_contigs(repeat_graph, repeat_regions, sequences):
    """
    Generate contigs from the repeat regions.

    Parameters:
        repeat_graph (nx.DiGraph): A directed graph representing the read overlaps.
        repeat_regions (list): A list of repeat regions, where each region is a set of read IDs.
        sequences (dict): A dictionary containing read IDs as keys and sequences as values.

    Returns:
        list: A list of contigs, each represented as a sequence of reads.
    """
    contigs = []

    for region in repeat_regions:
        valid_reads = [read_id for read_id in region if read_id in sequences]

        if len(valid_reads) < 2:
            continue  # Skip regions with less than two valid reads

        # Sort valid_reads based on their lengths in descending order
        valid_reads = sorted(valid_reads, key=lambda read_id: len(sequences[read_id]), reverse=True)

        # Generate the contig by concatenating reads based on overlaps
        contig_seq = sequences[valid_reads[0]]
        for i in range(1, len(valid_reads)):
            read_id = valid_reads[i]
            try:
                overlap_length = repeat_graph[valid_reads[i - 1]][read_id]['weight']
                contig_seq += sequences[read_id][overlap_length:]
            except KeyError:
                continue

        contigs.append(contig_seq)

    return contigs


def write_to_file(contigs, output_file):
    """
    Write contigs to a FASTA file.

    Parameters:
        contigs (list): A list of contigs, each represented as a sequence of reads.
        output_file (str): Path to the output FASTA file.

    Returns:
        None: The contigs are written to the output file.
    """
    records = []

    for i, contig in enumerate(contigs):
        # Create a sequence record for the contig
        contig_seq = "".join(contig)
        seq_record = SeqRecord(Seq(contig_seq), id=f"Contig{i + 1}", description="")

        records.append(seq_record)

    # Write the contig sequences to the output FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    print(f"Contigs written to {output_file}")


def main():
    """
    Main function to run the repeat graph assembler.

    This function reads a FASTQ file, constructs a repeat graph, identifies repeat regions,
    generates contigs, and writes the contigs to a FASTA file.
    """

    fastq_file = "../../Test_data/foo-reads.fq"
    output_file = "../../Test_data/repeat_graph_contigs.fasta"
    sequences = read_fastq_file(fastq_file)
    overlaps = find_overlaps(sequences)
    graph = construct_repeat_graph(sequences, overlaps)
    repeats = find_repeat_regions(graph)
    contigs = generate_contigs(graph, repeats, sequences)
    write_to_file(contigs, output_file)


if __name__ == "__main__":
    main()
