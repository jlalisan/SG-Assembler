#!/usr/bin/env python3

"""
Unfinished Script: Genome Assembly from Overlaps

This script aims to perform genome assembly by merging overlapping sequences obtained
from a FastQ file and overlaps information in a PAF file. It utilizes networkx and pandas
for graph-based operations and provides functions to read FastQ files, parse PAF files,
create a sparse graph, generate contigs, and merge similar contigs based on overlaps.

Note: This script is still unfinished and may require further development and testing.

"""
import time

import networkx as nx
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator



def read_fastq_file(fastq_file):
    """
    Read sequences from a FASTQ file and store them in a dictionary with sequence IDs as keys.

    Args:
        fastq_file (str): Path to the input FASTQ file.

    Returns:
        dict: A dictionary containing sequences with sequence IDs as keys.
    """
    sequences = {}
    with open(fastq_file) as f:
        for title, seq, qual in FastqGeneralIterator(f):
            sequences[title.split()[2]] = seq
    return sequences


def parse_paf(file_path):
    """
    Parse overlaps information from a PAF file and convert it into a DataFrame.

    Args:
        file_path (str): Path to the input PAF file.

    Returns:
        pandas.DataFrame: DataFrame containing the parsed overlaps information.
    """
    df = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], header=None)
    df.columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand', "target_id", 'target_length',
                  'target_start',
                  'target_end']
    return df


def sparsify_graph(sequences, overlaps):
    """
    Create a sparse graph from the sequences and overlaps information.

    Args:
        sequences (dict): Dictionary containing sequences with sequence IDs as keys.
        overlaps (pandas.DataFrame): DataFrame containing the parsed overlaps information.

    Returns:
        networkx.Graph: A sparse graph representing the overlaps between sequences.
    """
    graph = nx.Graph()
    for seq_id, seq in sequences.items():
        graph.add_node(seq_id, sequence=seq)
    for _, row in overlaps.iterrows():
        query_id = row['query_id']
        target_id = row['target_id']
        if query_id == target_id:
            continue
        query_end = row['query_end']
        target_start = row['target_start']
        if query_end > target_start:
            orientation = '+'
        else:
            orientation = '-'
        weight = query_end - target_start
        graph.add_edge(query_id, target_id, weight=weight, orientation=orientation)

    # Perform coverage-preserving sparsification
    tree = nx.minimum_spanning_tree(graph)
    for node in tree.nodes:
        neighbors = list(tree.neighbors(node))
        if len(neighbors) > 1:
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    if tree.has_edge(neighbors[i], neighbors[j]):
                        tree.remove_edge(neighbors[i], neighbors[j])

    return tree


def create_contigs(graph):
    """
    Create contigs from the sparse graph.

    Args:
        graph (networkx.Graph): A sparse graph representing the overlaps between sequences.

    Returns:
        list: A list of contigs, where each contig is represented as a list of tuples (node_id, sequence).
    """
    visited = set()
    contigs = []

    def dfs(visited, graph, node, contig):
        if node not in visited:
            visited.add(node)
            contig.append((node, graph.nodes[node]['sequence']))
            for neighbor in graph[node]:
                dfs(visited, graph, neighbor, contig)

    for node in graph:
        if node not in visited:
            contig = []
            dfs(visited, graph, node, contig)
            contigs.append(contig)

    return contigs


def merge_contigs(contigs, similarity_threshold=0.95):
    """
    Merge similar contigs based on a specified similarity threshold.

    Args:
        contigs (list): A list of contigs, where each contig is represented as a list of tuples (node_id, sequence).
        similarity_threshold (float, optional): The similarity threshold to consider two contigs as similar.
                                               Defaults to 0.95.

    Returns:
        list: A list of merged contigs, where each merged contig is represented as a tuple (node_id, merged_sequence).
    """
    merged_contigs = []
    contig_count = len(contigs)

    for i in range(contig_count):
        for j in range(i + 1, contig_count):
            contig_i = contigs[i][0][1]
            contig_j = contigs[j][0][1]

            # Calculate the similarity between contigs
            similarity = calculate_similarity(contig_i, contig_j)

            # Merge contigs if the similarity is above the threshold
            if similarity >= similarity_threshold:
                merged_contig = merge_based_on_overlap(contig_i, contig_j)
                merged_contigs.append((contigs[i][0][0], merged_contig))
    return merged_contigs


def jaccard_similarity(contig_i, contig_j):
    """
    Calculate the Jaccard similarity between two contigs.

    Args:
        contig_i (str): Sequence of the first contig.
        contig_j (str): Sequence of the second contig.

    Returns:
        float: Jaccard similarity score.
    """
    set_i = set(contig_i)
    set_j = set(contig_j)
    intersection = set_i.intersection(set_j)
    union = set_i.union(set_j)
    return len(intersection) / len(union)


def calculate_similarity(contig_i, contig_j):
    """
    Calculate the similarity between two contigs (wrapper function for Jaccard similarity).

    Args:
        contig_i (str): Sequence of the first contig.
        contig_j (str): Sequence of the second contig.

    Returns:
        float: Similarity score between 0 and 1.
    """
    return jaccard_similarity(contig_i, contig_j)


def merge_based_on_overlap(seq1, seq2, similarity_threshold=0.5):
    """
    Merge two sequences based on overlap if their similarity is above the threshold.

    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        similarity_threshold (float, optional): The similarity threshold to consider the sequences as overlapping.
                                               Defaults to 0.5.

    Returns:
        str or None: Merged sequence if similarity is above the threshold, else None.
    """
    similarity = calculate_similarity(seq1, seq2)

    if similarity >= similarity_threshold:
        overlap_start = 0
        for i in range(len(seq1)):
            if seq1[i:] == seq2[:len(seq1) - i]:
                overlap_start = i
                break

        merged_seq = seq1[:overlap_start] + seq2
        return merged_seq
    else:
        return None

def main():
    fastq_file_path = 'foo-reads.fq'
    paf_file_path = 'foo.paf'
    output_file_path = 'output.fasta'

    # Start the timer to measure execution time
    start = time.time()

    # Read sequences from the FastQ file
    sequences = read_fastq_file(fastq_file_path)

    # Parse overlaps information from the PAF file
    overlaps = parse_paf(paf_file_path)

    # Create a sparse graph representing the overlaps between sequences
    tree = sparsify_graph(sequences, overlaps)

    # Generate contigs from the sparse graph
    contigs = create_contigs(tree)

    # Merge similar contigs based on a specified similarity threshold
    merged_contigs = merge_contigs(contigs)

    # Print the number of merged contigs
    print(len(merged_contigs))

    # End the timer and calculate execution time
    end = time.time()
    execution_time = end - start
    print(f"Execution time: {execution_time:.2f} seconds")


if __name__ == '__main__':
    main()
