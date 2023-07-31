#!/usr/bin/env python3
"""
De Bruijn Graph-Based Assembler

This script reads FASTQ data, constructs a De Bruijn graph, merges overlapping k-mers,
traverses the graph using an Eulerian path, and writes the assembled contigs to a file.

Usage:
    1. Update the 'fastq_file' variable with the name of your input FASTQ file.
    2. Run the script using: Python de_bruijn_assembler.py

Author: Lisan Eisinga
"""

import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def read_fastq_file(fastq_file):
    """
    Read a FASTQ file and return a dictionary of sequences.

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
            # Use the first element in order to avoid possible duplicated read_ids
            sequences[title.split()[0]] = seq

    return sequences


def generate_kmers(sequence, k):
    """
    Generate k-mers from a given sequence.

    Parameters:
        sequence (str): The input DNA sequence.
        k (int): The length of k-mers.

    Yields:
        str: K-mers of length k.
    """
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]


def generate_de_bruijn_graph(sequences, k, threshold=2):
    """
    Generate a De Bruijn graph from a set of sequences.

    Parameters:
        sequences (dict): A dictionary with read IDs as keys and sequences as values.
        k (int): The length of k-mers for De Bruijn graph construction.
        threshold (int): Minimum count threshold to include a k-mer in the graph.

    Returns:
        dict: The De Bruijn graph represented as an adjacency list.
    """
    de_bruijn_graph = {}

    # Step 1: Count the occurrences of each k-mer
    kmer_counts = {}
    for read_id, sequence in sequences.items():
        for kmer in generate_kmers(sequence, k):
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1

    # Step 2: Build the sparse De Bruijn graph
    for kmer, count in kmer_counts.items():
        if count >= threshold:
            kmer_prefix = kmer[:-1]
            kmer_suffix = kmer[1:]
            if kmer_prefix not in de_bruijn_graph:
                de_bruijn_graph[kmer_prefix] = []
            de_bruijn_graph[kmer_prefix].append(kmer_suffix)

    return de_bruijn_graph


def merge_overlapping_kmers(de_bruijn_graph):
    """
    Merge overlapping k-mers in the De Bruijn graph.

    Parameters:
        de_bruijn_graph (dict): The De Bruijn graph represented as an adjacency list.

    Returns:
        dict: The condensed De Bruijn graph after merging overlapping k-mers.
    """
    condensed_graph = de_bruijn_graph.copy()

    # Create a set of all nodes with incoming edges (in-nodes)
    in_nodes = set(node for nodes in condensed_graph.values() for node in nodes)

    # Function to find all overlaps for a given k-mer
    def find_overlaps(kmer):
        overlaps = set()
        for node in condensed_graph[kmer]:
            for i in range(1, len(node) + 1):
                overlap = node[:i]
                if kmer.endswith(overlap):
                    overlaps.add((kmer, node, overlap))
        return overlaps

    # Iterate through the nodes in the De Bruijn graph
    for kmer in de_bruijn_graph:
        if kmer not in in_nodes:
            # Find all overlaps for the current k-mer
            overlaps = find_overlaps(kmer)

            while overlaps:
                # Merge the nodes with the longest overlap
                kmer, node, overlap = max(overlaps, key=lambda o: len(o[2]))

                # Merge the nodes by extending the k-mer with the node suffix
                merged_kmer = kmer + node[len(overlap):]

                # Update the graph to include the merged k-mer
                condensed_graph[kmer].remove(node)
                if not condensed_graph[kmer]:
                    del condensed_graph[kmer]
                if merged_kmer not in condensed_graph:
                    condensed_graph[merged_kmer] = []
                if node in condensed_graph:
                    condensed_graph[merged_kmer].extend(condensed_graph[node])
                    del condensed_graph[node]

                # Find new overlaps for the merged k-mer
                overlaps = find_overlaps(merged_kmer)

    return condensed_graph


def traverse_graph_eulerian_path(condensed_graph):
    """
    Traverse the De Bruijn graph using an Eulerian path and generate contigs.

    Parameters:
        condensed_graph (dict): The condensed De Bruijn graph after merging overlapping k-mers.

    Returns:
        list: A list of contigs (DNA sequences).
    """
    def find_eulerian_path(graph, start):
        stack = [start]
        eulerian_paths = []

        while stack:
            node = stack[-1]
            if graph[node]:
                stack.append(graph[node].pop())
            else:
                eulerian_paths.append(stack[:])
                stack.pop()

        return eulerian_paths

    # Find all possible Eulerian paths in the graph starting from each node
    eulerian_paths = []
    for start_node in condensed_graph:
        paths = find_eulerian_path(condensed_graph.copy(), start_node)
        eulerian_paths.extend(paths)

    # Combine nodes in each path to form contigs
    contigs = []
    for path in eulerian_paths:
        contig = path[0]
        for node in path[1:]:
            contig += node[-1]
        contigs.append(contig)

    return contigs


def write_contigs_to_file(contigs, output_file):
    """
    Write contigs to a file in FASTA format.

    Parameters:
        contigs (list): A list of contigs (DNA sequences).
        output_file (str): Path to the output file where contigs will be written.
    """
    with open(output_file, "w") as f:
        for i, contig in enumerate(contigs, start=1):
            f.write(f">Contig_{i}\n")
            f.write(contig + "\n")


def main():
    """
    Main function to assemble contigs from FASTQ data using the De Bruijn graph.

    Reads FASTQ data, constructs De Bruijn graph, merges overlapping k-mers,
    traverses the graph using an Eulerian path, and writes the contigs to a file.

    Change the 'fastq_file' and 'output_file' variables as needed.
    """
    k = 31
    input_file = '../../Test_data/foo-reads.fq'
    output_file = "../../Test_data/debruijn_contigs.fasta"

    sequences = read_fastq_file(input_file)
    de_bruijn_graph = generate_de_bruijn_graph(sequences, k)
    print("Before merging:", len(de_bruijn_graph))
    newgraph = merge_overlapping_kmers(de_bruijn_graph)
    print("After merging:", len(newgraph))
    contigs = traverse_graph_eulerian_path(newgraph)
    write_contigs_to_file(contigs, output_file)


if __name__ == "__main__":
    main()

