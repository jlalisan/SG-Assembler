#!/usr/bin/env python

"""
Usage:
  Assembly.py <fastq_file> [--paf_file=<paf_file>] [--output_file=<output_file>]

Options:
  -h --help                    Show this help message.
  --paf_file=<paf_file>         Path to the PAF file [default: None]
  --output_file=<output_file>   Path to the output FASTA file [default: None]

Author: Lisan Eisinga
Version: 7.3
Date of Completion: 2023-05-27
"""

import os
import subprocess
import sys
import time
import networkx as nx
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from docopt import docopt


def read_fastq_file(fastq_file):
    """
    Reads a FastQ file and returns a dictionary with the sequence titles as keys and the sequences as values.

    Parameters:
        fastq_file (str): The path to the FastQ file.

    Returns:
        dict: A dictionary with the sequence titles as keys and the sequences as values.
    """
    sequences = {}
    with open(fastq_file) as f:
        for title, seq, qual in FastqGeneralIterator(f):
            # Use the first element in order to avoid possible duplicated read_ids
            sequences[title.split()[0]] = seq
    return sequences


def create_paf(fastq_file):
    """
    Create a PAF file from a FastQ file using the 'minimap2' command-line tool.
    If the PAF file already exists, it will not be recreated.

    Parameters:
        fastq_file (str): The path to the FastQ file.

    Returns:
        str: The path to the PAF file.

    Raises:
        subprocess.CalledProcessError: If the 'minimap2' command fails.

    """
    # Generate the PAF file name based on the FastQ file name
    paf_file = f"{os.path.splitext(fastq_file)[0]}_overlaps.paf"

    # Check if the PAF file already exists
    if os.path.exists(paf_file):
        print(f"PAF file '{paf_file}' already exists.")
        return paf_file

    # Create the PAF file
    minimap_command = f"./minimap2/minimap2 -x --ava-ont {fastq_file} {fastq_file} > {paf_file}"
    subprocess.run(minimap_command, shell=True, check=True)
    print(f"PAF file '{paf_file}' created.")

    return paf_file


def parse_paf(paf_file):
    """
    Reads a PAF file and returns a pandas DataFrame with the relevant columns.

    Parameters:
        paf_file (str): The path to the PAF file.

    Returns:
        pandas.DataFrame: A DataFrame with columns 'query_id', 'query_length', 'query_start', 'query_end', 'strand',
        'target_id', 'target_length', 'target_start', and 'target_end'.
    """
    # Read the PAF file into a DataFrame using pandas
    df = pd.read_csv(paf_file, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], header=None)
    df.columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand',
                  'target_id', 'target_length', 'target_start', 'target_end']
    return df


def overlap_graph(sequences, overlaps):
    """
    Create a directed multigraph from sequences and overlaps.

    Parameters:
        sequences (dict): A dictionary of sequences.
        overlaps (pandas.DataFrame): A pandas DataFrame of overlaps.

    Returns:
        nx.MultiDiGraph: A directed multigraph.
    """
    # Create a directed multigraph
    graph = nx.MultiDiGraph()

    # Add nodes to the graph
    for seq_id in sequences:
        graph.add_node(seq_id, sequence=sequences[seq_id])

    # Add edges to the overlap graph
    for _, row in overlaps.iterrows():
        query_id = row['query_id']
        target_id = row['target_id']
        query_len = row['query_length']
        target_len = row['target_length']
        query_start = row['query_start']
        query_end = row['query_end']
        target_start = row['target_start']
        target_end = row['target_end']
        strand = row['strand']

        # Calculate overlap length based on the given information
        overlap_len = min(query_end - query_start + 1, target_end - target_start + 1)

        # Create a dictionary of attributes for the edge
        edge_attrs = {
            'query_start': query_start,
            'query_end': query_end,
            'target_start': target_start,
            'target_end': target_end,
            'strand': strand,
            'query_length': query_len,
            'target_length': target_len,
            'overlap_len': overlap_len  # Ensure 'overlap_len' attribute is present
        }

        if query_id != target_id:
            graph.add_edge(query_id, target_id, **edge_attrs)

    return graph


def remove_isolated_nodes(graph):
    """
    Removes isolated nodes from a graph.

    Parameters:
        graph (nx.MultiDiGraph): A directed multigraph.

    Returns:
        nx.MultiDiGraph: The graph with isolated nodes removed.
    """
    # Get a list of isolated nodes
    isolated_nodes = list(nx.isolates(graph))

    # Remove the isolated nodes from the graph
    graph.remove_nodes_from(isolated_nodes)

    return graph


def dfs(graph):
    """
    Perform depth-first search traversal on the graph.

    Parameters:
        graph (nx.MultiDiGraph): A directed multigraph.

    Returns:
        list: A list of contigs (paths) in the graph with correct read order.
    """
    def dfs_visit(node, path, visited):
        """
        Recursive helper function for depth-first search traversal.

        Parameters:
            node: The current node being visited.
            path: The current path being constructed.
            visited: A set of visited nodes.

        Returns:
            None
        """
        path.append(node)  # Append the current node to the path
        visited.add(node)  # Mark the current node as visited

        # Explore outgoing edges based on alignment information
        outgoing_edges = sorted(graph.out_edges(node, keys=True), key=lambda edge: graph.edges[edge]['target_start'])
        for _, child, _ in outgoing_edges:
            if child not in visited:
                dfs_visit(child, path, visited)

        # Explore incoming edges based on alignment information
        incoming_edges = sorted(graph.in_edges(node, keys=True), key=lambda edge: graph.edges[edge]['query_start'], reverse=True)
        for _, parent, _ in incoming_edges:
            if parent not in visited:
                dfs_visit(parent, path, visited)

    # Perform DFS traversal on each connected component
    contigs = []
    visited = set()
    for node in graph:
        if node not in visited:
            path = []
            dfs_visit(node, path, visited)
            contigs.append(path)

    return contigs


def get_consensus_sequences(graph, contigs):
    """
    Generate consensus sequences for each contig in the graph.

    Parameters:
        graph (nx.MultiDiGraph): A directed multigraph.
        contigs (list): A list of contigs (paths) in the graph with correct read order.

    Returns:
        list: A list of consensus sequences.

    Note:
        This function assumes that the overlaps in the graph have been calculated and stored as edge attributes.

    """

    consensus_sequences = []

    for contig in contigs:
        consensus_sequence = ""
        prev_end = None

        for node in contig:
            if prev_end is not None:
                overlapping_seq = None
                overlap_len = None

                # Retrieve outgoing edges from the previous end node
                outgoing_edges = list(graph.out_edges(prev_end, keys=True))

                # Retrieve incoming edges to the current node
                incoming_edges = list(graph.in_edges(node, keys=True))

                for _, target, key in outgoing_edges:
                    if (target, node, key) in incoming_edges:
                        # Retrieve the overlap length from the edge attribute
                        overlap_len = graph.edges[prev_end, target, key]['overlap_len']
                        target_seq = graph.nodes[target]['sequence']
                        overlapping_seq = target_seq[:overlap_len]
                        break

                if overlapping_seq is not None:
                    if consensus_sequence.endswith(overlapping_seq):
                        # Append the remaining portion of the current node's sequence after the overlap
                        consensus_sequence += graph.nodes[node]['sequence'][overlap_len:]
                    else:
                        # Concatenate the entire sequence of the current node
                        consensus_sequence += graph.nodes[node]['sequence']
                else:
                    # Concatenate the entire sequence of the current node
                    consensus_sequence += graph.nodes[node]['sequence']

            else:
                # Concatenate the entire sequence of the first node in the contig
                consensus_sequence += graph.nodes[node]['sequence']

            prev_end = node

        consensus_sequences.append(consensus_sequence)

    return consensus_sequences


def write_to_file(consensus_sequences, output_file):
    """
    Write consensus sequences to a FASTA file in the Flye format.

    Parameters:
        consensus_sequences (list): A list of consensus sequences.
        output_file (str): The path to the output FASTA file.

    Returns:
        None
    """

    with open(output_file, 'w') as f:
        for i, sequence in enumerate(consensus_sequences):
            header = f">contig_{i + 1}_length_{len(sequence)}\n"
            f.write(header)
            f.write(sequence)
            f.write('\n')


def main():
    """
    Perform the assembly process based on the provided command-line arguments.
    """
    arguments = docopt(__doc__)
    fastq_file = arguments['<fastq_file>']
    paf_file = arguments['--paf_file']
    output_file = arguments['--output_file']

    start = time.time()
    if not paf_file:
        paf_file = create_paf(fastq_file)

    print("Reading FastQ file...")
    sequences = read_fastq_file(fastq_file)

    print("Parsing PAF file...")
    overlaps = parse_paf(paf_file)

    print("Creating overlap graph...")
    mygraph = overlap_graph(sequences, overlaps)
    mygraph = remove_isolated_nodes(mygraph)

    if mygraph.number_of_nodes() >= 999:
        try:
            # Set recursion limit to avoid stack overflow error
            sys.setrecursionlimit(mygraph.number_of_nodes())
        except RecursionError:
            print("Error: Recursion limit cannot be set. The graph is too large for traversal.")
            return

    print("Traversing the graph...")
    try:
        contigs = dfs(mygraph)
    except RecursionError:
        print("Error: Recursion limit exceeded. Unable to traverse the graph.")
        return

    if contigs:
        print("Generating consensus sequences...")
        consensus_seqs = get_consensus_sequences(mygraph, contigs)

        if not output_file:
            base_name = os.path.splitext(fastq_file)[0]
            output_file = f"{base_name}_contigs.fasta"

        print(f"Writing contigs to {output_file}...")
        write_to_file(consensus_seqs, output_file)

    end = time.time()
    print(f"Assembly was finished in: {end - start} seconds")


if __name__ == '__main__':
    main()
