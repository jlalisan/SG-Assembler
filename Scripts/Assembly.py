#!/usr/bin/env python

"""
Usage:
  assembly.py <fastq_file> [--paf_file=<paf_file>] [--output_file=<output_file>] [--graphics]
  assembly.py (-h | --help)

Options:
  -h --help                     Show this help message.
  --paf_file=<paf_file>         Path to the PAF data set.
  --output_file=<output_file>   Path to the output Fasta data set.
  --graphics                    Display information using the graphics.py script. [default: False]

Author: Lisan Eisinga
Version: 7.5
Date of Completion: 11-06-2023

Input data set formats:
- fastq_file: FastQ data set containing biological sequence and quality score information (MinION).
"""

# Imports
import sys
import os
import subprocess
import logging
import cProfile

import networkx as nx
import pandas as pd
import psutil
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from docopt import docopt
import visualise


def read_fastq_file(fastq_file):
    """
    Reads a FastQ data set and returns a dictionary with the sequence titles as keys
    and the sequences as values.

    Parameters:
        fastq_file (str): The path to the FastQ data set.

    Returns:
        dict: A dictionary with the sequence titles as keys and the sequences as values.

    Raises:
        FileNotFoundError: If the FastQ data set does not exist.
        Exception: If there is an error while reading the FastQ data set.
    """
    if not os.path.exists(fastq_file):
        raise FileNotFoundError(f"FastQ data set {fastq_file} not found.")

    sequences = {}
    with open(fastq_file, "r", encoding="utf-8") as fastq:
        for title, seq, _ in FastqGeneralIterator(fastq):
            # Use the first element in order to avoid possible duplicated read_ids
            sequences[title.split()[0]] = seq

    return sequences


def create_paf(fastq_file):
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
        if not os.path.exists(fastq_file):
            raise FileNotFoundError(f"FastQ data set {fastq_file} not found.")

        # Generate the PAF data set name based on the FastQ data set name
        paf_file = f"{os.path.splitext(fastq_file)[0]}_overlaps.paf"

        # Check if the PAF data set already exists
        if os.path.exists(paf_file):
            logging.info("PAF data set %s already exists.", paf_file)
            return paf_file

        # Create the PAF data set
        minimap_command = f"./minimap2/minimap2 -x ava-ont {fastq_file} {fastq_file} > {paf_file}"
        subprocess.run(minimap_command, shell=True, check=True)
        logging.info("PAF data set %s created.", paf_file)

        return paf_file

    except FileNotFoundError as error:
        logging.error("FileNotFoundError: %s", str(error))
        raise

    except subprocess.CalledProcessError as error:
        logging.error("CalledProcessError: %s", str(error))
        raise


def parse_paf(paf_file):
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
        raise FileNotFoundError(f"PAF data set {paf_file} not found.")

    # Read the PAF data set into a DataFrame
    columns = ["query_id", "query_length", "query_start", "query_end", "strand",
               "target_id", "target_length", "target_start", "target_end",
               "alignment_block_length", "residue_matches", "mapping_quality"]
    paf_df = pd.read_csv(paf_file, sep="\t", header=None, usecols=range(12), names=columns)
    return paf_df


def overlap_graph(sequences, overlaps):
    """
    Create a directed multigraph from sequences and overlaps.

    Parameters:
        sequences (dict): A dictionary of sequences.
        overlaps (pandas.DataFrame): A pandas DataFrame of overlaps.

    Returns:
        nx.MultiDiGraph: A directed multigraph.

    Raises:
        KeyError: If any required columns are missing in the overlaps DataFrame.
        ValueError: If the sequence IDs in the overlaps DataFrame are not present
        in the sequences dictionary.
    """
    try:
        # Create a directed multigraph
        graph = nx.MultiDiGraph()

        # Add nodes to the graph
        for seq_id in sequences:
            graph.add_node(seq_id, sequence=sequences[seq_id])

        # Add edges to the overlap graph
        for _, row in overlaps.iterrows():
            # Check if all required columns are present in the overlaps DataFrame
            required_columns = ["query_id", "target_id", "query_length", "target_length",
                                "query_start", "query_end", "target_start", "target_end",
                                "strand", "alignment_block_length", "residue_matches",
                                "mapping_quality"]
            if not all(column in row for column in required_columns):
                raise KeyError("Missing required columns in the overlaps DataFrame.")

            query_id = row["query_id"]
            target_id = row["target_id"]

            # Check if the sequence IDs in the overlaps DataFrame are present in the dict
            if query_id not in sequences or target_id not in sequences:
                raise ValueError("Sequence ID not found in the sequences dictionary.")

            query_len = row["query_length"]
            target_len = row["target_length"]
            query_start = row["query_start"]
            query_end = row["query_end"]
            target_start = row["target_start"]
            target_end = row["target_end"]
            strand = row["strand"]
            alignment_block_length = row["alignment_block_length"]
            residue_matches = row["residue_matches"]
            mapping_quality = row["mapping_quality"]

            # Calculate overlap length based on the given information and strand
            overlap_len = 0

            if strand == "+":
                overlap_len = min(query_end - query_start + 1, target_end - target_start + 1)
            elif strand == "-":
                overlap_len = min(query_end - query_start + 1, target_start - target_end + 1)

            # Adjust overlap length if negative or greater than the minimum of query/target lengths
            overlap_len = max(0, min(overlap_len, min(query_len, target_len)))

            # Create a dictionary of attributes for the edge
            edge_attrs = {
                "query_start": query_start,
                "query_end": query_end,
                "target_start": target_start,
                "target_end": target_end,
                "strand": strand,
                "query_length": query_len,
                "target_length": target_len,
                "overlap_len": overlap_len,
                "alignment_block_length": alignment_block_length,
                "residue_matches": residue_matches,
                "mapping_quality": mapping_quality
            }

            if query_id != target_id:
                if not graph.has_edge(query_id, target_id):
                    graph.add_edge(query_id, target_id)
                graph.edges[query_id, target_id, 0].update(edge_attrs)

        return graph

    except KeyError as k_e:
        raise KeyError(f"Missing required columns in the overlaps DataFrame: {str(k_e)}") from k_e
    except ValueError as v_e:
        raise ValueError(f"Sequence ID not found in the sequences dictionary: {str(v_e)}") from v_e


def remove_isolated_nodes(graph):
    """
    Removes isolated nodes from a graph.

    Parameters:
        graph (nx.MultiDiGraph): A directed multigraph.

    Returns:
        nx.MultiDiGraph: The graph with isolated nodes removed.

    Raises:
        TypeError: If the input graph is not an instance of nx.MultiDiGraph.
    """
    try:
        if not isinstance(graph, nx.MultiDiGraph):
            raise TypeError("Input graph must be an instance of nx.MultiDiGraph.")

        # Get a list of isolated nodes
        isolated_nodes = list(nx.isolates(graph))

        # Remove the isolated nodes from the graph
        graph.remove_nodes_from(isolated_nodes)

        return graph

    except TypeError as t_e:
        raise TypeError(f"Invalid graph type: {str(t_e)}") from t_e


def dfs(graph):
    """
    Perform depth-first search traversal on the graph.

    Parameters:
        graph (nx.MultiDiGraph): A directed multi graph.

    Returns:
        list: A list of contigs (paths) in the graph with correct read order.
    """
    try:
        if not isinstance(graph, nx.MultiDiGraph):
            raise TypeError("graph must be of type nx.MultiDiGraph")

        def dfs_visit(node, path, visited):
            # Append the current node to the path
            path.append(node)
            # Mark the current node as visited
            visited.add(node)

            # Explore outgoing edges based on alignment information
            outgoing_edges = sorted(graph.out_edges(node, keys=True),
                                    key=lambda edge: (
                                        graph.edges[edge]["target_start"],
                                        -graph.edges[edge]["target_end"],
                                        -graph.edges[edge]["query_start"],
                                        graph.edges[edge]["query_end"]
                                    ))

            for _, child, _ in outgoing_edges:
                if child not in visited:
                    dfs_visit(child, path, visited)

            # Explore incoming edges based on alignment information
            incoming_edges = sorted(graph.in_edges(node, keys=True),
                                    key=lambda edge: (
                                        -graph.edges[edge]["query_start"],
                                        graph.edges[edge]["query_end"],
                                        graph.edges[edge]["target_start"],
                                        -graph.edges[edge]["target_end"]
                                    ))
            for parent, _, _ in incoming_edges:
                if parent not in visited:
                    dfs_visit(parent, path, visited)

        contigs = []
        # Track visited nodes to handle disconnected components
        visited = set()

        while len(visited) < len(graph.nodes):
            # Find the node with the most outgoing edges among unvisited nodes
            max_outgoing_edges = 0
            start_node = None
            for node in graph.nodes:
                if node not in visited:
                    outgoing_edges = graph.out_edges(node)
                    num_outgoing_edges = len(outgoing_edges)
                    if num_outgoing_edges > max_outgoing_edges:
                        max_outgoing_edges = num_outgoing_edges
                        start_node = node

            if start_node is not None:
                # Store the current path for each component
                path = []
                dfs_visit(start_node, path, visited)
                # Reverse the path to get the correct read order
                contigs.append(path)

        return contigs

    except TypeError as t_e:
        raise TypeError(f"Invalid graph type: {str(t_e)}") from t_e


def generate_sequence(graph, contigs):
    """
    Generate sequences based on the overlap graph and DFS contigs.

    Parameters:
        graph (nx.MultiDiGraph): A directed multigraph representing the overlap graph.
        contigs (list): A list of contigs (paths) in the graph with correct read order.

    Returns:
        list: A list of generated sequences.
    """
    try:
        if not isinstance(graph, nx.MultiDiGraph):
            raise TypeError("graph must be of type nx.MultiDiGraph")

        if not isinstance(contigs, list):
            raise TypeError("contigs must be a list")

        sequences = []
        for contig in contigs:
            sequence = ""
            prev_node = None
            for node in contig:
                # Check if the node is covered by another node
                if any(graph.out_edges(node)):
                    continue  # Skip if the node has outgoing edges (covered by another node)

                # Get the sequence of the node
                node_sequence = graph.nodes[node]["sequence"]
                overlap_len = 0
                alignment_block_length = 0

                # Get the overlap length and alignment block length with the previous node
                if prev_node is not None:
                    edge_attrs = graph.get_edge_data(prev_node, node)
                    if edge_attrs is not None:
                        overlap_len = edge_attrs[0].get("overlap_len", 0)
                        alignment_block_length = edge_attrs[0].get("alignment_block_length", 0)

                        # Adjust overlap length for reverse complement sequences
                        if edge_attrs[0].get("strand") == "-":
                            overlap_len = -(overlap_len + alignment_block_length)

                # Append the relevant portion of the node's sequence to the overall sequence
                if alignment_block_length > 0:
                    sequence += node_sequence[overlap_len:overlap_len + alignment_block_length]
                else:
                    sequence += node_sequence

                # Update the previous node
                prev_node = node

            sequences.append(sequence)

        return sequences

    except TypeError as t_e:
        raise TypeError(f"Invalid argument type: {str(t_e)}") from t_e


def write_to_file(filename, contigs):
    """
    Writes the contigs to a Fasta data set with individual contig numbering.

    Parameters:
        filename (str): The path to the output file.
        contigs (list): A list of contigs to be written.

    Raises:
        IOError: If there is an error while writing to the data set.

    Returns:
        None
    """
    try:
        with open(filename, "w", encoding="utf-8") as file:
            for i, contig in enumerate(contigs):
                header = f">contig_{i + 1}"
                sequence = "".join(contig)
                file.write(f"{header}\n{sequence}\n")
    except IOError:
        logging.info("Error: Unable to write to data set %s", filename)


def main(args):
    """
    Perform the assembly process based on the provided command-line arguments.
    """
    process = psutil.Process()
    memory_before = process.memory_info().rss / (1024 * 1024)

    fastq_file = args["<fastq_file>"]
    paf_file = args["--paf_file"]
    output_file = args["--output_file"]
    graphics = args["--graphics"]

    if not paf_file:
        paf_file = create_paf(fastq_file)

    logging.info("Reading FastQ data set...")
    sequences = read_fastq_file(fastq_file)

    logging.info("Parsing PAF data set...")
    overlaps = parse_paf(paf_file)

    logging.info("Creating overlap graph...")
    mygraph = overlap_graph(sequences, overlaps)
    mygraph = remove_isolated_nodes(mygraph)

    if mygraph.number_of_nodes() >= 999:
        try:
            # Set recursion limit to avoid stack overflow error
            sys.setrecursionlimit(mygraph.number_of_nodes())
        except RecursionError:
            logging.info("Error: Recursion limit cant be set. The graph is too big for traversal.")
            return

    logging.info("Traversing the graph...")
    try:
        contigs = dfs(mygraph)
    except RecursionError:
        logging.info("Error: Recursion limit exceeded. Unable to traverse the graph.")
        return

    if contigs:
        logging.info("Generating consensus sequences...")
        consensus_seqs = generate_sequence(mygraph, contigs)
        if not output_file:
            base_name = os.path.splitext(fastq_file)[0]
            output_file = f"{base_name}_contigs.fasta"
        logging.info("Writing contigs to %s...", output_file)
        write_to_file(output_file, consensus_seqs)

    if graphics:
        # Call the graphics functions from graphics.py
        visualise.plot_read_lengths(fastq_file)
        visualise.plot_gc_content(fastq_file)
        visualise.plot_quality_scores(fastq_file)
        visualise.plot_sequence_complexity(fastq_file)
        visualise.count_duplicates(fastq_file)

    memory_after = process.memory_info().rss / (1024 * 1024)
    logging.info("Memory usage: %s MB", round((memory_after - memory_before), 2))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    profiler = cProfile.Profile()
    profiler.enable()
    arguments = docopt(__doc__)
    main(arguments)
    profiler.disable()
    profiler.dump_stats("performance.prof")
    logging.info("Performance saved in performance.prof")
