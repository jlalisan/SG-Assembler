#!/usr/bin/env python

"""
Usage:
  Assembly.py <fastq_file> [--paf_file=<paf_file>] [--output_file=<output_file>] [--graphics]

Options:
  -h --help                    Show this help message.
  --paf_file=<paf_file>         Path to the PAF file [default: None]
  --output_file=<output_file>   Path to the output FASTA file [default: None]
  --graphics                   Display information using the graphics.py script [default: False]

Author: Lisan Eisinga
Version: 7.5
Date of Completion: 2023-05-27

Input file formats:
- fastq_file: FastQ file containing biological sequence and quality score information.
"""

# Imports
import sys
import os
import subprocess
import time
import logging
import networkx as nx
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from docopt import docopt
import visualise


def read_fastq_file(fastq_file):
    """
    Reads a FastQ file and returns a dictionary with the sequence titles as keys
    and the sequences as values.

    Parameters:
        fastq_file (str): The path to the FastQ file.

    Returns:
        dict: A dictionary with the sequence titles as keys and the sequences as values.

    Raises:
        FileNotFoundError: If the FastQ file does not exist.
        Exception: If there is an error while reading the FastQ file.
    """
    if not os.path.exists(fastq_file):
        raise FileNotFoundError(f"FastQ file {fastq_file} not found.")

    sequences = {}
    with open(fastq_file) as fastq:
        for title, seq, qual in FastqGeneralIterator(fastq):
            # Use the first element in order to avoid possible duplicated read_ids
            sequences[title.split()[0]] = seq

    return sequences


def create_paf(fastq_file):
    """
    Create a PAF file from a FastQ file using the "minimap2" command-line tool.
    If the PAF file already exists, it will not be recreated.

    Parameters:
        fastq_file (str): The path to the FastQ file.

    Returns:
        str: The path to the PAF file.

    Raises:
        FileNotFoundError: If the FastQ file does not exist.
        subprocess.CalledProcessError: If the "minimap2" command fails.
        Exception: If there is an error while creating the PAF file.
    """
    try:
        if not os.path.exists(fastq_file):
            raise FileNotFoundError(f"FastQ file {fastq_file} not found.")

        # Generate the PAF file name based on the FastQ file name
        paf_file = f"{os.path.splitext(fastq_file)[0]}_overlaps.paf"

        # Check if the PAF file already exists
        if os.path.exists(paf_file):
            logging.info("PAF file %s already exists.", paf_file)
            return paf_file

        # Create the PAF file
        minimap_command = f"./minimap2/minimap2 -x --ava-ont {fastq_file} {fastq_file} > {paf_file}"
        subprocess.run(minimap_command, shell=True, check=True)
        logging.info("PAF file %s created.", paf_file)

        return paf_file

    except FileNotFoundError as error:
        logging.error("FileNotFoundError: %s", str(error))
        raise

    except subprocess.CalledProcessError as error:
        logging.error("CalledProcessError: %s", str(error))
        raise


def parse_paf(paf_file):
    """
    Reads a PAF file and returns a pandas DataFrame with the relevant columns.

    Parameters:
        paf_file (str): The path to the PAF file.

    Returns:
        pandas.DataFrame: A DataFrame with columns "query_id", "query_length", "query_start",
        "query_end", "strand", "target_id", "target_length", "target_start", and "target_end".

    Raises:
        FileNotFoundError: If the PAF file does not exist.
        Exception: If there is an error while parsing the PAF file.
    """
    if not os.path.exists(paf_file):
        raise FileNotFoundError(f"PAF file {paf_file} not found.")

    # Read the PAF file into a DataFrame
    columns = ["query_id", "query_length", "query_start", "query_end", "strand",
               "target_id", "target_length", "target_start", "target_end"]
    paf_df = pd.read_csv(paf_file, sep="\t", header=None, usecols=range(9), names=columns)
    logging.info("PAF file %s parsed.", paf_file)

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
                                "query_start", "query_end", "target_start", "target_end", "strand"]
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

            # Calculate overlap length based on the given information
            overlap_len = min(query_end - query_start + 1, target_end - target_start + 1)

            # Create a dictionary of attributes for the edge
            edge_attrs = {
                "query_start": query_start,
                "query_end": query_end,
                "target_start": target_start,
                "target_end": target_end,
                "strand": strand,
                "query_length": query_len,
                "target_length": target_len,
                "overlap_len": overlap_len  # Ensure "overlap_len" attribute is present
            }

            if query_id != target_id:
                graph.add_edge(query_id, target_id, **edge_attrs)

        return graph

    except KeyError as error:
        raise KeyError(f"Missing required columns in the overlaps DataFrame: {str(error)}")
    except ValueError as error:
        raise ValueError(f"Sequence ID not found in the sequences dictionary: {str(error)}")


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

    except TypeError as error:
        raise TypeError(f"Invalid graph type: {str(error)}")


def dfs(graph):
    """
    Perform depth-first search traversal on the graph.

    Parameters:
        graph (nx.DiGraph or nx.MultiDiGraph): A directed graph.

    Returns:
        list: A list of contigs (paths) in the graph with correct read order.
    """
    try:
        if not isinstance(graph, nx.MultiDiGraph):
            raise TypeError("graph must be of type nx.MultiDiGraph")

        def dfs_visit(node, path, visited):
            path.append(node)  # Append the current node to the path
            visited.add(node)  # Mark the current node as visited

            # Explore outgoing edges based on alignment information
            outgoing_edges = sorted(graph.out_edges(node, keys=True),
                                    key=lambda edge: graph.edges[edge]['target_start'])
            for _, child, _ in outgoing_edges:
                if child not in visited:
                    dfs_visit(child, path, visited)

            # Explore incoming edges based on alignment information
            incoming_edges = sorted(graph.in_edges(node, keys=True),
                                    key=lambda edge: graph.edges[edge]['query_start'], reverse=True)
            for parent, _, _ in incoming_edges:
                if parent not in visited:
                    dfs_visit(parent, path, visited)

        contigs = []
        visited = set()  # Track visited nodes to handle disconnected components

        # Find nodes with no incoming edges
        start_nodes = [node for node in graph.nodes if not list(graph.in_edges(node))]

        for start_node in start_nodes:
            if start_node not in visited:
                path = []  # Store the current path for each component
                dfs_visit(start_node, path, visited)
                contigs.append(path)

        return contigs

    except TypeError as error:
        raise TypeError(f"Invalid graph type: {str(error)}")


def get_consensus_sequences(graph, contigs):
    """ Creates the contig sequences based on positions. """
    if not isinstance(graph, nx.MultiDiGraph):
        raise TypeError("Invalid graph type. Expected MultiDiGraph.")

    if not contigs:
        raise ValueError("Empty contigs list.")

    consensus_sequences = []
    for contig in contigs:
        sequence = []
        covered_positions = set()
        significant_nodes = []
        insignificant_nodes = []
        prev_node = None
        for i, node in enumerate(contig):
            if node in graph.nodes:
                node_sequence = graph.nodes[node]["sequence"]
                edge_data = None
                if i > 0:
                    edge_data = graph.get_edge_data(prev_node, node)
                    if edge_data is not None:
                        target_start = edge_data[0].get("target_start", 0)
                        target_end = edge_data[0].get("target_end", len(graph.nodes[prev_node]["sequence"]))
                        query_start = edge_data[0].get("query_start", 0)
                        query_end = edge_data[0].get("query_end", len(node_sequence))
                        if target_end >= target_start and query_end >= query_start:
                            significant_nodes.append((node, target_start, target_end))
                            if node == contig[0]:
                                uncovered_positions = range(target_end)
                            else:
                                uncovered_positions = [
                                    j
                                    for j in range(len(node_sequence))
                                    if j not in covered_positions
                                       and target_start <= j < target_end
                                       and query_start <= j < query_end
                                ]
                            sequence.extend(node_sequence[pos] for pos in uncovered_positions)
                            covered_positions.update(uncovered_positions)
                        else:
                            insignificant_nodes.append(node)
                    else:
                        insignificant_nodes.append(node)
                else:
                    significant_nodes.append((node, 0, len(node_sequence)))
                    sequence.extend(node_sequence)
                    covered_positions.update(range(len(node_sequence)))
            else:
                continue

            prev_node = node  # Update prev_node for the next iteration

        significant_nodes.sort(key=lambda x: x[1])  # Sort significant nodes based on target start

        for node, target_start, target_end in significant_nodes:
            edge_data = graph.get_edge_data(prev_node, node)
            if edge_data is not None:
                query_start = edge_data[0].get("query_start", len(graph.nodes[node]["sequence"]))
                query_end = edge_data[0].get("query_end", len(graph.nodes[node]["sequence"]))
                if node == contig[0]:
                    significant_sequence = graph.nodes[node]["sequence"][:query_end]
                    sequence = list(significant_sequence) + sequence
                else:
                    significant_sequence = graph.nodes[node]["sequence"][target_start:query_start]
                    insert_index = 0 if prev_node == contig[0] else sequence.index(
                        graph.nodes[prev_node]["sequence"][-1]) + target_start
                    sequence = list(significant_sequence) + sequence[:insert_index] + sequence[insert_index:]

        if sequence:
            consensus_sequences.append("".join(sequence))

    return consensus_sequences


def write_to_file(filename, contigs):
    """ Writes the contigs to a Fasta file."""
    try:
        with open(filename, 'w') as file:
            for i, contig in enumerate(contigs):
                header = f">contig_{i + 1}"
                sequence = contig
                file.write(f"{header}\n{sequence}\n")
    except IOError:
        logging.info("Error: Unable to write to file %s", filename)


def main():
    """
    Perform the assembly process based on the provided command-line arguments.
    """
    arguments = docopt(__doc__)
    fastq_file = arguments["<fastq_file>"]
    paf_file = arguments["--paf_file"]
    output_file = arguments["--output_file"]
    graphics = arguments["--graphics"]

    start = time.time()
    if not paf_file:
        paf_file = create_paf(fastq_file)

    logging.info("Reading FastQ file...")
    sequences = read_fastq_file(fastq_file)

    logging.info("Parsing PAF file...")
    overlaps = parse_paf(paf_file)

    logging.info("Creating overlap graph...")
    mygraph = overlap_graph(sequences, overlaps)
    mygraph = remove_isolated_nodes(mygraph)

    if mygraph.number_of_nodes() >= 999:
        try:
            # Set recursion limit to avoid stack overflow error
            sys.setrecursionlimit(mygraph.number_of_nodes())
        except RecursionError:
            logging.info("Error: Recursion limit cannot be set. The graph is too large for traversal.")
            return

    logging.info("Traversing the graph...")
    try:
        contigs = dfs(mygraph)
    except RecursionError:
        logging.info("Error: Recursion limit exceeded. Unable to traverse the graph.")
        return

    if contigs:
        logging.info("Generating consensus sequences...")
        consensus_seqs = get_consensus_sequences(mygraph, contigs)

        if not output_file:
            base_name = os.path.splitext(fastq_file)[0]
            output_file = f"{base_name}_contigs.fasta"

        logging.info("Writing contigs to %s...", output_file)
        write_to_file(consensus_seqs, output_file)

    if graphics:
        # Call the graphics functions from graphics.py
        visualise.plot_read_lengths(fastq_file)
        visualise.plot_gc_content(fastq_file)
        visualise.plot_quality_scores(fastq_file)
        visualise.plot_sequence_complexity(fastq_file)
        visualise.count_duplicates(fastq_file)

    end = time.time()
    logging.info("Assembly was finished in: %s seconds",  end - start)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
