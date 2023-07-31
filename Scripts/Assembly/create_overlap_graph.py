#!/usr/bin/env python3
"""
Overlap Graph Creation and Isolated Nodes Removal

This module contains functions to create a directed multigraph representing overlaps
between sequences and to remove isolated nodes from the graph.

Author: Lisan Eisinga
Date: 26-07-2023
"""

import networkx as nx
import pandas as pd


def overlap_graph(sequences: dict, overlaps: pd.DataFrame) -> nx.MultiDiGraph:
    """
    Create a directed multigraph from sequences and overlaps.

    Parameters:
        sequences: A dictionary of sequences.
        overlaps: A pandas DataFrame of overlaps.

    Returns:
        A directed multigraph representing the overlaps between sequences.

    Raises:
        KeyError: If any required columns are missing in the overlaps DataFrame.
        ValueError: If the sequence IDs in the overlaps DataFrame
                    are not present in the sequences dictionary.
    """
    try:
        # Create a directed multigraph
        graph = nx.MultiDiGraph()

        # Add nodes to the graph
        for seq_id, sequence in sequences.items():
            graph.add_node(seq_id, sequence=sequence)

        # Add edges to the overlap graph
        for _, row in overlaps.iterrows():
            # Check if all required columns are present in the overlaps DataFrame
            required_columns = [
                "query_id", "target_id", "query_length", "target_length",
                "query_start", "query_end", "target_start", "target_end",
                "strand", "alignment_block_length", "residue_matches",
                "mapping_quality"
            ]
            if not all(column in row for column in required_columns):
                raise KeyError("Missing required columns in the overlaps DataFrame.")

            query_id = row["query_id"]
            target_id = row["target_id"]

            # Check if the sequence IDs in the overlaps DataFrame are present in the dictionary
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
            overlap_length = 0

            if strand == "+":
                overlap_length = min(query_end - query_start + 1, target_end - target_start + 1)
            elif strand == "-":
                overlap_length = min(query_end - query_start + 1, target_start - target_end + 1)

            # Adjust overlap length if negative or greater than the minimum of query/target lengths
            overlap_length = max(0, min(overlap_length, min(query_len, target_len)))

            # Create a dictionary of attributes for the edge
            edge_attrs = {
                "query_start": query_start,
                "query_end": query_end,
                "target_start": target_start,
                "target_end": target_end,
                "strand": strand,
                "query_length": query_len,
                "target_length": target_len,
                "overlap_length": overlap_length,
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


def remove_isolated_nodes(graph: nx.MultiDiGraph) -> nx.MultiDiGraph:
    """
    Removes isolated nodes from a graph.

    Parameters:
        graph: A directed multigraph.

    Returns:
        The graph with isolated nodes removed.

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
