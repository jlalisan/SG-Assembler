"""
Consensus Sequence Generation

This script provides a function to generate consensus sequences based on the overlap graph,
and the DFS contigs.

Author: Lisan Eisinga
Date: July 26, 2023
"""

import networkx as nx


def generate_consensus_sequences(graph: nx.MultiDiGraph, contigs: list[list]) -> list[str]:
    """
    Generate consensus sequences based on the overlap graph and DFS contigs.

    Args:
        graph: A directed multigraph representing the overlap graph (instance of nx.MultiDiGraph).
        contigs: A list of contigs (paths) in the graph with correct read order.

    Returns:
        A list of generated consensus sequences.

    Raises:
        TypeError: If the input graph is not an instance of nx.MultiDiGraph,
        or if contigs is not a list.
    """
    try:
        if not isinstance(graph, nx.MultiDiGraph):
            raise TypeError("graph must be an instance of nx.MultiDiGraph")

        if not isinstance(contigs, list):
            raise TypeError("contigs must be a list")

        sequences = []
        for contig in contigs:
            consensus_sequence = ""
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

                # Append the relevant portion of the node's sequence to the consensus sequence
                if alignment_block_length > 0:
                    consensus_sequence += node_sequence[overlap_len:overlap_len + alignment_block_length]
                else:
                    consensus_sequence += node_sequence

                # Update the previous node
                prev_node = node

            sequences.append(consensus_sequence)

        return sequences

    except TypeError as type_error:
        raise TypeError(f"Invalid argument type: {str(type_error)}") from type_error
