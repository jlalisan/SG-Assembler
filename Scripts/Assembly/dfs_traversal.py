#!/usr/bin/env python3

"""
Depth-First Search (DFS) Traversal

This module provides a depth-first search (DFS) traversal function for a directed multigraph.

Author: Lisan Eisinga
Date: 26-07-2023
"""

import networkx as nx


def dfs(graph: nx.MultiDiGraph) -> list[list]:
    """
    Perform depth-first search traversal on the graph.

    Args:
        graph: A directed multigraph (instance of nx.MultiDiGraph).

    Returns:
        A list of contigs (paths) in the graph with the correct read order.

    Raises:
        TypeError: If the input graph is not an instance of nx.MultiDiGraph.
    """
    try:
        if not isinstance(graph, nx.MultiDiGraph):
            raise TypeError("graph must be an instance of nx.MultiDiGraph")

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

        # Ensure the graph is not empty
        if not graph:
            return []

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

    except TypeError as type_error:
        raise TypeError(f"Invalid graph type: {str(type_error)}") from type_error
