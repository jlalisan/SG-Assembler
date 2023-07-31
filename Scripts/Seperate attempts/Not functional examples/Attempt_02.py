#!/usr/bin/env python3
"""
Unfinished Python Script

This script contains a collection of functions to process sequence data, create an overlap graph,
simplify the graph, create contigs, and merge similar contigs.

Note: This script is unfinished and may require further modification and testing for a complete and
functional implementation.
"""

import time
import networkx as nx
import numpy as np
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import community as community_louvain
from matplotlib import pyplot as plt
import igraph as ig


def read_fastq_file(fastq_file):
    """
    Read sequences from a FASTQ file and store them in a dictionary.

    Parameters:
        fastq_file (str): The path to the FASTQ file.

    Returns:
        dict: A dictionary where keys are sequence IDs and values are sequences.
    """
    sequences = { }
    with open(fastq_file) as f:
        for title, seq, qual in FastqGeneralIterator(f):
            sequences[title.split()[2]] = seq
    return sequences


def parse_paf(file_path):
    """
    Parse a PAF file and convert it to a DataFrame.

    Parameters:
        file_path (str): The path to the PAF file.

    Returns:
        pd.DataFrame: The DataFrame containing the parsed PAF data.
    """
    df = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], header=None)
    df.columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand', "target_id", 'target_length',
                  'target_start', 'target_end']
    return df


def remove_self_loops(graph):
    """
    Remove self-loops from the graph.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The graph with self-loops removed.
    """
    graph.remove_edges_from(nx.selfloop_edges(graph))
    return graph


def remove_isolated_nodes(graph):
    """
    Remove isolated nodes from the graph.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The graph with isolated nodes removed.
    """
    isolated_nodes = list(nx.isolates(graph))
    graph.remove_nodes_from(isolated_nodes)
    return graph


def merge_nodes_with_single_edge(graph):
    """
    Merge nodes in the graph that have only a single edge connecting them.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The graph with nodes merged based on single edges.
    """
    nodes_to_merge = [node for node in graph.nodes() if graph.degree(node) == 1]
    while nodes_to_merge:
        node = nodes_to_merge.pop()
        neighbors = list(graph.neighbors(node))
        if len(neighbors) == 0:
            continue
        neighbor = neighbors[0]
        graph.nodes[neighbor]['sequence'] += graph.nodes[node]['sequence']
        graph = nx.contracted_nodes(graph, node, neighbor, self_loops=False)
        nodes_to_merge = [n for n in nodes_to_merge if n != neighbor]
    return graph


def remove_transitive_edges(graph):
    """
    Remove transitive edges from the graph.

    Parameters:
        graph (nx.DiGraph): The input directed graph.

    Returns:
        nx.DiGraph: The graph with transitive edges removed.
    """
    reduced_graph = nx.transitive_reduction(nx.to_directed(nx.DiGraph(graph)))
    transitive_edges = list(set(graph.edges()) - set(reduced_graph.edges()))
    graph.remove_edges_from(transitive_edges)
    return graph


def remove_single_connection_edges(graph):
    """
    Remove edges that connect nodes with a single connection each.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The graph with single connection edges removed.
    """
    single_connection_edges = [edge for edge in graph.edges() if
                               graph.degree(edge[0]) == 1 and graph.degree(edge[1]) == 1]
    graph.remove_edges_from(single_connection_edges)
    return graph


def collapse_linear_paths(graph):
    """
    Collapse linear paths in the graph into a single edge.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The graph with linear paths collapsed.
    """
    collapsed_graph = graph.copy()
    nodes_to_remove = []

    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))

        if len(neighbors) == 2:
            prev_node, next_node = neighbors

            if len(list(graph.neighbors(prev_node))) == 1 and len(list(graph.neighbors(next_node))) == 1:
                new_sequence = graph.nodes[prev_node]['sequence'] + graph.nodes[node]['sequence'] + \
                               graph.nodes[next_node]['sequence']

                collapsed_graph.add_edge(prev_node, next_node, sequence=new_sequence)
                nodes_to_remove.append(node)

    for node in nodes_to_remove:
        collapsed_graph.remove_node(node)

    return collapsed_graph


def remove_non_intersection_nodes(graph):
    """
    Remove non-intersection nodes from the graph.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The graph with non-intersection nodes removed.
    """
    non_intersection_nodes = [node for node in graph.nodes() if graph.in_degree(node) == 2]
    for node in non_intersection_nodes:
        neighbors = list(graph.neighbors(node))
        if len((neighbors)) >= 2:
            graph.remove_node(node)
            graph.add_edge(neighbors[0], neighbors[1], weight=1)
    return graph


def remove_cycles(graph):
    """
    Remove cycles from the graph.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The graph with cycles removed.
    """
    cycles = list(nx.simple_cycles(graph))
    for cycle in cycles:
        graph.remove_edges_from(cycle)
    return graph


def find_communities(graph):
    """
    Find communities in the graph using the Louvain method.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        list: A list of communities, where each community is represented as a list of nodes.
    """
    # Convert the networkx graph to an igraph graph for community detection
    ig_graph = ig.Graph.TupleList(nx.to_undirected(graph).edges(), directed=False)

    # Use the Louvain method to detect communities
    partition = ig_graph.community_leiden()

    # Convert the detected communities back to a list of node lists
    communities = []
    for com in set(partition.membership):
        members = [nodes for nodes in partition.graph.vs['name'] if
                   partition.membership[partition.graph.vs.find(name=nodes).index] == com]
        communities.append(members)

    return communities


def overlap_graph(sequences, overlaps):
    """
    Create an overlap graph from sequences and their overlaps.

    Parameters:
        sequences (dict): A dictionary where keys are sequence IDs and values are sequences.
        overlaps (pd.DataFrame): The DataFrame containing overlap information.

    Returns:
        nx.DiGraph: The constructed overlap graph.
    """
    graph = nx.DiGraph()
    for seq_id, seq in sequences.items():
        graph.add_node(seq_id, sequence=seq)
    for _, row in overlaps.iterrows():
        query_id = row['query_id']
        target_id = row['target_id']
        if query_id == target_id:
            continue
        query_start = row['query_start']
        query_end = row['query_end']
        target_start = row['target_start']
        target_end = row['target_end']
        if query_start < query_end and target_start < target_end:
            orientation = '++'
        elif query_start < query_end and target_start > target_end:
            orientation = '+-'
        elif query_start > query_end and target_start < target_end:
            orientation = '-+'
        else:
            orientation = '--'
        weight = min(query_end - query_start, target_end - target_start)
        graph.add_edge(query_id, target_id, weight=weight, orientation=orientation)
    return graph


def merge_communities(graph):
    """
    Merge communities in the graph.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.DiGraph: The graph with communities merged.
    """
    # Convert the directed graph to an undirected graph for partition computation
    undirected_graph = graph.to_undirected()

    # Compute the best partition using the Louvain method
    partition = community_louvain.best_partition(undirected_graph)

    # Create a new graph to represent the merged communities
    merged_graph = nx.DiGraph()

    # Add nodes to the merged graph, one for each community
    for community_id in set(partition.values()):
        merged_graph.add_node(community_id, sequence='')

    # Add the sequences of the original nodes to the corresponding community nodes
    for node, community_id in partition.items():
        merged_graph.nodes[community_id]['sequence'] += graph.nodes[node]['sequence']

    # Add edges between communities in the merged graph
    for edge in graph.edges():
        community1 = partition[edge[0]]
        community2 = partition[edge[1]]
        if community1 != community2:
            if merged_graph.has_edge(community1, community2):
                merged_graph[community1][community2]['weight'] += 1
            else:
                merged_graph.add_edge(community1, community2, weight=1)

    # Concatenate the nucleotides into a single sequence for each community
    for node in merged_graph:
        merged_graph.nodes[node]['sequence'] = ''.join(merged_graph.nodes[node]['sequence'])

    return merged_graph


def sparsify_graph(graph, sparsity):
    """
    Sparsify the graph by removing low-weight edges.

    Parameters:
        graph (nx.Graph): The input graph.
        sparsity (float): The percentage of edges to keep, represented as a value between 0 and 100.

    Returns:
        nx.Graph: The sparsified graph.
    """
    # Compute edge weights
    weights = np.array([data['weight'] for _, _, data in graph.edges(data=True)])

    # Compute threshold for edge selection based on sparsity percentage
    threshold = np.percentile(weights, sparsity)

    # Create a new sparse graph
    sparse_graph = nx.DiGraph()

    # Add nodes and edges to the sparse graph based on the threshold
    for u, v, data in graph.edges(data=True):
        if data['weight'] >= threshold:
            u_seq = graph.nodes[u]['sequence']
            v_seq = graph.nodes[v]['sequence']
            sparse_graph.add_node(u, sequence=u_seq)
            sparse_graph.add_node(v, sequence=v_seq)
            sparse_graph.add_edge(u, v, weight=data['weight'], orientation=data['orientation'])

    return sparse_graph


def simplify_graph(graph):
    """
    Simplify the graph by applying various graph simplification techniques.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        nx.Graph: The simplified graph.
    """
    # Sparsify the graph by removing low-weight edges
    graph = sparsify_graph(graph, sparsity=0.01)

    # Remove transitive edges from the graph
    graph = remove_transitive_edges(graph)

    # Remove cycles from the graph
    graph = remove_cycles(graph)

    # Remove edges that connect nodes with a single connection each
    graph = remove_single_connection_edges(graph)

    # Remove self-loops from the graph
    graph = remove_self_loops(graph)

    # Remove isolated nodes from the graph
    graph = remove_isolated_nodes(graph)

    # Merge communities in the graph (Optional: Uncomment the following line if desired)
    # graph = merge_communities(graph)

    # Collapse linear paths in the graph (Optional: Uncomment the following line if desired)
    # graph = collapse_linear_paths(graph)

    # Remove non-intersection nodes from the graph (Optional: Uncomment the following line if desired)
    # graph = remove_non_intersection_nodes(graph)

    return graph


def write_contigs_to_file(contigs, output_file):
    """
    Write contigs to a file.

    Parameters:
        contigs (list): A list of contigs, where each contig is represented as a list of nodes and sequences.
        output_file (str): The path to the output file.
    """
    with open(output_file, 'w') as f:
        count = 0
        for i, contig in enumerate(contigs):
            count += len(contig[0][1])
        print("Length of all contigs together:", count)


def draw_overlap_graph(graph):
    """
    Draw the overlap graph with labeled edges and edge weights.

    Parameters:
        graph: The input directed graph representing the overlap relationships.
    """
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True, node_color='lightblue', node_size=3000, font_size=12, font_weight='bold',
            arrowsize=20)
    edge_labels = { (u, v): d['weight'] for u, v, d in graph.edges(data=True) }
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, font_size=12)
    plt.show()


def create_contigs(graph):
    """
    Create contigs from the graph.

    Parameters:
        graph (nx.Graph): The input graph.

    Returns:
        list: A list of contigs, where each contig is represented as a list of nodes and sequences.
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


def merge_contigs(contigs, similarity_threshold=0.8):
    """
    Merge similar contigs into longer contigs.

    Parameters:
        contigs (list): A list of contigs, where each contig is represented as a list of nodes and sequences.
        similarity_threshold (float, optional): The minimum similarity required to merge contigs.
                                               Default is 0.8 (80% similarity).

    Returns:
        list: A list of merged contigs, where each contig is represented as a list of nodes and sequences.
    """
    merged_contigs = []
    visited = [False] * len(contigs)

    def merge_sequences(seq1, seq2):
        for i in range(len(seq1)):
            if seq1[i:] == seq2[:len(seq1) - i]:
                return seq1 + seq2[len(seq1) - i:]
        return None

    def find_similar_contig(contig, contigs):
        for i, other_contig in enumerate(contigs):
            if not visited[i]:
                merged_seq = merge_sequences(contig[0][1], other_contig[0][1])
                if merged_seq:
                    similarity = len(merged_seq) / (len(contig[0][1]) + len(other_contig[0][1]))
                    if similarity >= similarity_threshold:
                        return i, merged_seq
        return None, None

    for i, contig in enumerate(contigs):
        if not visited[i]:
            visited[i] = True
            current_contig = contig
            while True:
                similar_index, merged_seq = find_similar_contig(current_contig, contigs)
                if similar_index is not None:
                    visited[similar_index] = True
                    current_contig = [(current_contig[0][0], merged_seq)]
                else:
                    break
            merged_contigs.append(current_contig)

    return merged_contigs


def main():
    """
    Orchestrate the sequence assembly process (Not finished).
    Usage:
        1. Set the paths to the input and output files (fastq_file_path, paf_file_path, output_file_path).
        2. Uncomment the optional step to visualize the overlap graph if desired.
        3. Execute the script to initiate the sequence assembly process.
    """

    # Set the paths to the input and output files
    fastq_file_path = 'foo-reads.fq'
    paf_file_path = 'foo.paf'
    output_file_path = 'assembly.fasta'

    # Step 1: Read sequences from the FASTQ file
    sequences = read_fastq_file(fastq_file_path)

    # Step 2: Parse overlap information from the PAF file and create the overlap graph
    overlaps = parse_paf(paf_file_path)
    graph = overlap_graph(sequences, overlaps)

    # Step 3: Simplify the overlap graph
    simplified_graph = simplify_graph(graph)

    # Step 4: Create contigs from the simplified graph
    contigs = create_contigs(simplified_graph)

    # Step 5: Merge similar contigs into longer contigs
    merged_contigs = merge_contigs(contigs)

    # Step 6: Write the merged contigs to the output file
    write_contigs_to_file(merged_contigs, output_file_path)

    # Optionally, visualize the overlap graph (uncomment the following line if desired)
    draw_overlap_graph(simplified_graph)


if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    print("Execution time:", end_time - start_time, "seconds")
