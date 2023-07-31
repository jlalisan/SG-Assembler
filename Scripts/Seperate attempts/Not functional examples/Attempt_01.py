#!/usr/bin/env python3
"""
Unfinished Script - Not Functional

This script is intended to process sequence data, identify overlaps between sequences,
and construct contigs.

Note: This script is unfinished and not fully functional.
It requires additional development and testing to work as intended.
"""

import time
import networkx as nx
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity


def read_fastq_file(fastq_file):
    """
    Read a FastQ file and extract sequences.

    Parameters:
        fastq_file (str): The path to the FastQ file.

    Returns:
        dict: A dictionary where the sequence titles are keys and the sequences are values.
    """
    sequences = {}
    with open(fastq_file) as f:
        for title, seq, qual in FastqGeneralIterator(f):
            sequences[title.split()[2]] = seq
    return sequences


def parse_paf(file_path):
    """
    Parse a PAF file and create a DataFrame containing relevant information.

    Parameters:
        file_path (str): The path to the PAF file.

    Returns:
        pandas.DataFrame: A DataFrame containing columns: 'query_id', 'query_length',
        'query_start', 'query_end', 'strand', 'target_id', 'target_length', 'target_start',
        'target_end'.
    """
    df = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], header=None)
    df.columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand', "target_id",
                  'target_length', 'target_start', 'target_end']
    return df


def overlap_graph(sequences, overlaps):
    """
    Create an overlap graph using sequences and overlaps.

    Parameters:
        sequences (dict): A dictionary where the sequence titles are keys
        and the sequences are values.
        overlaps (pandas.DataFrame): DataFrame containing information about overlaps
        between sequences.

    Returns:
        networkx.MultiDiGraph: A directed multi-graph representing overlaps between sequences.
    """
    graph = nx.MultiDiGraph()
    for seq_id, seq in sequences.items():
        graph.add_node(seq_id, sequence=seq)
    for _, row in overlaps.iterrows():
        query_id = row['query_id']
        target_id = row['target_id']
        if query_id == target_id:
            continue
        weight = row['query_end'] - row['query_start']
        graph.add_edge(query_id, target_id, weight=weight)
    return graph


# ... Other functions go here ...

def identify_transitive_edges(graph):
    """
    Identify and remove transitive edges from a graph using transitive reduction.

    Parameters:
        graph (networkx.MultiDiGraph): The original graph.

    Returns:
        networkx.MultiDiGraph: The reduced graph after removing transitive edges.
    """
    reduced_graph = nx.transitive_reduction(nx.MultiDiGraph(graph))
    return reduced_graph


def remove_transitive_edges(graph):
    """
    Remove transitive edges from the graph.

    Parameters:
        graph (networkx.MultiDiGraph): The original graph.

    Returns:
        networkx.MultiDiGraph: The graph with transitive edges removed.
    """
    reduced_graph = identify_transitive_edges(graph)
    transitive_edges = list(set(graph.edges()) - set(reduced_graph.edges()))
    graph.remove_edges_from(transitive_edges)
    return graph


def remove_small_bubbles(graph, threshold=2):
    """
    Remove small bubbles from the graph.

    A small bubble is a node with one predecessor and one successor,
    and both connected nodes also have only one
    predecessor and one successor.
    The sequence of the node in the bubble and its connected nodes is concatenated,
    and the bubble node is removed from the graph.

    Parameters:
        graph (networkx.MultiDiGraph): The original graph.
        threshold (int, optional): The maximum size of the bubble to be removed. Defaults to 2.

    Returns:
        networkx.MultiDiGraph: The graph with small bubbles removed.
    """
    for node in list(graph.nodes()):
        preds = list(graph.predecessors(node))
        succs = list(graph.successors(node))
        if len(preds) == 1 and len(succs) == 1:
            pred = str(preds)
            succ = str(succs)
            if graph.has_node(pred) and graph.has_node(succ) and len(list(graph.successors(pred))) == 1 and len(
                    list(graph.predecessors(succ))) == 1:
                pred_seq = graph.nodes[pred]['sequence']
                node_seq = graph.nodes[node]['sequence']
                succ_seq = graph.nodes[succ]['sequence']
                if len(node_seq) <= threshold and len(pred_seq) + len(node_seq) + len(succ_seq) <= threshold + 2:
                    graph.remove_node(node)
                    graph.add_edge(pred, succ, sequence=pred_seq + node_seq + succ_seq)

    return graph


def simplify_graph(unitig_graph):
    """
    Simplify the graph by removing small bubbles and construct unitigs.

    Parameters:
        unitig_graph (networkx.MultiDiGraph): The original graph.

    Returns:
        networkx.MultiDiGraph: The simplified graph containing unitigs.
    """
    unitig_graph = remove_small_bubbles(unitig_graph)

    unitigs = []
    for node in list(unitig_graph.nodes()):
        preds = list(unitig_graph.predecessors(node))
        succs = list(unitig_graph.successors(node))
        if len(preds) == 1 and len(succs) == 1:
            unitig = [node]
            while len(succs) == 1:
                node = "".join(succs)
                unitig.append(node)
                succs = list(unitig_graph.successors(node))
            unitigs.append(unitig)

    simplified_graph = nx.MultiDiGraph()
    for unitig in unitigs:
        seq = ''
        for i in range(len(unitig)):
            node = unitig[i]
            seq += unitig_graph.nodes[node]['sequence']
            if i < len(unitig) - 1:
                next_node = unitig[i + 1]
                edges = unitig_graph.get_edge_data(node, next_node)
                for _, edge in edges.items():
                    simplified_graph.add_edge(node, next_node, **edge)
        for node in unitig:
            simplified_graph.nodes[node]['sequence'] = seq

    return simplified_graph


def extract_contigs(graph):
    """
    Extract contigs from the simplified graph.

    Parameters:
        graph (networkx.MultiDiGraph): The simplified graph containing unitigs.

    Returns:
        list: A list of contigs (sequences).
    """
    contigs = []
    graph = remove_transitive_edges(graph)
    graph = graph.to_directed()
    while graph:
        path = nx.algorithms.dag.dag_longest_path(graph, weight='weight')
        sequence = ''.join([graph.nodes[node]['sequence'] for node in path])
        contigs.append(sequence)
        graph.remove_nodes_from(path)
        graph = remove_transitive_edges(graph)
    return contigs


def calculate_contig_similarity(contigs):
    """
    Calculate the pairwise similarity between contigs using TF-IDF vectorization and cosine similarity.

    Parameters:
        contigs (list): A list of contigs (sequences).

    Returns:
        numpy.ndarray: A similarity matrix representing the pairwise similarity between contigs.
    """
    contig_strings = [''.join(contig) for contig in contigs]
    vectorizer = TfidfVectorizer(analyzer='char', ngram_range=(4, 4))
    contig_vectors = vectorizer.fit_transform(contig_strings)
    similarity_matrix = cosine_similarity(contig_vectors)
    return similarity_matrix


def merge_contigs_based_on_similarity(contigs, similarity_matrix, threshold=0.8):
    """
    Merge similar contigs based on their pairwise similarity scores.

    Parameters:
        contigs (list): A list of contigs (sequences).
        similarity_matrix (numpy.ndarray): A similarity matrix representing the pairwise similarity between contigs.
        threshold (float, optional): The similarity threshold for merging contigs. Defaults to 0.8.

    Returns:
        list: A list of merged contigs.
    """
    merged_contigs = []
    merged_indices = set()

    for i in range(len(contigs)):
        if i in merged_indices:
            continue

        similar_contigs = [contigs[i]]
        for j in range(i + 1, len(contigs)):
            if similarity_matrix[i][j] >= threshold:
                similar_contigs.append(contigs[j])
                merged_indices.add(j)

        merged_contig = ''.join(similar_contigs)
        merged_contigs.append(merged_contig)

    return merged_contigs


def main():
    """
    Main function to process sequence data, identify overlaps, and construct contigs.

    Reads input files, constructs the overlap graph, simplifies the graph, extracts contigs,
    calculates contig similarity, and finally merges similar contigs. Prints various statistics
    related to the processed data.
    """
    fastq_file_path = 'foo-reads.fq'
    paf_file_path = 'foo.paf'
    output_file_path = 'output.fasta'

    start = time.time()
    sequences = read_fastq_file(fastq_file_path)
    overlaps = parse_paf(paf_file_path)
    graph = overlap_graph(sequences, overlaps)
    removed_graph = remove_transitive_edges(graph)
    bubbles = remove_small_bubbles(removed_graph)
    simplified_graph = simplify_graph(removed_graph)
    contigs = extract_contigs(simplified_graph)
    similarity_matrix = calculate_contig_similarity(contigs)
    merged_contigs = merge_contigs_based_on_similarity(contigs, similarity_matrix)

    print(simplified_graph)

    total_length = 0
    for contig in sequences.values():
        total_length += len(contig)
    print("Total sequence length before processing:", total_length)

    total = 0
    for node in simplified_graph.nodes(data=True):
        total += len(node[1]['sequence'])
    print("Total sequence length after processing:", total)

    contig_length = 0
    for contig in merged_contigs:
        contig_length += len(contig)
    print("Total contig length after merging:", contig_length)
    print("Total number of merged contigs:", len(merged_contigs))

    end = time.time()
    print("Total elapsed time:", end - start)

if __name__ == '__main__':
    main()
