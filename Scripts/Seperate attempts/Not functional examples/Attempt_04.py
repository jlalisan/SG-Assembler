#!/usr/bin/env python3
"""
Unfinished and Non-Functional Assembly Graph Construction Script

This script is intended to construct an assembly graph from DNA sequences and overlap data. However, it is still
incomplete and not fully functional.

It defines functions to read DNA sequences from a FastQ file, parse overlap information from a PAF file, and
construct an overlap graph using the A-Bruijn graph approach. The graph is then simplified by removing unique edges,
and a depth-first traversal of the graph is performed to find visited nodes.

Please note that this is an unfinished script and should not be used for any real-world applications without proper
modification and completion.
"""
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def read_fastq_file(fastq_file):
    """Read sequences from a FastQ file and store them in a dictionary.

    Args:
        fastq_file (str): Path to the FastQ file.

    Returns:
        dict: A dictionary where keys are read IDs and values are sequences.
    """
    sequence = {}
    with open(fastq_file) as f:
        for title, seq, qual in FastqGeneralIterator(f):
            sequence[title.split()[2]] = seq
    return sequence


def parse_paf(file_path):
    """Parse PAF (Pairwise mApping Format) file and convert it into a pandas DataFrame.

    Args:
        file_path (str): Path to the PAF file.

    Returns:
        pd.DataFrame: A DataFrame containing the PAF data with labeled columns.
    """
    df = pd.read_csv(file_path, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], header=None)
    df.columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand', "target_id",
                  'target_length', 'target_start', 'target_end']
    return df


class Graph:
    def __init__(self):
        """Initialize an assembly graph."""
        self.nodes = {}
        self.edges = {}

    def add_node(self, node_id, sequence):
        """Add a node to the assembly graph.

        Args:
            node_id (str): ID of the node.
            sequence (str): DNA sequence of the node.
        """
        if str(node_id) not in self.nodes:
            self.nodes[str(node_id)] = sequence

    def add_edge(self, source, target, overlap_data):
        """Add an edge to the assembly graph.

        Args:
            source (str): ID of the source node.
            target (str): ID of the target node.
            overlap_data (dict): Data containing overlap information between the nodes.
        """
        if source not in self.edges:
            self.edges[source] = {}
        if target not in self.edges:
            self.edges[target] = {}
        self.edges[source][target] = overlap_data

    def simplify_graph(self):
        """Simplify the assembly graph by removing unique edges."""
        # Classify all edges in the graph into unique and repeat edges
        unique_edges = set()
        repeat_edges = set()
        for source in self.edges:
            for target in self.edges[source]:
                overlap_data = self.edges[source][target]
                if source == self.nodes[target]:
                    unique_edges.add((source, target))
                else:
                    repeat_edges.add((source, target))

        # Remove all unique edges from the assembly graph
        for edge in unique_edges:
            del self.edges[edge[0]][edge[1]]

    def find_starting_node(self):
        """Find the starting node for depth-first traversal.

        Returns:
            str: ID of the starting node.
        """
        return max(self.edges, key=lambda node: len(self.edges[node]))

    def dfs_traversal(self):
        """Perform depth-first traversal on the assembly graph.

        Returns:
            set: Set of visited node IDs during the traversal.
        """
        visited = set()
        starting_node = self.find_starting_node()
        self._dfs_helper(starting_node, visited)
        return visited

    def _dfs_helper(self, node, visited):
        """Helper function for depth-first traversal.

        Args:
            node (str): ID of the current node.
            visited (set): Set to track visited node IDs.
        """
        if node not in visited:
            visited.add(node)
            for neighbor in self.edges[node]:
                self._dfs_helper(neighbor, visited)


def overlap_graph(sequences, overlaps):
    """Create an overlap graph using sequences and overlaps data.

    Args:
        sequences (dict): A dictionary with read IDs as keys and DNA sequences as values.
        overlaps (pd.DataFrame): A DataFrame containing overlap data between reads.

    Returns:
        Graph: An instance of the Graph class representing the overlap graph.
    """
    mygraph = Graph()
    # Sort the input sequences before processing
    sorted_sequences = sorted(sequences.items(), key=lambda x: x[0])
    # Add nodes to the graph
    for readId, seq in sorted_sequences:
        mygraph.add_node(f'{readId}', seq)  # add the 'read=' prefix to the node ID
    # Add edges to the graph using the A-Bruijn graph approach
    for _, row in overlaps.iterrows():
        mygraph.add_edge(f'{row["query_id"]}', f'{row["target_id"]}', {
            'query_start': row['query_start'],
            'query_end': row['query_end'],
            'target_start': row['target_start'],
            'target_end': row['target_end'],
            'strand': row['strand'],
            'overlap_length': row['query_end'] - row['query_start'] + 1
        })

    return mygraph

def main():
    """
    Main function to orchestrate the assembly graph construction process.

    1. Reads DNA sequences from a FastQ file.
    2. Parses overlap information from a PAF file.
    3. Constructs the overlap graph using A-Bruijn graph approach.
    4. Simplifies the graph by removing unique edges.
    5. Performs depth-first traversal on the graph to find visited nodes.

    Please ensure that the paths to the input FastQ and PAF files are correctly specified.
    """
    fastq_file_path = 'foo-reads.fq'
    paf_file_path = 'foo.paf'

    # Read DNA sequences from FastQ file
    sequences = read_fastq_file(fastq_file_path)

    # Parse overlap information from PAF file
    overlaps = parse_paf(paf_file_path)

    # Construct the overlap graph
    mygraph = overlap_graph(sequences, overlaps)

    # Simplify the graph
    mygraph.simplify_graph()

    # Perform depth-first traversal on the graph
    visited_nodes = mygraph.dfs_traversal()

    print("DFS traversal:", visited_nodes)


if __name__ == '__main__':
    main()
