#!/usr/bin/env python3

"""
Unit test script for the assembly.py
"""

__author__ = "Lisan Eisinga"
__version__ = "1.3.0"
__date__ = "31-05-2023"


# Imports
import os
import unittest
import pandas as pd
import networkx as nx
from Assembly import read_fastq_file, create_paf, \
    parse_paf, overlap_graph, remove_isolated_nodes, \
    dfs, get_consensus_sequences, write_to_file


class ReadFastqFileTest(unittest.TestCase):
    """Unit tests for the read_fastq_file function."""

    def setUp(self):
        """Set up the test environment."""
        self.valid_fastq_file = 'valid.fastq'
        self.invalid_fastq_file = 'invalid.fastq'

        # Create a valid FastQ file for testing
        with open(self.valid_fastq_file, 'w') as f:
            f.write('@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nJJJJ\n')

    def tearDown(self):
        """Clean up the test environment."""
        # Remove the temporary FastQ file
        os.remove(self.valid_fastq_file)

    def test_read_fastq_file_valid_file(self):
        """Test reading a valid FastQ file."""
        expected_sequences = {'seq1': 'ACGT', 'seq2': 'TGCA'}
        result = read_fastq_file(self.valid_fastq_file)
        self.assertEqual(result, expected_sequences)

    def test_read_fastq_file_invalid_file(self):
        """Test reading a non-existent FastQ file."""
        with self.assertRaises(FileNotFoundError):
            read_fastq_file(self.invalid_fastq_file)


class TestCreatePAF(unittest.TestCase):
    """Unit tests for the create_paf function."""

    def test_create_paf(self):
        """Test creating a PAF file from a FastQ file."""
        # Create a temporary FastQ file
        fastq_file = "unittest.fastq"
        with open(fastq_file, "w") as f:
            f.write("@read1\nATCG\n+\nAAAA\n")

        # Call the function to create the PAF file
        paf_file = create_paf(fastq_file)

        # Assert that the PAF file exists
        self.assertTrue(os.path.exists(paf_file))

        # Cleanup: Delete the temporary files
        os.remove(fastq_file)
        os.remove(paf_file)

    def test_create_paf_non_existing_file(self):
        """Test creating a PAF file from a non-existent FastQ file."""
        non_existing_file = "non_existing.fastq"
        # Assert that a FileNotFoundError is raised
        with self.assertRaises(FileNotFoundError):
            create_paf(non_existing_file)


class ParsePafTest(unittest.TestCase):
    """Unit tests for the parse_paf function."""

    def setUp(self):
        """Set up the test environment."""
        self.valid_paf_file = 'valid.paf'
        self.invalid_paf_file = 'invalid.paf'

        # Create a valid PAF file for testing
        with open(self.valid_paf_file, 'w') as f:
            f.write('query1\t100\t0\t50\t+\ttarget1\t200\t50\t100\nquery2\t150\t0\t100\t-\ttarget2\t300\t0\t150\n')

    def tearDown(self):
        """Clean up the test environment."""
        # Remove the temporary PAF file
        os.remove(self.valid_paf_file)

    def test_parse_paf_valid_file(self):
        """Test parsing a valid PAF file."""
        expected_columns = ['query_id', 'query_length', 'query_start', 'query_end', 'strand', 'target_id',
                            'target_length', 'target_start', 'target_end']
        expected_data = [
            ['query1', 100, 0, 50, '+', 'target1', 200, 50, 100],
            ['query2', 150, 0, 100, '-', 'target2', 300, 0, 150]
        ]
        expected_df = pd.DataFrame(expected_data, columns=expected_columns)

        result = parse_paf(self.valid_paf_file)
        pd.testing.assert_frame_equal(result, expected_df)

    def test_parse_paf_invalid_file(self):
        """Test parsing a non-existent PAF file."""
        with self.assertRaises(FileNotFoundError):
            parse_paf(self.invalid_paf_file)


class OverlapGraphTest(unittest.TestCase):
    """Unit tests for the overlap_graph function."""

    def test_overlap_graph_valid_input(self):
        """Test creating an overlap graph with valid input."""
        sequences = {
            'seq1': 'ACGT',
            'seq2': 'TGCA',
            'seq3': 'AAAA'
        }

        overlaps_data = {
            'query_id': ['seq1', 'seq2', 'seq2'],
            'target_id': ['seq2', 'seq1', 'seq3'],
            'query_length': [4, 4, 4],
            'target_length': [4, 4, 4],
            'query_start': [0, 0, 0],
            'query_end': [3, 3, 3],
            'target_start': [0, 0, 0],
            'target_end': [3, 3, 3],
            'strand': ['+', '-', '+']
        }
        overlaps_df = pd.DataFrame(overlaps_data)

        expected_graph = nx.MultiDiGraph()
        expected_graph.add_node('seq1', sequence='ACGT')
        expected_graph.add_node('seq2', sequence='TGCA')
        expected_graph.add_node('seq3', sequence='AAAA')
        expected_graph.add_edge('seq1', 'seq2', query_start=0, query_end=3, target_start=0,
                                target_end=3, strand='+', query_length=4, target_length=4,
                                overlap_len=4)
        expected_graph.add_edge('seq2', 'seq1', query_start=0, query_end=3, target_start=0,
                                target_end=3, strand='-', query_length=4, target_length=4,
                                overlap_len=4)
        expected_graph.add_edge('seq2', 'seq3', query_start=0, query_end=3, target_start=0,
                                target_end=3, strand='+', query_length=4, target_length=4,
                                overlap_len=4)

        result_graph = overlap_graph(sequences, overlaps_df)
        self.assertEqual(result_graph.nodes, expected_graph.nodes)
        self.assertEqual(result_graph.edges, expected_graph.edges)

    def test_overlap_graph_missing_columns(self):
        """Test handling missing columns in the overlaps DataFrame."""
        sequences = {
            'seq1': 'ACGT',
            'seq2': 'TGCA'
        }

        overlaps_data = {
            'query_id': ['seq1', 'seq2'],
            'target_id': ['seq2', 'seq1'],
            'query_length': [4, 4]
            # Missing columns: 'target_length', 'query_start', 'query_end', 'target_start', 'target_end', 'strand'
        }
        overlaps_df = pd.DataFrame(overlaps_data)

        with self.assertRaises(KeyError):
            overlap_graph(sequences, overlaps_df)

    def test_overlap_graph_invalid_sequence_id(self):
        """Test handling an invalid sequence ID in the overlaps DataFrame."""
        sequences = {
            'seq1': 'ACGT',
            'seq2': 'TGCA'
        }

        overlaps_data = {
            'query_id': ['seq1', 'seq2', 'seq3'],  # seq3 is not present in the sequences dictionary
            'target_id': ['seq2', 'seq1', 'seq2'],
            'query_length': [4, 4, 4],
            'target_length': [4, 4, 4],
            'query_start': [0, 0, 0],
            'query_end': [3, 3, 3],
            'target_start': [0, 0, 0],
            'target_end': [3, 3, 3],
            'strand': ['+', '-', '+']
        }
        overlaps_df = pd.DataFrame(overlaps_data)

        with self.assertRaises(ValueError):
            overlap_graph(sequences, overlaps_df)


class RemoveIsolatedNodesTest(unittest.TestCase):
    """Unit tests for the remove_isolated_nodes function."""

    def test_remove_isolated_nodes_valid_graph(self):
        """Test removing isolated nodes from a valid graph."""
        graph = nx.MultiDiGraph()
        graph.add_node('A')
        graph.add_node('B')
        graph.add_edge('A', 'B')

        expected_graph = nx.MultiDiGraph()
        expected_graph.add_node('A')
        expected_graph.add_node('B')
        expected_graph.add_edge('A', 'B')

        result_graph = remove_isolated_nodes(graph)
        self.assertEqual(result_graph.nodes, expected_graph.nodes)
        self.assertEqual(result_graph.edges, expected_graph.edges)

    def test_remove_isolated_nodes_empty_graph(self):
        """Test removing isolated nodes from an empty graph."""
        graph = nx.MultiDiGraph()

        expected_graph = nx.MultiDiGraph()

        result_graph = remove_isolated_nodes(graph)
        self.assertEqual(result_graph.nodes, expected_graph.nodes)
        self.assertEqual(result_graph.edges, expected_graph.edges)

    def test_remove_isolated_nodes_invalid_graph_type(self):
        """Test handling an invalid graph type."""
        graph = nx.Graph()  # Using Graph instead of MultiDiGraph

        with self.assertRaises(TypeError):
            remove_isolated_nodes(graph)


class DFSTest(unittest.TestCase):
    """Unit tests for the dfs function."""

    def test_dfs_valid_graph(self):
        """Test finding connected components (contigs) in a valid graph."""
        graph = nx.MultiDiGraph()
        graph.add_node('A')
        graph.add_node('B')
        graph.add_node('C')
        graph.add_node('D')
        graph.add_edge('A', 'B', target_start=0, query_start=0)
        graph.add_edge('B', 'C', target_start=1, query_start=1)
        graph.add_edge('A', 'D', target_start=2, query_start=2)

        expected_contigs = [['A', 'B', 'C', 'D']]

        result_contigs = dfs(graph)
        self.assertEqual(result_contigs, expected_contigs)

    def test_dfs_multiple_contigs(self):
        """Test finding multiple contigs in the graph."""
        graph = nx.MultiDiGraph()
        graph.add_node('A')
        graph.add_node('B')
        graph.add_node('C')
        graph.add_node('D')
        graph.add_edge('A', 'B', target_start=0, query_start=0)
        graph.add_edge('C', 'D', target_start=1, query_start=0)

        expected_contigs = [['A', 'B'], ['C', 'D']]

        result_contigs = dfs(graph)
        self.assertEqual(result_contigs, expected_contigs)

    def test_dfs_empty_graph(self):
        """Test finding contigs in an empty graph."""
        graph = nx.MultiDiGraph()

        expected_contigs = []

        result_contigs = dfs(graph)
        self.assertEqual(result_contigs, expected_contigs)

    def test_dfs_invalid_graph_type(self):
        """Test handling an invalid graph type."""
        graph = nx.Graph()  # Using Graph instead of MultiDiGraph

        with self.assertRaises(TypeError):
            dfs(graph)


class ConsensusSequenceTest(unittest.TestCase):
    """Unit tests for the get_consensus_sequences function."""

    def test_get_consensus_sequences(self):
        """Test generating consensus sequences from a graph and contigs."""
        # Create a test graph
        graph = nx.MultiDiGraph()

        graph.add_node("A", sequence="ATCGTAGAT")
        graph.add_node("B", sequence="ACCGTAG")
        graph.add_node("C", sequence="TAATCGT")

        graph.add_edge("A", "C", query_start=2, query_end=8,
                       target_start=0, target_end=3)
        graph.add_edge("B", "C", query_start=4, query_end=5,
                       target_start=0, target_end=5)

        # Define contigs
        contigs = [["A", "C", "B"]]

        # Define expected consensus sequences
        expected_consensus_sequences = ["TAATATCGTAGAT"]

        # Generate consensus sequences
        result_consensus_sequences = get_consensus_sequences(graph, contigs)

        # Compare the expected and result consensus sequences
        self.assertEqual(expected_consensus_sequences, result_consensus_sequences)

    def test_get_consensus_sequences_empty_contigs(self):
        """Test handling an empty contigs list."""
        # Create a test graph
        graph = nx.MultiDiGraph()

        # Define an empty contigs list
        contigs = []

        # Test for expected ValueError
        with self.assertRaises(ValueError):
            get_consensus_sequences(graph, contigs)

    def test_get_consensus_sequences_invalid_graph_type(self):
        """Test handling an invalid graph type."""
        # Create an invalid graph type (using Graph instead of MultiDiGraph)
        graph = nx.Graph()
        graph.add_node('A', sequence='ACGT')

        # Define a contig
        contigs = [['A']]

        # Test for expected TypeError
        with self.assertRaises(TypeError):
            get_consensus_sequences(graph, contigs)


class WriteToFileTest(unittest.TestCase):
    """Unit tests for the write_to_file function."""

    def test_write_to_file(self):
        """Test writing contigs to a file."""
        # Test data
        contigs = ['ACGT', 'GCTA', 'TACG']
        filename = 'unittest.fasta'

        # Call the function
        write_to_file(filename, contigs)

        # Check if the file was created
        self.assertTrue(os.path.exists(filename))

        # Read the file content
        with open(filename, 'r') as file:
            file_content = file.read()

        # Check if the file content matches the expected output
        expected_output = ">contig_1\nACGT\n>contig_2\nGCTA\n>contig_3\nTACG\n"
        self.assertEqual(file_content, expected_output)

        # Clean up - delete the test file
        os.remove(filename)


if __name__ == '__main__':
    unittest.main()
