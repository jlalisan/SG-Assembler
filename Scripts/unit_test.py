#!/usr/bin/env python3

"""
Unit test script for the assembly.py.
Used to make sure the output is correct if the code is changed.
"""

__author__ = "Lisan Eisinga"
__version__ = "1.5.0"
__date__ = "11-06-2023"

# Imports
import os
import unittest
import pandas as pd
import networkx as nx
from assembly import read_fastq_file, create_paf, \
    parse_paf, overlap_graph, remove_isolated_nodes, \
    dfs, generate_sequence, write_to_file


class ReadFastqFileTest(unittest.TestCase):
    """Unit tests for the read_fastq_file function."""

    def setUp(self):
        """Set up the test environment."""
        self.valid_fastq_file = "valid.fastq"
        self.invalid_fastq_file = "invalid.fastq"

        # Create a valid FastQ file for testing
        with open(self.valid_fastq_file, "w", encoding="utf-8") as file_handle:
            file_handle.write("@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nJJJJ\n")

    def tearDown(self):
        """Clean up the test environment."""
        # Remove the temporary FastQ file
        os.remove(self.valid_fastq_file)

    def test_read_fastq_file_valid_file(self):
        """Test reading a valid FastQ file."""
        expected_sequences = {"seq1": "ACGT", "seq2": "TGCA"}
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
        with open(fastq_file, "w", encoding="utf-8") as file_handle:
            file_handle.write("@read1\nATCG\n+\nAAAA\n")

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
        self.valid_paf_file = "valid.paf"
        self.invalid_paf_file = "invalid.paf"

        # Create a valid PAF file for testing
        with open(self.valid_paf_file, "w", encoding="utf-8") as file_handle:
            file_handle.write("query1\t100\t0\t50\t+\ttarget1\t200\t50\t100\t300\t43\t0\n"
                              "query2\t150\t0\t100\t-\ttarget2\t300\t0\t150\t400\t31\t0\n")

    def tearDown(self):
        """Clean up the test environment."""
        # Remove the temporary PAF file
        os.remove(self.valid_paf_file)

    def test_parse_paf_valid_file(self):
        """Test parsing a valid PAF file."""
        expected_columns = ["query_id", "query_length", "query_start", "query_end", "strand",
                            "target_id", "target_length", "target_start", "target_end",
                            "alignment_block_length", "residue_matches", "mapping_quality"]
        expected_data = [
            ["query1", 100, 0, 50, "+", "target1", 200, 50, 100, 300, 43, 0],
            ["query2", 150, 0, 100, "-", "target2", 300, 0, 150, 400, 31, 0]
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
        # Define the sequences and overlaps data for the test case
        sequences = {
            "seq1": "ACGT",
            "seq2": "TGCA",
            "seq3": "AAAA"
        }

        overlaps_data = {
            "query_id": ["seq1", "seq2", "seq2"],
            "target_id": ["seq2", "seq1", "seq3"],
            "query_length": [4, 4, 4],
            "target_length": [4, 4, 4],
            "query_start": [0, 0, 0],
            "query_end": [3, 3, 3],
            "target_start": [0, 0, 0],
            "target_end": [3, 3, 3],
            "strand": ["+", "-", "+"],
            "alignment_block_length": [4, 4, 4],
            "residue_matches": [4, 4, 4],
            "mapping_quality": [0, 0, 0]
        }
        overlaps_df = pd.DataFrame(overlaps_data)

        # Define the expected overlap graph
        expected_graph = nx.MultiDiGraph()
        expected_graph.add_node("seq1", sequence="ACGT")
        expected_graph.add_node("seq2", sequence="TGCA")
        expected_graph.add_node("seq3", sequence="AAAA")
        expected_graph.add_edge("seq1", "seq2", query_start=0, query_end=3, target_start=0,
                                target_end=3, strand="+", query_length=4, target_length=4,
                                overlap_len=4, alignment_block_length=4, residue_matches=4,
                                mapping_quality=0)
        expected_graph.add_edge("seq2", "seq1", query_start=0, query_end=3, target_start=0,
                                target_end=3, strand="-", query_length=4, target_length=4,
                                overlap_len=4, alignment_block_length=4, residue_matches=4,
                                mapping_quality=0)
        expected_graph.add_edge("seq2", "seq3", query_start=0, query_end=3, target_start=0,
                                target_end=3, strand="+", query_length=4, target_length=4,
                                overlap_len=4, alignment_block_length=4, residue_matches=4,
                                mapping_quality=0)

        # Create the overlap graph and verify the result
        result_graph = overlap_graph(sequences, overlaps_df)
        self.assertEqual(result_graph.nodes, expected_graph.nodes)
        self.assertEqual(result_graph.edges, expected_graph.edges)

    def test_overlap_graph_missing_columns(self):
        """Test handling missing columns in the overlaps DataFrame."""
        # Define the sequences and overlaps data for the test case
        sequences = {
            "seq1": "ACGT",
            "seq2": "TGCA"
        }

        overlaps_data = {
            "query_id": ["seq1", "seq2"],
            "target_id": ["seq2", "seq1"],
            "query_length": [4, 4],
            "alignment_block_length": [4, 4],
            "residue_matches": [4, 4],
            "mapping_quality": [0, 0]
        }
        overlaps_df = pd.DataFrame(overlaps_data)

        # Verify that a KeyError is raised when handling missing columns
        with self.assertRaises(KeyError):
            overlap_graph(sequences, overlaps_df)

    def test_overlap_graph_invalid_sequence_id(self):
        """Test handling an invalid sequence ID in the overlaps DataFrame."""
        # Define the sequences and overlaps data for the test case
        sequences = {
            "seq1": "ACGT",
            "seq2": "TGCA"
        }

        overlaps_data = {
            "query_id": ["seq1", "seq2", "seq3"],
            "target_id": ["seq2", "seq1", "seq2"],
            "query_length": [4, 4, 4],
            "target_length": [4, 4, 4],
            "query_start": [0, 0, 0],
            "query_end": [3, 3, 3],
            "target_start": [0, 0, 0],
            "target_end": [3, 3, 3],
            "strand": ["+", "-", "+"],
            "alignment_block_length": [4, 4, 4],
            "residue_matches": [4, 4, 4],
            "mapping_quality": [0, 0, 0]
        }
        overlaps_df = pd.DataFrame(overlaps_data)

        # Verify that a ValueError is raised when handling an invalid sequence ID
        with self.assertRaises(ValueError):
            overlap_graph(sequences, overlaps_df)


class RemoveIsolatedNodesTest(unittest.TestCase):
    """Unit tests for the remove_isolated_nodes function."""

    def test_remove_isolated_nodes_valid_graph(self):
        """Test removing isolated nodes from a valid graph."""
        # Create a valid graph with nodes and an edge
        graph = nx.MultiDiGraph()
        graph.add_node("A")
        graph.add_node("B")
        graph.add_edge("A", "B")

        # Define the expected graph after removing isolated nodes
        expected_graph = nx.MultiDiGraph()
        expected_graph.add_node("A")
        expected_graph.add_node("B")
        expected_graph.add_edge("A", "B")

        # Perform the removal of isolated nodes and verify the result
        result_graph = remove_isolated_nodes(graph)
        self.assertEqual(result_graph.nodes, expected_graph.nodes)
        self.assertEqual(result_graph.edges, expected_graph.edges)

    def test_remove_isolated_nodes_empty_graph(self):
        """Test removing isolated nodes from an empty graph."""
        # Create an empty graph
        graph = nx.MultiDiGraph()

        # Define the expected graph (empty graph)
        expected_graph = nx.MultiDiGraph()

        # Perform the removal of isolated nodes and verify the result
        result_graph = remove_isolated_nodes(graph)
        self.assertEqual(result_graph.nodes, expected_graph.nodes)
        self.assertEqual(result_graph.edges, expected_graph.edges)

    def test_remove_isolated_nodes_invalid_graph_type(self):
        """Test handling an invalid graph type."""
        # Create an invalid graph type (Graph instead of MultiDiGraph)
        graph = nx.Graph()

        # Verify that a TypeError is raised when removing isolated nodes from an invalid graph type
        with self.assertRaises(TypeError):
            remove_isolated_nodes(graph)


class DFSTest(unittest.TestCase):
    """Unit tests for the dfs function."""

    def test_dfs_valid_graph(self):
        """Test finding connected components (contigs) in a valid graph."""
        # Create a valid graph with nodes and edges
        # Each edge contains alignment information for sorting
        graph = nx.MultiDiGraph()
        graph.add_node("A")
        graph.add_node("B")
        graph.add_node("C")
        graph.add_node("D")
        graph.add_edge("A", "B", target_start=0, query_start=0, query_end=2, target_end=2)
        graph.add_edge("B", "C", target_start=1, query_start=1, query_end=3, target_end=3)
        graph.add_edge("A", "D", target_start=2, query_start=2, query_end=4, target_end=4)

        # Define the expected contigs in the graph
        expected_contigs = [["A", "B", "C", "D"]]

        # Perform the depth-first search (DFS) and verify the result
        result_contigs = dfs(graph)
        self.assertEqual(expected_contigs, result_contigs)

    def test_dfs_multiple_contigs(self):
        """Test finding multiple contigs in the graph."""
        # Create a graph with multiple contigs
        graph = nx.MultiDiGraph()
        graph.add_node("A")
        graph.add_node("B")
        graph.add_node("C")
        graph.add_node("D")
        graph.add_edge("A", "B", target_start=0, query_start=0, query_end=3, target_end=3)
        graph.add_edge("C", "D", target_start=1, query_start=0, query_end=4, target_end=3)

        # Define the expected contigs in the graph
        expected_contigs = [["A", "B"], ["C", "D"]]

        # Perform the depth-first search (DFS) and verify the result
        result_contigs = dfs(graph)
        self.assertEqual(result_contigs, expected_contigs)

    def test_dfs_empty_graph(self):
        """Test finding contigs in an empty graph."""
        # Create an empty graph
        graph = nx.MultiDiGraph()

        # Define the expected contigs (empty list)
        expected_contigs = []

        # Perform the depth-first search (DFS) and verify the result
        result_contigs = dfs(graph)
        self.assertEqual(result_contigs, expected_contigs)

    def test_dfs_invalid_graph_type(self):
        """Test handling an invalid graph type."""
        # Create an invalid graph type (Graph instead of MultiDiGraph)
        graph = nx.Graph()

        # Verify that a TypeError is raised when performing DFS on an invalid graph type
        with self.assertRaises(TypeError):
            dfs(graph)


class GenerateSequenceTest(unittest.TestCase):
    """ Unit test for the generation of sequences. """

    def test_generate_sequence(self):
        """ Tests to see it does not use the same node twice."""
        # Create a graph
        graph = nx.MultiDiGraph()
        graph.add_node("A", sequence="AT")
        graph.add_node("B", sequence="TC")
        graph.add_node("C", sequence="CG")
        graph.add_edge("A", "B", alignment_block_length=2)
        graph.add_edge("B", "C", alignment_block_length=2)

        # Create draft contigs
        draft_contigs = [["A", "B", "C"], ["A", "B", "B", "C"]]

        # Call the generate_sequence function
        accurate_contigs = generate_sequence(graph, draft_contigs)

        # Check the output
        expected_contigs = ["ATCGCG", "ATCGCG"]
        self.assertEqual(expected_contigs, accurate_contigs)


class WriteToFileTest(unittest.TestCase):
    """Unit tests for the write_to_file function."""

    def test_write_to_file(self):
        """Test writing contigs to a file."""
        # Test data
        contigs = ["ACGT", "GCTA", "TACG"]
        filename = "unittest.fasta"

        # Call the function
        write_to_file(filename, contigs)

        # Check if the file was created
        self.assertTrue(os.path.exists(filename))

        # Read the file content
        with open(filename, "r", encoding="utf-8") as file:
            file_content = file.read()

        # Check if the file content matches the expected output
        expected_output = ">contig_1\nACGT\n>contig_2\nGCTA\n>contig_3\nTACG\n"
        self.assertEqual(file_content, expected_output)

        # Clean up - delete the test file
        os.remove(filename)


if __name__ == "__main__":
    unittest.main()
