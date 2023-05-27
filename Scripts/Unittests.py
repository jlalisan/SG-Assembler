import os
import pandas as pd
import networkx as nx
import unittest
from Attempt_7 import read_fastq_file, parse_paf, overlap_graph

class ReadFastqFileTest(unittest.TestCase):
    def setUp(self):
        self.valid_fastq_file = 'valid.fastq'
        self.invalid_fastq_file = 'invalid.fastq'

        # Create a valid FastQ file for testing
        with open(self.valid_fastq_file, 'w') as f:
            f.write('@seq1\nACGT\n+\nIIII\n@seq2\nTGCA\n+\nJJJJ\n')

    def tearDown(self):
        # Remove the temporary FastQ file
        os.remove(self.valid_fastq_file)

    def test_read_fastq_file_valid_file(self):
        # Test case for a valid FastQ file
        expected_sequences = {'seq1': 'ACGT', 'seq2': 'TGCA'}
        result = read_fastq_file(self.valid_fastq_file)
        self.assertEqual(result, expected_sequences)

    def test_read_fastq_file_invalid_file(self):
        # Test case for a non-existent FastQ file
        with self.assertRaises(FileNotFoundError):
            read_fastq_file(self.invalid_fastq_file)


class ParsePafTest(unittest.TestCase):
    def setUp(self):
        self.valid_paf_file = 'valid.paf'
        self.invalid_paf_file = 'invalid.paf'

        # Create a valid PAF file for testing
        with open(self.valid_paf_file, 'w') as f:
            f.write('query1\t100\t0\t50\t+\ttarget1\t200\t50\t100\nquery2\t150\t0\t100\t-\ttarget2\t300\t0\t150\n')

    def tearDown(self):
        # Remove the temporary PAF file
        os.remove(self.valid_paf_file)

    def test_parse_paf_valid_file(self):
        # Test case for a valid PAF file
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
        # Test case for a non-existent PAF file
        with self.assertRaises(FileNotFoundError):
            parse_paf(self.invalid_paf_file)

class OverlapGraphTest(unittest.TestCase):
    def test_overlap_graph_valid_input(self):
        # Test case for valid input
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
        expected_graph.add_edge('seq1', 'seq2', query_start=0, query_end=3, target_start=0, target_end=3,
                                strand='+', query_length=4, target_length=4, overlap_len=4)
        expected_graph.add_edge('seq2', 'seq1', query_start=0, query_end=3, target_start=0, target_end=3,
                                strand='-', query_length=4, target_length=4, overlap_len=4)
        expected_graph.add_edge('seq2', 'seq3', query_start=0, query_end=3, target_start=0, target_end=3,
                                strand='+', query_length=4, target_length=4, overlap_len=4)

        result_graph = overlap_graph(sequences, overlaps_df)
        self.assertEqual(result_graph.nodes, expected_graph.nodes)
        self.assertEqual(result_graph.edges, expected_graph.edges)

    def test_overlap_graph_missing_columns(self):
        # Test case for missing columns in overlaps DataFrame
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
        # Test case for invalid sequence ID in overlaps DataFrame
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

if __name__ == '__main__':
    unittest.main()
