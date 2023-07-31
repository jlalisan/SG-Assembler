#!/usr/bin/env python3
"""
Genome Assembly Pipeline

This script performs the genome assembly process based on the provided command-line arguments.

Usage:
  main.py <fastq_file> [--paf_file=<paf_file>] [--output_file=<output_file>] [--visualize]
  main.py (-h | --help)

Options:
  -h --help                        Show this help message.
  -p --paf_file=<paf_file>         Path to the PAF data set.
  -o --output_file=<output_file>   Path to the output Fasta data set.
  -v --visualize                   Display information using visualization.py. [default: False]

Author: Lisan Eisinga
Version: 2.0
Date of Completion: 26-07-2023

Input data set formats:
- fastq_file: FastQ data set containing biological sequence and quality score information (MinION).
"""

import logging
import os
import sys
from docopt import docopt
from Assembly import fastq_parser, paf_parser, create_overlap_graph, \
    dfs_traversal, consensus_generator, file_writer
from Visualization import visualization as vis


def assemble_genome(fastq_file: str, paf_file: str, output_file: str, visualize: str) -> None:
    """
    Perform the genome assembly process based on the provided command-line arguments.

    Parameters:
        fastq_file: The path to the FastQ file.
        paf_file: The path to the PAF file.
        output_file: The path to the output Fasta file.
        visualize: A string indicating whether to display information using visualization.py.

    Returns:
        None
    """
    try:
        # Read the FastQ data set
        logging.info("Reading FastQ data set...")
        sequences = fastq_parser.read_fastq_file(fastq_file)

        # Parse the PAF data set
        logging.info("Parsing PAF data set...")
        overlaps = paf_parser.parse_paf(paf_file)

        # Create the overlap graph and remove isolated nodes
        logging.info("Creating overlap graph...")
        overlap_graph = create_overlap_graph.overlap_graph(sequences, overlaps)
        overlap_graph = create_overlap_graph.remove_isolated_nodes(overlap_graph)

        # Set recursion limit to avoid stack overflow error
        if overlap_graph.number_of_nodes() >= 999:
            sys.setrecursionlimit(overlap_graph.number_of_nodes())

        # Traverse the overlap graph to generate contigs
        logging.info("Traversing the graph...")
        contigs = dfs_traversal.dfs(overlap_graph)
        # Generate consensus sequences from the contigs
        if contigs:
            logging.info("Generating consensus sequences...")
            consensus_seq = consensus_generator.generate_consensus_sequences(overlap_graph, contigs)

            # If output file is not provided, create it from the FastQ file
            if not output_file:
                base_name = os.path.splitext(fastq_file)[0]
                output_file = f"{base_name}_contigs.fasta"

            # Write the consensus sequences to the output file
            logging.info("Writing contigs to %s...", output_file)
            file_writer.write_to_file(output_file, consensus_seq)

        # Visualize the FastQ data set if requested
        if visualize:
            logging.info("Visualizing FastQ data set...")
            vis.plot_read_lengths_histogram(fastq_file)
            vis.plot_gc_content_histogram(fastq_file)
            vis.plot_rolling_mean_quality(fastq_file)
            vis.plot_sequence_complexity(fastq_file)
            logging.info("All visualizations are saved in the current working directory.")

    except FileNotFoundError as fnf_error:
        logging.error("File not found: %s", str(fnf_error))
    except ValueError as value_error:
        logging.error("Value error: %s", str(value_error))


def main(args: dict[str, str]) -> None:
    """
    Parse command-line arguments and perform the genome assembly process.

    Parameters:
        args: A dictionary of command-line arguments.

    Returns:
        None
    """
    try:
        # Parse command-line arguments
        fastq_file = args["<fastq_file>"]
        paf_file = args["--paf_file"]
        output_file = args["--output_file"]
        visualize = args["--visualize"]

        # If PAF file is not provided, create it from the FastQ file
        if not paf_file:
            logging.info("Creating PAF data set...")
            paf_file = paf_parser.create_paf(fastq_file)
            logging.info("PAF data set created: %s", paf_file)

        # Perform the genome assembly process
        assemble_genome(fastq_file, paf_file, output_file, visualize)

    except FileNotFoundError as fnf_error:
        logging.error("File not found: %s", str(fnf_error))


if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Parse command-line arguments and run the main function
    arguments = docopt(__doc__)
    main(arguments)
