#!/usr/bin/env python3
"""
Contig Writing

This module contains a function to write contigs to a Fasta data set,
with individual contig numbering.

Author: Lisan Eisinga
Date: 26-07-2023
"""

import logging


def write_to_file(filename: str, contigs: list) -> None:
    """
    Writes the contigs to a Fasta data set with individual contig numbering.

    Parameters:
        filename: The path to the output file.
        contigs: A list of contigs to be written.

    Raises:
        IOError: If there is an error while writing to the data set.

    Returns:
        None
    """
    try:
        with open(filename, "w", encoding="utf-8") as file:
            for i, contig in enumerate(contigs):
                header = f">contig_{i + 1}"
                sequence = "".join(contig)
                file.write(f"{header}\n{sequence}\n")
    except IOError as error:
        logging.error("Error: Unable to write to data set %s. %s", filename, str(error))
