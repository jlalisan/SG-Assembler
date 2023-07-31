#!/usr/bin/env python3

"""Generate an artificial FastQ file with random data.

Usage:
  generate_artificial_fastq.py <num_reads> <read_length> [--output FILE]

Options:
  -h --help            Show this help message.
  -o --output FILE     Specify the output file name [default: artificial_file.fastq].
"""

import random
from docopt import docopt


def generate_random_sequence(length: int) -> str:
    """
    Generate a random nucleotide sequence of a given length.

    Parameters:
        length (int): The length of the sequence to generate.

    Returns:
        str: A random nucleotide sequence.
    """
    nucleotides = 'ACGT'
    return ''.join(random.choice(nucleotides) for _ in range(length))


def generate_minion_quality(length: int) -> str:
    """Generate MinION-like quality scores for a sequence of a given length.

    Args:
        length (int): The length of the quality scores to generate.

    Returns:
        str: MinION-like quality scores as ASCII characters.
    """
    return ''.join(chr(random.randint(0, 93) + 33) for _ in range(length))


def generate_artificial_fastq(num_reads: int, read_length: int) -> list:
    """
    Generate artificial FastQ data.

    Parameters:
        num_reads (int): The number of reads to generate.
        read_length (int): The length of each read.

    Returns:
        list: A list of artificial FastQ entries.
    """
    artificial_fastq = []

    for i in range(num_reads):
        read_name = f'@Read{i + 1}'
        sequence = generate_random_sequence(read_length)
        quality = generate_minion_quality(read_length)

        fastq_entry = f"{read_name}\n{sequence}\n+\n{quality}\n"
        artificial_fastq.append(fastq_entry)

    return artificial_fastq


def main():
    """
    Generate an artificial FastQ file with the given number of reads and read length.

    The script generates random FastQ data with the specified number of reads and read length.
    The resulting data is saved to a FastQ file specified by the optional '--output' argument.

    Arguments:
        <num_reads> (int): The number of reads to generate.
        <read_length> (int): The length of each read.

    Options:
        --output FILE (str): Specify the output file name [default: artificial_file.fastq].
    """
    args = docopt(__doc__)

    num_reads = int(args['<num_reads>'])
    read_length = int(args['<read_length>'])
    output_file = args['--output']

    artificial_fastq_data = generate_artificial_fastq(num_reads, read_length)

    with open(output_file, "w", encoding="utf-8") as file_handle:
        file_handle.write(''.join(artificial_fastq_data))

    print(f"Generated {num_reads} artificial reads of length {read_length} in '{output_file}'.")


if __name__ == "__main__":
    main()
