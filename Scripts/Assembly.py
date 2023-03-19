#!usr/bin/env python3

"""
Assembler script
"""

# Imports
import re
import subprocess
import sys
import argparse as ap

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from paf_reader import parse_paf

# Code

reads = []
sequences = []


# Build in possible class if needed
def file_reader(input_file):
    """ Reads through the Fastq file and extracts the sequences and read numbers. """
    for line in open(input_file):
        # Isolates the sequences.
        sequence = re.compile("^[A-Z]{5,}")
        # Isolates the read numbers.
        read = re.compile("read=\d*")
        # Finds all the matches.
        match = (read.findall(line.strip()))
        seqmatch = (sequence.findall(line.strip()))
        # Removes empty matches
        if match != [] or seqmatch != []:
            reads.append(''.join(match))
            sequences.append(''.join(seqmatch))
        while '' in reads and '' in sequences:
            reads.remove("")
            sequences.remove("")
    # Orders the matches with each other.
    for index in range(len(reads)):
        myread = reads[index]
        # Splits sequence in lines of 70 nucleotides.
        myseq = (re.sub("(.{70})", "\\1\n", sequences[index], 0, re.DOTALL))

        # Prints the sequence (REMOVE LATER)
        print(f"{myread} \n {myseq}".replace(' ', ''))


def Alignment(input_file):
    """ Uses minimap two to align the reads """
    filenames = []
    for files in args.fastq_files:
        filenames.append(files.split('.')[0])

    # Runs the overlaps for all files.
    for index in range(len(filenames)):
        myfile = filenames[index]
        subprocess.call(f"minimap2/minimap2 -x ava-ont {input_file} {input_file} > {myfile}_overlaps.paf", shell=True)

        # Creates a CSV file of the overlaps with the proper row names.
        with open(f"{myfile}_overlaps.paf", "r+") as handle:
            df = parse_paf(handle, dataframe=True)
            df.to_csv(f"{myfile}_overlaps.csv", index=False)

def msa(file, output):
    sequences = []
    for line in open(file):
        if line.startswith("A") or line.startswith("C") or line.startswith("G") or line.startswith("T"):
            sequences.append(line.strip())

    longest_length = max(len(s) for s in sequences)
    padded_sequences = [s.ljust(longest_length, '-') for s in sequences]
    records = (SeqRecord(Seq(s), id="alpha") for s in padded_sequences)

    SeqIO.write(records, output, "fasta")

    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=file, out=output)
    print(cline)

def contig_creating():
    """ Works on the contigs of the reads """
    # Make reads into contigs
    # Extend read if it is a repeat to make it anchor somewhere for better contigs
    # Create super contigs
    # Create scaffolds


def Assembly():
    """ Builds the genome assembly """
    # Use the contigs and map against the original 'genome' OR
    # Recreate the original genome through the contigs (research needed)


def main():
    for files in args.fastq_files:
        Alignment(files)
    #msa("b1_1.fq", "output.fasta")

    # file_reader('test.txt')


if __name__ == '__main__':
    argparser = ap.ArgumentParser(description="Arguments for the Assembly")
    argparser.add_argument("-n", action="store",
                           dest="n", required=True, type=int,
                           help="Amount of cores to be used")
    argparser.add_argument("fastq_files", action="store",
                           nargs='+', help="At least one Minion file")
    argparser.add_argument("-k", action="store", dest="k", required=False, type=int,
                           default=5, help="Size of the k-mers")
    args = argparser.parse_args()

    main()
else:
    sys.exit("Program is ending")
