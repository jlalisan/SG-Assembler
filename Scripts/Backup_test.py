# Test file and backup file.
import random
import subprocess
import sys
import argparse as ap


def create_fake_genome(length):
    """ Creates a random sequence with a fake genome header"""
    fake_genome = ""
    with open("fake_genome.txt", "w") as writeto:
        writeto.write(">Fake_header_for_fake_genome" + "\n")
        for i in range(length):
            fake_genome += random.choice(["A", "C", "T", "G"])
            # Writes it to the new fake genome file in lines of length 50
        writeto.write('\n'.join(fake_genome[i:i + 50] for i in range(0, len(fake_genome), 50)))


create_fake_genome(30000)


def Alignment():
    """ Uses minimap two to align the reads """
    subprocess.call("minimap2/minimap2 -x ava-ont reads_1.fasta reads_2.fasta > overlaps.paf ", shell=True)
    # Need to check which reads go in
    # What is inside the overlaps.paf file

    # Extract useful information from the file
    # Prepare for read extraction and adjust some reads to make them longer in case of repeats


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
    file_reader('test.txt')


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
