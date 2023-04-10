# Test file and backup file.
from __future__ import annotations
from __future__ import annotations

import random
import subprocess
import sys
import argparse as ap

import difflib
import itertools
import os
import random
import subprocess
import re
from statistics import mode

import pandas
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline

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



# Converts Fastq files into Fasta files
def convert_to_fasta(inputfile, outputfile):
    with open(inputfile, "rU") as input_handle:
        with open(outputfile, "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "fastq")
            count = SeqIO.write(sequences, output_handle, "fasta")

    print("Converted %i records" % count)
    return output_handle


# convert_to_fasta("b1_1.fq", "reads_1.fasta")
def fix_headers(inputfile):
    read = re.compile("read=\d*")
    matched = []
    seqs = []
    with open("testing.txt", "w") as testfile:
        with open(inputfile, "r+") as rr:
            for line in rr:
                if "read" in line:
                    match = read.findall(line)
                    matched.append(line)
                    testfile.write(line.replace(line, ">" + "".join(match)))
                    testfile.write("\n")
                else:
                    seqs.append(line)
                    testfile.write(line)


# fix_headers("reads_1.fasta")
def create_fake_genome(length):
    fake_genome = ""
    with open("fake_genome.fa", "w") as writeto:
        writeto.write(">Fake_header_for_fake_genome" + "\n")
        for i in range(length):
            fake_genome += random.choice(["A", "C", "T", "G"])
        writeto.write('\n'.join(fake_genome[i:i + 50] for i in range(0, len(fake_genome), 50)))


# create_fake_genome(30000)
def minimapper():
    """ Uses minimap two to align the reads """
    subprocess.call("minimap2/minimap2 -x ava-ont testing.txt testing.txt > overlaps.paf ", shell=True)


# minimapper()
def paffer():
    from paf_reader import parse_paf
    with open("overlaps.paf", "r+") as handle:
        df = parse_paf(handle, dataframe=True)
        df.to_csv("overlaps.csv", index=False)


# paffer()
def basepaircalling(inputfile):
    count = { "A": 0, "C": 0, "T": 0, "G": 0 }
    longest_read = 0
    reads = []
    temp = ""
    for line in open(inputfile):
        if line.startswith("A") or line.startswith("C") or line.startswith("T") or line.startswith("G"):
            reads.append(line)
            temp = [len(ele) for ele in reads]

            if len(line) > longest_read:
                longest_read = len(line)
            for item in line.strip():
                count[item] += 1

    totalbasepairs = int(count['A'] + count['C'] + count['T'] + count['G'])
    GC = (int(count['G'] + count['C']) / int(totalbasepairs)) * 100
    res = [0 if len(temp) == 0 else (float(sum(temp)) / len(temp))]

    print(f"The Average length of reads is: {''.join(str(res))}")
    print(f"The longest read in the file is: {longest_read} basepairs long")
    print(f"The total amount of basepairs is: {totalbasepairs}")
    print(f"The percentage of GC is: {GC}%")


def get_sequences(file):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    sequences_list = []
    read_line = []
    with open(file) as f:
        for line in f:
            if line.startswith("A") or line.startswith("T") or line.startswith("C") or line.startswith("G"):
                sequences_list.append(line.strip())
            else:
                if line.startswith("@"):
                    read = re.compile("read=\d*")
                    read_line.append(("".join(read.findall(line))))

    for index in range(len(read_line)):
        longest_length = max(len(s) for s in sequences_list)
        padded_sequences = [s.ljust(longest_length, '-') for s in sequences_list]
        records = (SeqRecord(Seq(s), read_line[index], description=f"reads for line {read_line[index]}") for s in
                   padded_sequences)

        SeqIO.write(records, "msa.txt", "fasta")

    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input="testing.txt", out="msa.txt")
    print(cline)


# get_sequences("test.txt")
# basepaircalling("b1_1.fq")
def aligning_test(input_file):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align.Applications import MuscleCommandline

    sequences = []

    longest_length = max(len(s) for s in sequences)
    padded_sequences = [s.ljust(longest_length, '-') for s in sequences]
    records = (SeqRecord(Seq(s), id="alpha") for s in padded_sequences)

    SeqIO.write(records, "msa.txt", "fasta")

    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input="testing.txt", out="msa.txt")
    print(cline)


# aligning_test('test.txt')

def msa(file, output):
    sequences = []
    ids = []
    read = re.compile("read=\d*")
    for line in open(file):
        if line.startswith("A") or line.startswith("C") or line.startswith("G") or line.startswith("T"):
            sequences.append(line.strip())

        for item in read.findall(line):
            ids.append(''.join(read.findall(item)))
    # Needs fixing
    for index in range(len(ids)):
        alpha = ids[index]

    longest_length = max(len(s) for s in sequences)
    padded_sequences = [s.ljust(longest_length, '-') for s in sequences]
    records = (SeqRecord(Seq(s), id=f"{ids[index]}") for s in padded_sequences)
    #
    SeqIO.write(records, output, "fasta")
    #
    from Bio.Align.Applications import MuscleCommandline
    cline = MuscleCommandline(input=file, out=output)
    print(cline)


# msa("b1_1.fq", "output.fasta")

from pathlib import Path
from itertools import chain
import Bio
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Align.AlignInfo import SummaryInfo


#
# SeqRecord = Bio.SeqRecord.SeqRecord


def get_consensus_seq(filename: Path | str) -> SeqRecord:
    common_alignment = MultipleSeqAlignment(
        chain(*AlignIO.parse(filename, "fasta"))
    )
    summary = SummaryInfo(common_alignment)
    consensus = summary.dumb_consensus(0.7, "N")
    with open("consensus.fasta", 'w') as test:
        test.write(str(consensus))
    return consensus


# get_consensus_seq('output.fasta')

def make_file(fastq):
    readnum = []
    reads = []
    with open("reads_1.fasta", "w") as writeto:
        for line in open(fastq):
            readnumber = re.compile("read=\d*")
            readnum.append("".join(readnumber.findall(line)))
            read = re.compile("^[A-Z]{5,}")
            reads.append("".join(read.findall(line)))

        readnum = list(filter(None, readnum))
        reads = list(filter(None, reads))

        for index in range(len(readnum)):
            writeto.write(readnum[index] + "\n")
            writeto.write(reads[index] + "\n")
    writeto.close()

    myfile = open("reads_1.fasta")
    for index in range(len(readnum)):
        readnumber = myfile.readline()
        readseq = myfile.readline()


# make_file("b1_1.fq")


def sequence_ripper(file1, file2):
    for line in open(file1):
        startpos = line.split("\t")[2].strip()
        endpos = line.split("\t")[3].strip()
        targetpos = line.split("\t")[7].strip()
        targetend = line.split("\t")[8].strip()
        readnumquery = line.split('\t')[0].strip()
        readnumtarget = line.split('\t')[5].strip()
        totalreads = []

        readybeedy = re.compile("read=\d*")
        totalreads.append("".join(readybeedy.findall(line)))
        totalreads = list(filter(None, totalreads))

        myfile = open("reads_1.fasta")
        step = 2
        for lineno, line in enumerate(myfile):
            if lineno % step == 0:
                readline = line.strip()
            if lineno % step == 1:
                readseq = line.strip()
                if readline == readnumquery.strip():
                    print(readnumquery)
                elif readline == readnumtarget.strip():
                    print(readnumtarget)


sequence_ripper("overlaps.paf", "reads_1.fasta")


def overlapmatch(overlaps, fastq):
    for line in open(overlaps):
        startpos = line.split("\t")[2].strip()
        endpos = line.split("\t")[3].strip()
        targetpos = line.split("\t")[7].strip()
        targetend = line.split("\t")[8].strip()
        readnumquery = line.split('\t')[0].strip()
        readnumtarget = line.split('\t')[5].strip()
        for_rev = line.split("\t")[4].strip()

    for reads in open(fastq):
        # if read num starts with readnum query:
        # if for_rev = +:
        # sequence.append (READNAME) readnumquary[startpos:endpos] to readnumtarget[startpos:endpos]
        # else:
        # sequence.append (READNAME) readnumquary[startpos:endpos] to readnumtarget.reverse[startpos:endpos]
        # Possibly add them to a sequence list, but it may get long
        pass

    # If sequence is + then append to map in forward way
    # Else append to map in reverse way
    # Make sure to append from the starting point in the sequence up until the ending point
    # Only keep appending to read until the first readnum is different, then start new contig
    # Keep first readnum name as contig name
    #

    """
    THIS WORKS BUT WITHOUT THE FORWARD OR REVERSE BUSINESS.
    sequences = ["ATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACC",
                 "ATCATCAACTGGTGGTGAAATGACTGGGCAAGTGCTGTTGGTGCTGGAAATAA",
                 "AGTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAAGACACCGACAACTTTCTT",
                 "AAGTCTTTTGTCCTTCCTTCTTTCCAAAAGCATCTGACTTCTTAACTAG"]
    # Tranpose the list to group characters based on their position in words
    # Remove None values and sort alphabetically and use the mode function from statistics module
    totalmatch = "".join((list(map(lambda x: mode(sorted(filter(lambda v: v != None, x))), (itertools.zip_longest(*sequences))))))
    print(totalmatch)
    """


# overlapmatch()


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

def overlapmatch(overlaps, fastq):
    """
    Okay, ik moet de start positie pakken en die opzoeken in de fastq, vanaf de start tot de end de sequentie eruit trekken
    dan zelfde voor de target, en die samenvoegen onder read=73: ATCTGT enzo en dat tot de volgende read begin
    """

    endpos = []
    targetpos = []
    targetend = []
    readnumq = []
    readnumt = []
    for_rev = []
    startpos = []

    sequences = []
    complemented = {}
    shortseqs = {}

    for line in open(overlaps):
        line = line.strip().split("\t")
        endpos.append(line[3])
        targetpos.append(line[7])
        targetend.append(line[8])
        readnumq.append(line[0])
        readnumt.append(line[5])
        for_rev.append(line[4])
        startpos.append(line[2])

    tada = open(fastq).readlines()

    for index in range(len(startpos)):
        for i, l in enumerate(tada):
            if readnumq[index] in l and for_rev[index] == "+":
                shortseqs[readnumq[index]] = tada[i+1].strip()
            if readnumq[index] in l and for_rev[index] == "-":
                complemented[readnumq[index]] = tada[i+1][::-1].strip()
        for line

    #print(shortseqs)
    #print(complemented)
    # if readnumquary is the same as readnum(in reads):
    # read
    # if read num starts with readnum query:
    # if for_rev = +:
    # sequence.append (READNAME) readnumquary[startpos:endpos] to readnumtarget[startpos:endpos]
    # else:
    # sequence.append (READNAME) readnumquary[startpos:endpos] to readnumtarget.reverse[startpos:endpos]
    # Possibly add them to a sequence list, but it may get long
    pass

    # If sequence is + then append to map in forward way
    # Else append to map in reverse way
    # Make sure to append from the starting point in the sequence up until the ending point
    # Only keep appending to read until the first readnum is different, then start new contig
    # Keep first readnum name as contig name
    #

    # THIS WORKS BUT WITHOUT THE FORWARD OR REVERSE BUSINESS.
    # sequences = ["ATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACC",
    #             "ATCATCAACTGGTGGTGAAATGACTGGGCAAGTGCTGTTGGTGCTGGAAATAA",
    #             "AGTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAAGACACCGACAACTTTCTT",
    #             "AAGTCTTTTGTCCTTCCTTCTTTCCAAAAGCATCTGACTTCTTAACTAG"]
    ## Tranpose the list to group characters based on their position in words
    ## Remove None values and sort alphabetically and use the mode function from statistics module
    # totalmatch = "".join((list(map(lambda x: mode(sorted(filter(lambda v: v is not None, x))), (itertools.zip_longest(*sequences))))))
    # print(totalmatch)


overlapmatch("overlaps.paf", "b1_1.fq")