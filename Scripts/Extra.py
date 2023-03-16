#!usr/bin/env python3

"""
Script for the gathering of information for the report
"""


def fileinformation(inputfile):
    """ Finds the amount of bases, GC percentage and the longest read. """
    count = { "A": 0, "C": 0, "T": 0, "G": 0 }
    longest_read = 0
    for line in open(inputfile):
        if line.startswith("A") or line.startswith("C") or line.startswith("T") or line.startswith("G"):
            if len(line) > longest_read:
                longest_read = len(line)
            for item in line.strip():
                count[item] += 1

    totalbasepairs = int(count['A'] + count['C'] + count['T'] + count['G'])
    GC = (int(count['G'] + count['C']) / int(totalbasepairs)) * 100
    print(f"The longest read in the file is {longest_read} basepairs long")
    print(f"The total amount of basepairs is: {totalbasepairs}")
    print(f"The total GC percentage is: {GC}%")
