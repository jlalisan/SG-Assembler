#!usr/bin/env python3

"""
Script for the gathering of information for the report
"""
import math
import os
import re

import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
from numpy import mean


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

def fastq_to_dataframe(filename, size=100000000):
    """Convert fastq to dataframe.
        size: limit to the first reads of total size
        Returns: dataframe with reads
    """

    ext = os.path.splitext(filename)[1]
    if ext=='.fastq' or ext=='.gz':
        fastq_parser = SeqIO.parse(open(filename, "rt"), "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    i=0
    res=[]
    for fastq_rec in fastq_parser:
        #print (fastq_rec.seq)
        i+=1
        if i>size:
            break
        res.append([fastq_rec.id, str(fastq_rec.seq)])
    df = pd.DataFrame(res, columns=['id','seq'])
    df['length'] = df.seq.str.len()
    return df

def normpdf(x, mean, sd):
    """sample a normal distribution at given point"""

    var = float(sd)**2
    denom = (2*math.pi*var)**.5
    num = math.exp(-(float(x)-float(mean))**2/(2*var))
    return (num/denom)

def plot_fastq_gc_content(filename, ax=None, limit=5000000):
    from Bio.SeqUtils import GC
    if ax is None:
        f,ax=plt.subplots(figsize=(12,5))
    df = fastq_to_dataframe(filename, size=limit)
    gc = df.seq.apply(lambda x: GC(x))
    gc.hist(ax=ax,bins=150,color='black',grid=False,histtype='step',lw=2)
    ax.set_xlim((0,100))
    x=np.arange(1,100,.1)
    f = [normpdf(i, gc.mean(), gc.std()) for i in x]
    ax2=ax.twinx()
    ax2.plot(x,f)
    ax2.set_ylim(0,max(f))
    ax.set_title('GC content',size=15)
    plt.show()
    return

def ave_qual(fastq):
    allscores = []
    for record in SeqIO.parse(fastq, "fastq"):
        pattern = re.compile("read=\d*")
        readnum = re.findall(pattern, record.description)
        allscores.append(record.letter_annotations['phred_quality'])
        #print(f"The average score of {''.join(readnum)} is {mean(record.letter_annotations['phred_quality'])}")

    totalavg = 0.0
    for scores in allscores:
        score = mean(scores)
        totalavg += score

    print(f"the average score for file {fastq} is {totalavg / len(allscores)}")