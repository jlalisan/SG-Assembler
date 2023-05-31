#!/usr/bin/env python3

"""
visualise.py script to plot read lengths, GC content, quality scores, 
and sequence complexity of fastq files.
"""

__author__ = "Lisan Eisinga"
__version__ = "2.4.0"
__date__ = "27-05-2023"

import collections
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC
from matplotlib import pyplot as plt
from scipy.stats import norm
import seaborn as sns


def plot_read_lengths(fastq_file):
    """
    Plots a histogram of the read lengths in a FastQ file.

    Parameters:
    fastq_file (str): The path to the FastQ file.

    Returns:
    None
    """
    # Extracts the read lengths from the FastQ file
    read_lengths = [len(record.seq) for record in SeqIO.parse(fastq_file, "fastq")]

    # Creates a histogram of the read lengths
    _, axis = plt.subplots(figsize=(8, 6))
    axis.hist(read_lengths, bins=50, color="#1f77b4", alpha=0.8, edgecolor="black")

    # Sets the labels and title of the plot
    axis.set_xlabel("Read Length", fontsize=14)
    axis.set_ylabel("Count", fontsize=14)
    axis.set_title(f"Read length distribution of {fastq_file}", fontsize=16)

    # Sets the tick size and removes the top and right spines
    axis.tick_params(axis="both", which="major", labelsize=12)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    # Sets the y-axis to a logarithmic scale
    axis.set_yscale("log")

    # Saves the plot with a name based on the FastQ file name and displays it
    plt.savefig(f"read_lengths_{fastq_file}.png")
    plt.show()


def plot_gc_content(filename):
    """
    Plots the GC content distribution of a FastQ file.

    Parameters:
    filename (str): The path to the FastQ file.

    Returns:
    None
    """
    # Calculates the GC percentage for each record in the FastQ file.
    gc_content = []
    for record in SeqIO.parse(filename, "fastq"):
        gc_content.append(GC(record.seq))

    # Calculates the average GC percentage.
    avg_gc = sum(gc_content) / len(gc_content)

    # Plots the histogram of GC content distribution.
    _, axis = plt.subplots()
    axis.hist(gc_content, bins=np.arange(0, 101, 1), density=True,
              edgecolor="black", linewidth=0.5)
    axis.set_xticks(np.arange(0, 101, 5))
    axis.set_xlabel("GC Percentage")
    axis.set_ylabel("Density")
    axis.set_title(f"GC Content Distribution of {filename}")
    axis.set_ylim(0, max(axis.get_yticks()) * 1.1)

    # Creates the normal distribution.
    mu, std = norm.fit(gc_content)
    x = np.linspace(mu - 3 * std, mu + 3 * std, 1000)
    p = norm.pdf(x, mu, std)
    axis.plot(x, p, "k", linewidth=2)
    axis.set_xlim(20, 50)

    # Saves the plot as a PNG file and displays it.
    plt.savefig(f"gc_content_{filename}.png")
    plt.show()

    # Prints the average GC percentage.
    print("Average GC percentage:", avg_gc)


def plot_quality_scores(fastq):
    """
    Plots the quality score frequencies of a FastQ file.

    Parameters:
    fastq (str): The path to the FastQ file.

    Returns:
    None
    """
    # Initialize variables
    counts = {}
    total_score = 0
    total_bases = 0

    # Parse the FastQ file and count the quality scores
    for seqrecord in SeqIO.parse(fastq, "fastq"):
        for qual in seqrecord.letter_annotations["phred_quality"]:
            if qual not in counts:
                counts[qual] = qual
            else:
                counts[qual] += 1
            total_score += qual
            total_bases += 1

    # Sort the quality scores and plot them
    mydict = collections.OrderedDict(sorted(counts.items()))
    _, axis = plt.subplots(figsize=(18, 5))
    clrs = plt.cm.Greys(np.linspace(0, 1, len(mydict)))
    clrs[0:25] = (0.2, 0.2, 0.2, 1.0)
    barlist = axis.bar(range(len(sorted(mydict))), list(mydict.values()),
                       align="center", color=clrs)

    # Color the bars based on their quality score
    for i in range(len(barlist)):
        if i < 20:
            barlist[i].set_color("black")
        else:
            barlist[i].set_color("green")

    # Add a legend and labels to the plot
    accepted_patch = plt.Rectangle((0, 0), 1, 1, fc="green")
    rejected_patch = plt.Rectangle((0, 0), 1, 1, fc="black")
    axis.legend([accepted_patch, rejected_patch], ["Score above 20", "Score below 20"], loc="upper right", fontsize=12)
    axis.set_xticks(range(len(mydict)))
    axis.set_xticklabels(list(mydict.keys()), fontsize=12)
    axis.set_title(f"The quality score frequencies of the {fastq} file", fontsize=16)
    axis.set_xlabel("Scores of the bases", fontsize=14)
    axis.set_ylabel("Total amount of bases", fontsize=14)
    axis.grid(True)

    # Save and show the plot
    plt.savefig(f"Quality_{fastq}.png")
    plt.show()


def calculate_sequence_complexity(seq):
    """
    Calculate the sequence complexity of a DNA sequence.

    Parameters:
    seq (str): A DNA sequence.

    Returns:
    float: The sequence complexity, which is the number of unique 4-mers
    in the sequence divided by the length of the sequence minus 3.
    """
    kmers = set()
    for i in range(len(seq) - 3):
        kmers.add(seq[i:i + 4])
    return len(kmers) / (len(seq) - 3)


def plot_sequence_complexity(fastq):
    """
    Calculate and plot the sequence complexity of each sequence in a fastq file.

    Parameters:
    fastq (str): The path of a fastq file.

    Returns:
    None
    """
    complexities = []
    for seqrecord in SeqIO.parse(fastq, "fastq"):
        complexities.append(calculate_sequence_complexity(str(seqrecord.seq)))

    # Plot histogram of sequence complexities
    sns.histplot(complexities, bins=50, kde=True, color="black", edgecolor="black",
                 linewidth=0.5)
    plt.title(f"The sequence complexity of the {fastq} file", fontsize=16)
    plt.xlabel("Sequence complexity", fontsize=14)
    plt.ylabel("Density", fontsize=14)
    sns.set_style("whitegrid")

    # Save plot to file and display it
    plt.savefig(f"sequence_complexity_{fastq}.png")
    plt.show()

    # Print average sequence complexity
    avg_complexity = np.mean(complexities)
    print(f"The average sequence complexity is {avg_complexity:.2f}")


def count_duplicates(fastq_file):
    """
    Count the number of duplicate sequences in a FASTQ file.

    Parameters:
        fastq_file (str): Path to the FASTQ file.

    Returns:
        int: Number of duplicate sequences in the FASTQ file.
    """
    # Count the number of occurrences of each sequence
    seq_counts = collections.defaultdict(int)
    for record in SeqIO.parse(fastq_file, "fastq"):
        seq_counts[str(record.seq)] += 1

    # Count the number of duplicate sequences
    num_duplicates = sum(count for count in seq_counts.values() if count > 1)

    # Print the result
    print(f"Number of duplicate sequences in {fastq_file}: {num_duplicates}")
    return num_duplicates
