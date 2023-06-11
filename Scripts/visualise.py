#!/usr/bin/env python3

"""
visualise.py script to plot read lengths, GC content, quality scores, 
and sequence complexity of fastq files.
"""

__author__ = "Lisan Eisinga"
__version__ = "2.4.2"
__date__ = "11-06-2023"

import collections
import logging
import warnings
import numpy as np
from Bio import SeqIO, BiopythonDeprecationWarning
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
    # Create a new figure and axes for each plot
    plt.figure()

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
    logging.info("Read length plot of %s saved", fastq_file)


def plot_gc_content(fastq_file):
    """
    Plots the GC content distribution of a FastQ file.

    Parameters:
    fastq_file (str): The path to the FastQ file.

    Returns:
    None
    """
    # Create a new figure and axes for each plot
    plt.figure()

    # Ignores the DeprecationWarning
    warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)

    # Calculates the GC percentage for each record in the FastQ file.
    gc_content = []
    for record in SeqIO.parse(fastq_file, "fastq"):
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
    axis.set_title(f"GC Content Distribution of {fastq_file}")
    axis.set_ylim(0, max(axis.get_yticks()) * 1.1)

    # Creates the normal distribution.
    mean, std = norm.fit(gc_content)
    probability_density = np.linspace(mean - 3 * std, mean + 3 * std, 1000)
    probability = norm.pdf(probability_density, mean, std)
    axis.plot(probability_density, probability, "k", linewidth=2)
    axis.set_xlim(20, 50)

    # Saves the plot as a PNG file and displays it.
    plt.savefig(f"gc_content_{fastq_file}.png")
    logging.info("GC plot of %s saved", fastq_file)

    # Prints the average GC percentage.
    logging.info("Average GC percentage: %s", round(avg_gc, 2))


def plot_quality_scores(fastq_file):
    """
    Plots the quality score frequencies of a FastQ file.

    Parameters:
    fastq_file (str): The path to the FastQ file.

    Returns:
    None
    """
    # Create a new figure and axes for each plot
    plt.figure()

    # Initialize variables
    counts = {}
    total_score = 0
    total_bases = 0

    # Parse the FastQ file and count the quality scores
    for seqrecord in SeqIO.parse(fastq_file, "fastq"):
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
    clrs = plt.cm.get_cmap("gray")(np.linspace(0, 1, len(mydict)))
    clrs[0:25] = (0.2, 0.2, 0.2, 1.0)
    barlist = axis.bar(range(len(sorted(mydict))), list(mydict.values()),
                       align="center", color=clrs)

    # Color the bars based on their quality score
    for i, mybar in enumerate(barlist):
        if i < 20:
            mybar.set_color("black")
        else:
            mybar.set_color("green")

    # Add a legend and labels to the plot
    accepted_patch = plt.Rectangle((0, 0), 1, 1, fc="green")
    rejected_patch = plt.Rectangle((0, 0), 1, 1, fc="black")
    axis.legend([accepted_patch, rejected_patch], ["Score above 20", "Score below 20"],
                loc="upper right", fontsize=12)
    axis.set_xticks(range(len(mydict)))
    axis.set_xticklabels(list(mydict.keys()), fontsize=12)
    axis.set_title(f"The quality score frequencies of {fastq_file} (MinION quality score)",
                   fontsize=16)
    axis.set_xlabel("Scores of the bases", fontsize=14)
    axis.set_ylabel("Total amount of bases", fontsize=14)
    axis.grid(True)

    # Save and show the plot
    plt.savefig(f"Quality_{fastq_file}.png")
    logging.info("Quality plot of %s saved", fastq_file)


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


def plot_sequence_complexity(fastq_file):
    """
    Calculate and plot the sequence complexity of each sequence in a fastq file.

    Parameters:
    fastq_file (str): The path of a fastq file.

    Returns:
    None
    """
    # Create a new figure and axes for each plot
    plt.figure()
    complexities = []
    for seqrecord in SeqIO.parse(fastq_file, "fastq"):
        complexities.append(calculate_sequence_complexity(str(seqrecord.seq)))

    # Plot histogram of sequence complexities
    sns.histplot(complexities, bins=50, kde=True, color="black", edgecolor="black",
                 linewidth=0.5)
    plt.title(f"The sequence complexity of {fastq_file}", fontsize=16)
    plt.xlabel("Sequence complexity", fontsize=14)
    plt.ylabel("Density", fontsize=14)
    sns.set_style("whitegrid")

    # Save plot to file and display it
    plt.savefig(f"sequence_complexity_{fastq_file}.png")
    logging.info("Sequence complexity plot of %s saved", fastq_file)

    # Print average sequence complexity
    avg_complexity = np.mean(complexities)
    logging.info("The average sequence complexity is %s", round(float(avg_complexity), 2))


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
    logging.info("Number of duplicate sequences in %s: %s", fastq_file, num_duplicates)
    return num_duplicates
