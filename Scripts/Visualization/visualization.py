#!/usr/bin/env python3

"""
MinION Data Analysis

This script contains functions to perform various analyses on MinION sequencing data,
including plotting histograms of read lengths and GC content, and calculating sequence
complexity using kernel density estimation.

Author: Lisan Eisinga
Version: 1.09
Date of Completion: 26-07-2023

Input data set formats:
- fastq_file: FastQ data set containing biological sequence and quality score information (MinION).

"""

import logging
import numpy as np
import seaborn as sns
from Bio import SeqIO
from matplotlib import pyplot as plt, ticker
from scipy.stats import norm


def plot_read_lengths_histogram(filename):
    """
    Plot the histogram of read lengths with a log scale for the x-axis.

    This function reads a fastQ file, calculates the read lengths for each read,
    and then plots a histogram showing the distribution of read lengths.
    The histogram is displayed with a log scale for the x-axis.

    Parameters:
        filename (str): The filename of the fastQ file containing the sequence reads.

    Returns:
        None.

    Note:
        The function saves the plot as a high-resolution image,
        in the current working directory named 'minion_read_lengths_histogram.png'.
    """

    # Open the FastQ file and calculate the read lengths
    with open(filename, 'r', encoding="utf-8") as filehandle:
        read_lengths = []
        for i, line in enumerate(filehandle):
            if i % 4 == 1:
                read_length = len(line.strip())
                read_lengths.append(read_length)

    # Plot the histogram of read lengths with a log scale for the x-axis
    plt.figure(figsize=(8, 6))
    plt.hist(read_lengths, bins=50, color='steelblue', edgecolor='black', alpha=0.7, log=True)

    # Set axis labels and title with a bold font
    plt.xlabel('Read Length (Base pairs)', fontsize=12, fontweight='bold')
    plt.ylabel('Frequency', fontsize=12, fontweight='bold')
    plt.title('MinION Read Length Distribution', fontsize=14, fontweight='bold')

    # Add grid lines
    plt.grid(color='lightgray', linestyle='--', linewidth=0.5, alpha=0.7)

    # Customize plot border
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)

    # Format read length tick labels with comma separator for thousands
    plt.xticks(fontsize=10)
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter
                                        (lambda temp, pos: '{:,.0f}'.format(temp)))

    # Save the plot as a high-resolution image.
    plt.savefig('minion_read_lengths_histogram.png', dpi=300, bbox_inches='tight')

    plt.clf()


def plot_gc_content_histogram(filename):
    """
    Plot the histogram of GC content for reads in a fastQ file.

    This function reads a fastQ file, calculates the GC content for each read,
    and then plots a histogram showing the distribution of GC content.
    The histogram is fitted with a normal distribution curve.

    Parameters:
        filename (str): The filename of the fastQ file containing the sequence reads.

    Returns:
        None.

    Note:
        The function saves the plot as a high-resolution image in the current working directory
        named 'minion_gc_content_histogram.png'.
    """

    # Open the FastQ file and calculate the GC content for each read
    with open(filename, 'r', encoding="utf-8") as filehandle:
        gc_contents = []
        for i, line in enumerate(filehandle):
            if i % 4 == 1:
                sequence = line.strip()
                gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
                gc_contents.append(gc_content)

    # Plot the histogram of GC content
    plt.figure(figsize=(8, 6))
    plt.hist(gc_contents, bins=25, color='steelblue', edgecolor='black',
             alpha=0.7, density=True, linewidth=1.2)

    # Fit the histogram with a normal distribution
    mean_gc = np.mean(gc_contents)
    std_gc = np.std(gc_contents)
    linespace = np.linspace(min(gc_contents), max(gc_contents), 100)
    fit_line = norm.pdf(linespace, mean_gc, std_gc)
    plt.plot(linespace, fit_line, color='black', linewidth=2, linestyle='--', label='Normal Fit')

    # Set axis labels and title with a bold font
    plt.xlabel('GC Content (%)', fontsize=12, fontweight='bold')
    plt.ylabel('Density', fontsize=12, fontweight='bold')
    plt.title('MinION GC Content Histogram', fontsize=14, fontweight='bold')

    # Add grid lines
    plt.grid(color='lightgray', linestyle='--', linewidth=0.5, alpha=0.7)

    # Customize plot border
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)

    # Format GC content tick labels with one decimal place
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    # Make the histogram bars slightly transparent
    plt.gca().patches[0].set_alpha(0.7)

    # Set appropriate axis limits
    plt.xlim(min(gc_contents) - 2, max(gc_contents) + 2)
    plt.ylim(0, max(fit_line) + 0.03)

    # Remove legend frame
    legend = plt.legend(fontsize=10)
    legend.get_frame().set_linewidth(0)

    # Save the plot as a high-resolution image
    plt.savefig('minion_gc_content_histogram.png', dpi=300, bbox_inches='tight')

    plt.clf()


def plot_rolling_mean_quality(filename, window_size=100):
    """
    Plot the rolling mean quality scores of reads in a fastQ file.

    This function reads a fastQ file, extracts the quality scores for each read,
    calculates the mean quality score for each read,
    and then applies a rolling average to smoothen the plot.
    The rolling mean quality scores are then plotted using Matplotlib.

    Parameters:
        filename (str): The filename of the fastQ file containing the quality scores.
        window_size (int, optional): The size of the window used for the rolling average.
        Default is 100.

    Returns:
        None.

    Note:
        The function saves the plot as a high-resolution image in the current working directory
        named 'minion_quality_plot.png'.
    """

    # Open the FastQ file and extract the quality scores
    with open(filename, 'r', encoding="utf-8") as filehandle:
        quality_scores = []
        for i, line in enumerate(filehandle):
            if i % 4 == 3:
                quality_scores.append([ord(temp) - 33 for temp in line.strip()])

    # Calculate the mean quality score for each read
    mean_quality_scores = [sum(scores) / len(scores) for scores in quality_scores]

    # Apply rolling average to smoothen the plot
    rolling_mean = np.convolve(mean_quality_scores, np.ones(window_size) / window_size,
                               mode='valid')

    # Plot the rolling mean quality scores using matplotlib
    if len(rolling_mean) == 0:
        print("Error: No quality scores found.")
    else:
        plt.figure(figsize=(8, 6))
        plt.plot(rolling_mean, color='steelblue', linewidth=1.5)

        # Customize grid style
        plt.grid(color='lightgray', linestyle='--', linewidth=0.5, alpha=0.7)

        # Set axis labels and title with a bold font
        plt.xlabel('Read', fontsize=12, fontweight='bold')
        plt.ylabel('Mean Quality Score', fontsize=12, fontweight='bold')
        plt.title('MinION Quality Scores (Rolling Average)', fontsize=14, fontweight='bold')

        # Set axis tick font size and format x-axis tick labels with comma separator for thousands
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)

        # Customize plot border
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_linewidth(0.5)
        plt.gca().spines['left'].set_linewidth(0.5)

        # Save the plot as a high-resolution image
        plt.savefig('minion_quality_plot.png', dpi=300, bbox_inches='tight')

        plt.clf()


def calculate_sequence_complexity(seq):
    """
    Calculate the sequence complexity of a DNA sequence.

    Parameters:
        seq (str): A DNA sequence.

    Returns:
        float: The sequence complexity, which is the number of unique 4-mers
        in the sequence divided by the length of the sequence.
    """
    kmers = set()
    for i in range(len(seq) - 3):
        kmers.add(seq[i:i + 4])
    return len(kmers) / len(seq)


def plot_sequence_complexity(fastq_file):
    """
    Calculate and plot the sequence complexity of each sequence in a fastq file.

    Parameters:
        fastq_file (str): The path of a fastq file.

    Returns:
        None

    Note:
        The function saves the plot as an image in the current working directory named
        'sequence_complexity_<fastq_file>.png'.
        Additionally, the average sequence complexity is printed to the console.
    """
    # Create a new figure and axes for the plot
    plt.figure(figsize=(8, 6))

    complexities = []
    for seqrecord in SeqIO.parse(fastq_file, "fastq"):
        complexities.append(calculate_sequence_complexity(str(seqrecord.seq)))

    # Plot kernel density estimate (KDE) of sequence complexities
    sns.kdeplot(complexities, color="steelblue", linewidth=2, fill=True, alpha=0.6)

    # Add vertical line at the mean complexity
    mean_complexity = np.mean(complexities)
    plt.axvline(x=mean_complexity, color='red', linestyle='--', label='Mean Complexity')

    plt.title("KDE Plot: The Sequence Complexity Distribution", fontsize=14, fontweight='bold')
    plt.xlabel("Sequence Complexity", fontsize=12, fontweight='bold')
    plt.ylabel("Estimated Density", fontsize=12, fontweight='bold')
    plt.grid(color='lightgray', linestyle='--', linewidth=0.5, alpha=0.7)

    # Set the x-axis lower limit to 0
    plt.xlim(0, 1)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_linewidth(0.5)
    plt.gca().spines['left'].set_linewidth(0.5)
    # Add legend
    plt.legend()

    # Save plot to file and display it
    plt.savefig("sequence_complexity.png", dpi=300, bbox_inches='tight')

    # Clear the plot to avoid overlap with subsequent plots
    plt.clf()

    # Print average sequence complexity
    avg_complexity = np.mean(complexities)
    logging.info("The average sequence complexity is %.2f", avg_complexity)
