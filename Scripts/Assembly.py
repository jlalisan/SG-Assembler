from Bio import SeqIO
from collections import defaultdict


def read_fastq(filename):
    sequences = []
    for record in SeqIO.parse(filename, "fastq"):
        sequences.append(str(record.seq))
    return sequences


def create_kmers(sequence, k):
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]


def build_abruijn_graph(sequences, k):
    graph = defaultdict(set)
    for sequence in sequences:
        kmers = create_kmers(sequence, k)
        for kmer in kmers:
            prefix = kmer[:-1]
            suffix = kmer[1:]
            graph[prefix].add(suffix)
    return graph


def generate_disjointigs(graph):
    disjointigs = []
    graph_items = list(graph.items())  # Create a copy of the dictionary items

    for node, edges in graph_items:
        if len(edges) == 1:
            next_node = next(iter(edges))
            if len(graph[next_node]) == 1:
                disjointig = node + next_node[-1]
                disjointigs.append(disjointig)
                if node in graph:
                    del graph[node]
                if next_node in graph:
                    del graph[next_node]

                # Extend the disjointig
                while next_node in graph and len(graph[next_node]) == 1:
                    next_node = next(iter(graph[next_node]))
                    disjointig += next_node[-1]
                    del graph[next_node]

                disjointigs.append(disjointig)

    return disjointigs


def concatenate_disjointigs(disjointigs):
    concatenated_string = ''.join(disjointigs)
    return concatenated_string


def count_kmers(sequence, k):
    kmers = { }
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers[kmer] = kmers.get(kmer, 0) + 1
    return kmers


def identify_repetitive_sequences(concatenated_disjointigs, k, threshold):
    kmers = count_kmers(concatenated_disjointigs, k)
    repetitive_kmers = {kmer: count for kmer, count in kmers.items() if count >= threshold}
    return repetitive_kmers

def resolve_bridged_repeats(repeat_graph, long_reads, k):
    for read in long_reads:
        kmers = list(create_kmers(read, k))
        for i in range(len(kmers) - 1):
            prefix = kmers[i][:-1]
            suffix = kmers[i + 1][1:]
            if prefix in repeat_graph and suffix in repeat_graph[prefix]:
                repeat_graph[prefix].remove(suffix)
                repeat_graph[prefix].add(kmers[i + 1])


def create_repeat_graph(concatenated_disjointigs, k, threshold, long_reads):
    repetitive_kmers = identify_repetitive_sequences(concatenated_disjointigs, k, threshold)
    repeat_graph = defaultdict(set)
    for kmer, count in repetitive_kmers.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        repeat_graph[prefix].add(suffix)
    resolve_bridged_repeats(repeat_graph, long_reads, k)
    return repeat_graph

def traverse_repeat_graph(repeat_graph, long_reads):
    assembled_contigs = []
    contig_id = 1
    for prefix, suffixes in repeat_graph.items():
        for suffix in suffixes:
            for read in long_reads:
                if prefix in read and suffix in read:
                    start = read.index(prefix)
                    end = read.index(suffix) + len(suffix)
                    contig = read[start:end]
                    assembled_contigs.append((f"contig_{contig_id}", contig))
                    contig_id += 1
                    break
    return assembled_contigs

def write_fasta(contigs, output_file):
    with open(output_file, "w") as f:
        for contig_id, sequence in contigs:
            f.write(f">{contig_id}\n")
            f.write(f"{sequence}\n")

def main():
    filename = "test.fq"
    k = 31
    threshold = 12
    output_file = 'mytiggies.fasta'

    sequences = read_fastq(filename)
    abruijn_graph = build_abruijn_graph(sequences, k)
    disjointigs = generate_disjointigs(abruijn_graph)
    concatenated_disjointigs = concatenate_disjointigs(disjointigs)
    repeat_graph = create_repeat_graph(concatenated_disjointigs, k, threshold, sequences)
    assembled_contigs = traverse_repeat_graph(repeat_graph, sequences)
    write_fasta(assembled_contigs, output_file)


if __name__ == "__main__":
    main()
