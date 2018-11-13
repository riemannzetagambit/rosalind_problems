from math import ceil
import sys

from .rosalind_utils import get_rosalind_data, process_fasta_file


def _get_read_overlap(seq1, seq2, min_overlap=0.5):
    '''
    Given two sequences that overlap by at least half the shorter one's length,
    return the sequence gotten by merging them on the overlap

    Example:
        seq1 = 'ACGG'
        seq2 = 'GTACG'
        returns 'GTACGG'
    '''
    short_seq, long_seq = sorted([seq1, seq2], key=lambda x: len(x))
    short_len, long_len = len(short_seq), len(long_seq)
    # e.g. if getting half the length of the shortest read, round up if it is odd in length
    min_overlap_len = int(ceil(min_overlap*short_len)) # have to cast to int, py3 ceil returns float. oddly

    # this range runs through the short sequence from longest to shortest length
    # start with longest match possible, which obeys principle of parsimony
    for i in range(0, short_len - min_overlap_len):
        # check if the end of the short sequence aligns with the beginning of the long sequence
        # i.e. overlap reads as (short -> long)
        if short_seq[i:] == long_seq[:short_len-i]:
            # attach long_seq to the end of short_seq but remove overlap
            # currently this is the shortest overlap since we've gone first
            return short_seq + long_seq[short_len-i:]

        elif long_seq[(long_len-short_len)+i:] == short_seq[:short_len-i]:
            return long_seq + short_seq[short_len-i:]
        
    # return none if there is no match
    return None


def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]
    '''
    reads = process_fasta_file(sequence_data).values()

    # use the first read as the seed for the chromosome and add on from there
    # in principle you could start with any read as a seed but randomizing this process could lead to variability in
    # performance and consistency
    chromosome = list(reads)[0]
    # keep track of what reads we've already put into the chromosome
    reads_assembled = [list(reads)[0]]

    # NOTE(dstone): the only assumptions we are allowed to make: "there exists a unique way to reconstruct the entire
    # chromosome from these reads by gluing together pairs of reads that overlap by more than half their length."

    # assemble the chromosome by finding the (unique, thank god for this problem we don't have to deal with bubbles)
    # pair that it has overlap. Start with the first read as the seed, find the other read it overlaps with, and set the
    # chromosome to those combined reads. Then repeat the process with the newly updated chromosome. Continue to add on
    # unique bits to the chromosome until you have exhausted all reads

    # we assume that all reads are used in assembly and do not deal with the case of leftover unaligned reads
    while len(reads_assembled) != len(reads):
        # we have not used all the reads in our assembly, so keep assembling, running over all the reads for matches to
        # the current chromosome update each time

        overlaps_by_read = {r: _get_read_overlap(chromosome, r) for r in reads if r not in reads_assembled}
        # get rid of reads with no overlap
        overlaps_by_read = {read: overlap for read, overlap in overlaps_by_read.items() if overlap is not None}
        shortest_overlap_read = min(overlaps_by_read, key=lambda x: len(overlaps_by_read[x]))
        # set the chromosome to the shortest overlap
        chromosome = _get_read_overlap(chromosome, shortest_overlap_read)
        reads_assembled.append(shortest_overlap_read)


    print(chromosome)
    return chromosome


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
