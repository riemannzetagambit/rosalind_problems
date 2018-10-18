from math import ceil
import sys

from rosalind_utils import get_reverse_complement, get_rosalind_data, process_fasta_file


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
    short_len = len(short_seq)
    # e.g. if getting half the length of the shortest read, round up if it is odd in length
    min_overlap_len = int(ceil(min_overlap*short_len)) # have to cast to int, py3 ceil returns float. oddly

    # this range runs through the short sequence from longest to shortest length
    # start with longest match possible, which obeys principle of parsimony
    for i in range(short_len - min_overlap_len, short_len + 1, 1):
        # check if the end of the short sequence aligns with the beginning of the long sequence
        # i.e. overlap reads as (short -> long)
        if short_seq[i:] == long_seq[:short_len-i]:
            # attach long_seq to the end of short_seq but remove overlap
            return short_seq + long_seq[short_len-i:]
    # now instead go through long sequence start with checking last len(short_seq) to coincide with short_seq and
    # iteratively check shorter sequences
    # again be as greedy as possible and look for as much overlap as possible
    for i in range(0, min_overlap_len + 1, 1):
        # do the opposite: check if the end of the long sequence aligns with the beginning of the short sequence
        # i.e. overlap reads as (long -> short)
        # if len(short_seq) = N, check last N nts of long_seq to coincide with all of short_seq, then N-1 nts with first
        # N-1 nts of short_seq, etc.
        if long_seq[short_len-min_overlap_len+i:] == short_seq[:min_overlap-i+1]:
            # attach short_seq to the end of long_seq but remove overlap
            return long_seq + short_seq[min_overlap-i+1:]
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
    chromosome = reads[0]
    # keep track of what reads we've already put into the chromosome
    reads_assembled = [reads[0]]

    # assemble the chromosome by finding the (unique, thank god for this problem we don't have to deal with bubbles)
    # pair that it has overlap. Start with the first read as the seed, find the other read it overlaps with, and set the
    # chromosome to those combined reads. Then repeat the process with the newly updated chromosome. Continue to add on
    # unique bits to the chromosome until you have exhausted all reads

    # we assume that all reads are used in assembly and do not deal with the case of leftover unaligned reads
    while len(reads_assembled) != len(reads):
        # we have not used all the reads in our assembly, so keep assembling, running over all the reads for matches to
        # the current chromosome update each time
        for read in reads:
            if read in reads_assembled:
                continue
            # this is the heavy lifting done to get overlap with this read (if any)
            overlap = _get_read_overlap(chromosome, read)
            if overlap is not None:
                # we have a new extension to the chromosome, so update it
                chromosome = overlap
                # make sure we don't use this read again
                reads_assembled.append(read)

    test_chromosome_solution = 'ATTAGACCTGCCGGAATAC'
    print(chromosome)
    assert chromosome == test_chromosome_solution
    return chromosome


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
