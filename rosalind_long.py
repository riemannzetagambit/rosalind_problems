from math import ceil
import sys

from rosalind_utils import get_reverse_complement, get_rosalind_data, process_fasta_file


def _has_sufficient_overlap(seq1, seq2, min_overlap=0.5):
    '''
    Return true if at least half the shortest sequence out of seq[12] overlaps with the other sequence
    Otherwise return false
    '''
    short_seq, long_seq = sorted([seq1, seq2], key=lambda x: len(x)) 
    short_len = len(short)
    # look for overlap from both ends
    # we check for at least half overlap from the shortest read
    # first `any` checks for longest possible overlap (since by assumption of parsimony we want to find greatest overlap)
    # by checking for overlap from full length of short sequence down to half its length at its end, looking for overlap
    # with the beginning of the long sequence
    # second `any` swaps the two, checks for overlap of end of the long sequence with the beginning of the short sequence


    # TODO(dstone): instead make this function return the i for which there is overlap between the two sequences (along
    # with order), or None if no overlap
    if (any([short_seq[i:] == long_seq[:short_len-i] for i in range(0, short_len - min_overlap_len + 1, 1)]) 
        or any([long_seq[i:] == short_seq[:short_len-i] for i in range(0, short_len - min_overlap_len + 1, 1)])):
        return True
    else:
        return False

    for i in range(0, short_len - min_overlap_len + 1, 1):
        if blah:
            return i, 1
        elif blah:
            return i, -1 # -1 indicates short, long
    return None, NoneV

def _get_read_overlap(seq1, seq2, min_overlap):
    '''
    Given two sequences that overlap by at least half the shorter one's length, 
    return the sequence gotten by merging them on the overlap

    Example:
        seq1 = 'ACGG'
        seq2 = 'GTACG'
        returns 'GTACGG'
    '''
    short_seq, long_seq = sorted([seq1, seq2], key=lambda x: len(x)) 
    short_len = len(short)
    min_overlap_len = ceil(min_overlap*short_len)

    for i in range(0, short_len - min_overlap_len + 1, 1):
        if short_seq[i:] == long_seq[:short_len-i]:
            return short_seq + long_seq[short_len-i:]
        elif long_seq[i:] == short_seq[:short_len-i]:
        elif blah:
            return 
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
    chromosome = reads[0]
    # keep track of what reads we've already put into the chromosome
    reads_assembled = []


    print('\n'.join([' '.join(map(str, tup)) for tup in palindromes_loc_lens]))
    return palindromes_loc_lens


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
