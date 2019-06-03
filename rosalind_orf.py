import re
import sys

from rosalind_utils import get_open_reading_frames, get_reverse_complement, get_rosalind_data, process_fasta_file
from rosalind_utils import RNA_CODON_DICT as rcd

def _get_open_reading_frames(seq):
    '''
    :param seq: str
        DNA string to determine (unique) open reading frames from. A single DNA sequence can have multiple ORFs, so we
        return a unique set of all possible ORFs generated from seq

    :return: set
        set() of unique ORFs from seq
    '''
    # this function was ported over to rosalind_utils for more general usage in other problems. This function is a no-op
    # left for posterity's sake
    

def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]

    '''
    # we expect a single read, which is the first value in the process_fasta_file returned dict
    read = list(process_fasta_file(sequence_data).values())[0]

    # because we use the generator twice below, we need to cast this to a list lest we exhaust it
    orfs = list(get_open_reading_frames(read))

    print('\n'.join(orfs))
    return '\n'.join(orfs)


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
