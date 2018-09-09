import sys

from rosalind_utils import get_rosalind_data

def _get_hamming_distance(seq1, seq2):
    '''
    Given two sequences of equal length, return the Hamming distance, i.e. the number of differences between the two
    sequences (computed position by position)

    Raises if sequences are not of equal length
    '''
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must of equal length. '
                         'Instead got\n\tlen(seq1): {}\n\tlen(seq2): {}'.format(len(seq1), len(seq2)))
    return sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])


def solve_problem(sequence_data):
    '''
    Assumes only two sequences in file
    '''
    seq1, seq2 = sequence_data[0], sequence_data[1]
    hamming_distance = _get_hamming_distance(seq1, seq2)

    print(f'{hamming_distance}')
    return hamming_distance


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
