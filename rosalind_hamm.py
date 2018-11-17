import sys

from rosalind_utils import get_hamming_distance, get_rosalind_data

def solve_problem(sequence_data):
    '''
    Assumes only two sequences in file
    '''
    seq1, seq2 = sequence_data[0], sequence_data[1]
    hamming_distance = get_hamming_distance(seq1, seq2)

    print(f'{hamming_distance}')
    return hamming_distance


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
