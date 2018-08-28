from collections import Counter
import sys

from rosalind_utils import DNA_ALPHABET, get_rosalind_data

def solve_problem(sequence_data):
    sequence = sequence_data[0]
    nt_counts = dict(Counter(sequence).most_common()) 
    nt_counts_str = ' '.join([str(nt_counts[nt]) for nt in DNA_ALPHABET])
    print(nt_counts_str)
    return nt_counts_str


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
