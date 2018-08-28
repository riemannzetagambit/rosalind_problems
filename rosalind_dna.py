from collections import Counter
import sys

from .rosalind_utils import DNA_ALPHABET

def solve_problem(sequence_data):
   nt_counts = dict(Counter(random_sequence).most_common()) 
   print(' '.join([nt_counts[nt] for nt in DNA_ALPHABET]))


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    data_file = args[-1]
    with open(data_file, 'w') as f:
        sequence_data = f.readlines()
    solve_problem(sequence_data)
