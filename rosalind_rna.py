import sys

from rosalind_utils import DNA_ALPHABET, get_rosalind_data

def solve_problem(sequence_data):
    dna_sequence = sequence_data[0]
    rna_sequence = dna_sequence.replace('T', 'U')
    print(rna_sequence)
    return rna_sequence


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
