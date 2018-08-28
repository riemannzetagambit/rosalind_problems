from random import choice

DNA_ALPHABET = 'ACGT'
RNA_ALPHABET = 'ACGU'

RC_DICT = {'A': 'T',
           'T': 'A',
           'G': 'C',
           'C': 'G',
           'N': 'N'}


def get_rosalind_data(filename):
    with open(filename, 'r') as f:
        rosalind_data = f.readlines()

    return rosalind_data


def random_genetic_sequence(length=150, alphabet=DNA_ALPHABET):
    return ''.join([choice(alphabet) for _ in range(length)])


def get_reverse_complement(seq):
    return ''.join([RC_DICT[x] for x in seq[::-1]])
