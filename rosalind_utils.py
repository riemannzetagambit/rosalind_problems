from random import choice

DNA_ALPHABET = 'ACGT'
RNA_ALPHABET = 'ACGU'

RC_DICT = {'A': 'T',
           'T': 'A',
           'G': 'C',
           'C': 'G',
           'N': 'N'}

RNA_CODON_DICT = {
                  'UUU': 'F',
                  'UUC': 'F',
                  'UUA': 'L',
                  'UUG': 'L',
                  'UCU': 'S',
                  'UCC': 'S',
                  'UCA': 'S',
                  'UCG': 'S',
                  'UAU': 'Y',
                  'UAC': 'Y',
                  'UAA': 'Stop',
                  'UAG': 'Stop',
                  'UGU': 'C',
                  'UGC': 'C',
                  'UGA': 'Stop',
                  'UGG': 'W',
                  'CUU': 'L',
                  'CUC': 'L',
                  'CUA': 'L',
                  'CUG': 'L',
                  'CCU': 'P',
                  'CCC': 'P',
                  'CCA': 'P',
                  'CCG': 'P',
                  'CAU': 'H',
                  'CAC': 'H',
                  'CAA': 'Q',
                  'CAG': 'Q',
                  'CGU': 'R',
                  'CGC': 'R',
                  'CGA': 'R',
                  'CGG': 'R',
                  'AUU': 'I',
                  'AUC': 'I',
                  'AUA': 'I',
                  'AUG': 'M',
                  'ACU': 'T',
                  'ACC': 'T',
                  'ACA': 'T',
                  'ACG': 'T',
                  'AAU': 'N',
                  'AAC': 'N',
                  'AAA': 'K',
                  'AAG': 'K',
                  'AGU': 'S',
                  'AGC': 'S',
                  'AGA': 'R',
                  'AGG': 'R',
                  'GUU': 'V',
                  'GUC': 'V',
                  'GUA': 'V',
                  'GUG': 'V',
                  'GCU': 'A',
                  'GCC': 'A',
                  'GCA': 'A',
                  'GCG': 'A',
                  'GAU': 'D',
                  'GAC': 'D',
                  'GAA': 'E',
                  'GAG': 'E',
                  'GGU': 'G',
                  'GGC': 'G',
                  'GGA': 'G',
                  'GGG': 'G',
                 }

def get_rosalind_data(filename):
    with open(filename, 'r') as f:
        # dangerous for v. large files
        rosalind_data = f.readlines()

    return [rd.strip() for rd in rosalind_data]


def random_genetic_sequence(length=150, alphabet=DNA_ALPHABET):
    return ''.join([choice(alphabet) for _ in range(length)])


def get_reverse_complement(seq):
    return ''.join([RC_DICT[x] for x in seq[::-1]])
