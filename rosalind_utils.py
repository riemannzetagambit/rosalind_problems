from random import choice

DNA_ALPHABET = 'ACGT'
RNA_ALPHABET = 'ACGU'

def random_genetic_sequence(length=150, alphabet=DNA_ALPHABET):
    return ''.join([choice(alphabet) for _ in range(length)])
