from collections import Counter
import sys

from rosalind_utils import get_rosalind_data, RNA_CODON_DICT


def solve_problem(sequence_data):
    '''
    Assumes a single sequence
    '''
    seq = sequence_data[0]
    codon_len = 3
    # chunk the sequence up into 3 nts at a time; chops off anything less than 3 at the end
    coding_seqs = [seq[pos:pos + codon_len] for pos in range(0, len(seq), codon_len)]
    # need to remove the final stop codon, since it doesn't actually code for anything
    prot_seq = ''.join([RNA_CODON_DICT[coding_seq] for coding_seq in coding_seqs]).rstrip('Stop')

    print(f'{prot_seq}')
    return prot_seq


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
