import sys

from .rosalind_utils import get_rosalind_data, RC_DICT

def solve_problem(sequence_data):
    dna_seq = sequence_data[0]
    # reverse the string
    reverse_seq = dna_seq[::-1].strip()
    # complement it
    rc_seq = ''.join([RC_DICT[x] for x in reverse_seq])
    print(rc_seq)


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
