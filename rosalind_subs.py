import sys

from .rosalind_utils import get_rosalind_data

def _get_substring_positions(seq, motif):
    '''
    Find starting locations in seq where motif occurs
    '''
    motif_len = len(motif)

    # the i+1 is to match Rosalind conventions-- start from 1, not 0
    # this is not the fastest way to substring search
    return [i+1 for i in range(len(seq)) if seq[i:i + motif_len] == motif]


def solve_problem(sequence_data):
    '''
    Assumes only two sequences in file
    '''
    seq, motif = sequence_data[0], sequence_data[1]
    subs_positions = ' '.join(map(str, _get_substring_positions(seq, motif)))

    print(f'{subs_positions}')
    return subs_positions


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
