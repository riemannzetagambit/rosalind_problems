from typing import List
import sys

from rosalind_utils import get_rosalind_data, process_fasta_file

_TRANSITIONS = [('A', 'G'), ('C', 'T')]
_TRANSVERSIONS = [('A', 'C'), ('A', 'T'), ('C', 'G'), ('G', 'T')]


def solve_problem(sequence_data: list) -> str:
    '''
    Assumptions: input file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    ...
    '''
    reads_dict = process_fasta_file(sequence_data)
    reads = list(reads_dict.values())

    num_reads, len_s1, len_s2 = len(reads), len(reads[0]), len(reads[1])
    if num_reads != 2 or len_s1 != len_s2:
        raise ValueError('Expected only two reads/sequences to process, which must be of equal length. Got {} reads, '
                         'with first two having lengths {} and {}.'.format(num_reads, len_s1, len_s2))
    transitions, transversions = 0, 0
    for pair in zip(*reads):
        if pair[0] != pair[1]:
            # we have a mismatch at this element-- sort it so we can check if an element in transitions/versions
            sp = tuple(sorted(pair))
            if sp in _TRANSITIONS:
                transitions += 1
            elif sp in _TRANSVERSIONS:
                transversions += 1
            else:
                raise ValueError('Got mismatch that was not in DNA transitions/transversions: {}'.format(sp))

    print(transitions/transversions)
    return transitions/transversions


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
