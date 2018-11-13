from itertools import combinations
import sys

from .rosalind_utils import get_rosalind_data, process_fasta_file


def _has_k_overlap(string1, string2, k):
    return string1[-k:] == string2[:k]


def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]

    '''
    adjacency_list = []
    k = 3 # want to generate the adjacency list for O_3, so set overlap to k = 3
    read_dict = process_fasta_file(sequence_data)

    for read_pair in combinations(read_dict.keys(), 2):
        if _has_k_overlap(read_dict[read_pair[0]], read_dict[read_pair[1]], k):
            adjacency_list.append(read_pair)
        # check if the other direction gives you a match
        elif _has_k_overlap(read_dict[read_pair[1]], read_dict[read_pair[0]], k):
            adjacency_list.append(read_pair[::-1]) # since we matched in opposite order, reverse the tuple to append

    print('\n'.join([' '.join(pair) for pair in adjacency_list]))
    return adjacency_list


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
