import sys

from .rosalind_utils import get_reverse_complement, get_rosalind_data, process_fasta_file


def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]
    '''
    seq = process_fasta_file(sequence_data).values()[0]
    seq_len = len(seq)

    # problem asks us to restrict to between 4 and 12
    palindrome_min_length = 4
    palindrome_max_length = 12

    palindromes_loc_lens = []
    for k in range(palindrome_min_length, palindrome_max_length+1):
        for i in range(seq_len-k+1):
            if seq[i:i+k] == get_reverse_complement(seq[i:i+k]):
                palindromes_loc_lens.append((i+1, k)) # positions start at 1, not 0, for Rosalind
    print('\n'.join([' '.join(map(str, tup)) for tup in palindromes_loc_lens]))
    return palindromes_loc_lens


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
