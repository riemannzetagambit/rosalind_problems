import sys

from .rosalind_utils import get_rosalind_data, process_fasta_file

def _get_kmers_from_sequence(seq, k):
    '''
    Given a sequence seq, find all k-mers within the sequence.

    e.g. for 
        seq = 'ATCG', 
        k = 3, 
        this would return
        ['ATC', 'TCG']
    '''
    # only do range through length - k + 1 to stay within string bounds
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]


def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]
    '''
    reads_dict = process_fasta_file(sequence_data)
    reads = reads_dict.values()
    # use the shortest read as the seed for the motif, since the longest common motif can be at most the length of the
    # shortest sequence
    starting_motif = min(reads, key=lambda x: len(x))
    # generate k-mers from the starting motif, starting with k = len(starting motif),
    # and generate successively smaller k-mers and check to see if each k-mer is in all other sequences
    # return on the first motif that is in all other sequences (which will be the longest, 
    # since we start from longest motif and go down
    for k in range(len(starting_motif), 0, -1):
        for kmer in _get_kmers_from_sequence(starting_motif, k):
            if all(kmer in seq for seq in reads):
                print(kmer)
                return kmer


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
