from collections import Counter
import sys

from rosalind_utils import DNA_ALPHABET, get_rosalind_data, process_fasta_file


def _get_consensus_profile_matrix(sequence_data):
    '''

    '''
    position_counter_dict = {}
    profile_matrix_dict = {}
    read_dict = process_fasta_file(sequence_data)
    sequences = list(read_dict.values()) # we don't care about read IDs in this problem

    # assume we have a single length for all sequences
    seq_len = len(sequences[0])
    if not all([len(seq) == seq_len for seq in sequences]):
        raise ValueError("Assumption violated that all sequences are of equal length! Some sequences are not same len")
    # for each position, get the counts for everything in that position
    for i in range(seq_len):
        position_counter_dict[i] = {}
        nts_at_pos = ''.join([seq[i] for seq in sequences])
        nt_dict = dict(Counter(nts_at_pos).most_common())
        # need to account for zeroes
        position_counter_dict[i] = {nt: nt_dict.get(nt, 0) for nt in DNA_ALPHABET}

    # now convert position-based dictionary to nucleotide-based dict
    for nt in DNA_ALPHABET:
        profile_matrix_dict[nt] = [position_counter_dict[i][nt] for i in range(seq_len)]

    return position_counter_dict, profile_matrix_dict



def _get_consensus_sequence(sequence_data):
    '''

    '''
    profile_matrix_dict, _ = _get_consensus_profile_matrix(sequence_data)
    positions = sorted(map(int, profile_matrix_dict.keys())) # the keys are positions, so we get their ordered int values

    # for each position, find the nucleotdie with the max number of counts and use that as the seq nt in that position
    return ''.join([max(profile_matrix_dict[i], key=profile_matrix_dict[i].get) for i in positions])

def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]

    '''
    _, profile_matrix_dict = _get_consensus_profile_matrix(sequence_data)
    # convert dict to printable string
    profile_matrix = '\n'.join(['{nt}: {counts}'.format(nt=nt, counts=' '.join(map(str, counts))) 
                                for nt, counts in profile_matrix_dict.items()])
    consensus_seq = _get_consensus_sequence(sequence_data)

    print(f'{consensus_seq}\n{profile_matrix}')
    return consensus_seq, profile_matrix


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
