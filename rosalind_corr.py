from collections import Counter
from itertools import product
import sys

from .rosalind_utils import (get_hamming_distance, get_reverse_complement,
                            get_rosalind_data, process_fasta_file)


def _get_read_ordered_by_rc(read):
    '''
    Given a read, return it or its reverse complement, whichever comes first when sorted alphabetically (with `sored`)
    Allows me to unambiugously resolve equivalences when decided if two reads are the same up to reverse complement
    '''
    sorted_read = sorted([read, get_reverse_complement(read)])
    return sorted_read[0]

def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]

    if correctly sequenced, then appears in the dataset at least twice (possibly as a reverse complement)
    if incorrectly sequenced, then appears in the dataset exactly once, and its Hamming distance is 1 with respect to
        exactly one correct read in the dataset (or its reverse complement)
    '''
    reads = process_fasta_file(sequence_data).values()

    # error reads are unique, and below we use the _get_read_ordered_by_rc to unambiguously sort the collection of reads
    # we need a way to track the original error reads so as to identify which output to display, since rc equivalence
    # otherwise adds ambiguity; note this dict is useless for correct reads (which can include rc-equivalent sequences)
    rc_ordering_dict = {_get_read_ordered_by_rc(r): r for r in reads}
    # need to get all reads that show up at least twice, but need to account for reverse complement
    ordered_reads = Counter(map(_get_read_ordered_by_rc, reads))
    # correct reads have at least two appearances (modulo reverse complement) in the fasta reads
    correct_reads = [read for read, count in ordered_reads.items() if count >= 2]
    # revert error reads to the form found in the data if they were rc'ed in ordered_reads
    error_reads = [rc_ordering_dict[read] for read, count in ordered_reads.items() if count == 1]

    corrected_reads = []
    for (error_read, correct_read) in product(error_reads, correct_reads):
        # need to iterate through and check if distances between error and correct reads are 1, up to RC
        # note the 'ground truth' in this case will be determined by the error read, since it is unique
        # which is to say, when in doubt of whether to return the reverse complement or not, 
        # choose the one that is Hamming distance of 1 away from unique error read in the set
        if get_hamming_distance(error_read, correct_read) == 1:
            corrected_reads.append((error_read, correct_read))
        elif get_hamming_distance(error_read, get_reverse_complement(correct_read)) == 1:
            corrected_reads.append((error_read, get_reverse_complement(correct_read)))

    # catch errors
    for er in error_reads:
        if all([er != x and get_reverse_complement(er) != x for (x, _) in corrected_reads]):
            raise ValueError('No paired correct read found for error read:\n\t{}'.format(er))
    for cr in correct_reads:
        if all([cr != x and get_reverse_complement(cr) != x for (_, x) in corrected_reads]):
            raise ValueError('No paired error read found for correct read:\n\t{}'.format(cr))

    print('\n'.join(['{}->{}'.format(er, cr) for (er, cr) in corrected_reads]))
    # from now on I want to return the string solution that Rosalind expects, rather than my custom data structure used
    # along the way. Instead join the structure into the string format that Rosalind would agree with-- makes unit
    # testing easier
    return '\n'.join(['{}->{}'.format(er, cr) for (er, cr) in corrected_reads])


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
