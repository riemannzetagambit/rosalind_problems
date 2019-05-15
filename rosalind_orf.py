import re
import sys

from rosalind_utils import get_reverse_complement, get_rosalind_data, process_fasta_file
from rosalind_utils import RNA_CODON_DICT as rcd

def _get_reading_frames(seq):
    '''
    For a given sequence, return the 6 possible reading frames associated to it and its reverse complement,
    without regard for start/stop codons.

    These reading frames are translated into AA, so the return result is an array of 6 arrays of AA translations for
    each possible translation
    '''
    rc_seq = get_reverse_complement(seq)
    # convert to RNA for translation
    seq, rc_seq = seq.replace('T', 'U'), rc_seq.replace('T', 'U')

    rfs = []
    # go through 2x forward and reverse positions
    for tmp_seq in [seq, rc_seq]:
        seq_len = len(tmp_seq)
        # go through 3x places to start for unique codons
        for i in range(3):
            # iterate through the sequence from the starting point in steps of 3 until you hit last translatable codon
            rfs.append(tmp_seq[i: seq_len - (seq_len - i) % 3])

    # return an array with 6 contents, each itself an array with possible codon translation of that reading frame
    return rfs


def _get_translated_reads_frames(seq):
    '''
    Break out this function from _get_reading_frames and just do the translation here
    '''
    rfs = _get_reading_frames(seq)
    # _get_reading_frames already shifted the reading from for us, so we just iterate through each rf the same way
    # translate the reading frame with the RNA_CODON_DICT (rcd)
    # yield a generator of a generator
    for rf in rfs:
        # yield another generator that is the amino acids, proceeding through nucleotides 3 at a time
        yield (rcd[rf[j: j+3]] for j in range(0, len(rf), 3))

def _get_open_reading_frames(seq):
    '''
    :param seq: str
        DNA string to determine (unique) open reading frames from. A single DNA sequence can have multiple ORFs, so we
        return a unique set of all possible ORFs generated from seq

    :return: set
        set() of unique ORFs from seq
    '''
    # TODO: just use regex instead: scratch code from ipython session that is simpler than the dropwhile and takewhile
    # hoops I was jumping through. much simpler
    trfs = _get_translated_reads_frames(seq)
    orfs = []
    for trf in trfs:
        orfs.extend([m.group(1).replace('Stop', '') for m in re.finditer('(?=(M.*?Stop))', ''.join(trf))])
    # casting to set() gives only unique values
    return set(orfs)
    

def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]

    '''
    # we expect a single read, which is the first value in the process_fasta_file returned dict
    read = list(process_fasta_file(sequence_data).values())[0]

    orfs = _get_open_reading_frames(read)

    print('\n'.join(orfs))
    return '\n'.join(orfs)


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
