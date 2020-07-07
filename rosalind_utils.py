from collections import OrderedDict
from random import choice
import re
from typing import Generator, List, TextIO

AA_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
DNA_ALPHABET = 'ACGT'
RNA_ALPHABET = 'ACGU'

RC_DICT = {'A': 'T',
           'T': 'A',
           'G': 'C',
           'C': 'G',
           'N': 'N'}

RNA_CODON_DICT = {
                  'UUU': 'F',
                  'UUC': 'F',
                  'UUA': 'L',
                  'UUG': 'L',
                  'UCU': 'S',
                  'UCC': 'S',
                  'UCA': 'S',
                  'UCG': 'S',
                  'UAU': 'Y',
                  'UAC': 'Y',
                  'UAA': 'Stop',
                  'UAG': 'Stop',
                  'UGU': 'C',
                  'UGC': 'C',
                  'UGA': 'Stop',
                  'UGG': 'W',
                  'CUU': 'L',
                  'CUC': 'L',
                  'CUA': 'L',
                  'CUG': 'L',
                  'CCU': 'P',
                  'CCC': 'P',
                  'CCA': 'P',
                  'CCG': 'P',
                  'CAU': 'H',
                  'CAC': 'H',
                  'CAA': 'Q',
                  'CAG': 'Q',
                  'CGU': 'R',
                  'CGC': 'R',
                  'CGA': 'R',
                  'CGG': 'R',
                  'AUU': 'I',
                  'AUC': 'I',
                  'AUA': 'I',
                  'AUG': 'M',
                  'ACU': 'T',
                  'ACC': 'T',
                  'ACA': 'T',
                  'ACG': 'T',
                  'AAU': 'N',
                  'AAC': 'N',
                  'AAA': 'K',
                  'AAG': 'K',
                  'AGU': 'S',
                  'AGC': 'S',
                  'AGA': 'R',
                  'AGG': 'R',
                  'GUU': 'V',
                  'GUC': 'V',
                  'GUA': 'V',
                  'GUG': 'V',
                  'GCU': 'A',
                  'GCC': 'A',
                  'GCA': 'A',
                  'GCG': 'A',
                  'GAU': 'D',
                  'GAC': 'D',
                  'GAA': 'E',
                  'GAG': 'E',
                  'GGU': 'G',
                  'GGC': 'G',
                  'GGA': 'G',
                  'GGG': 'G',
                 }

def get_rosalind_data(filename):
    with open(filename, 'r') as f:
        # dangerous for v. large files
        rosalind_data = f.readlines()

    return [rd.strip() for rd in rosalind_data]


def process_fasta_file(file_handle: TextIO, read_id_delimiter: str='>') -> dict:
    '''
    Given a file in fasta format, return the read IDs and their associated reads (which can span over multiple lines) as
    a dictionary
    This method assumes that the read IDs in the fasta file are unique, that there are read IDs followed immediately but
    an unknown number of lines of sequence data, and no other gaps (i.e. no blank lines). (this is the standard fasta
    format)
    :param file_handle: list
        List of the lines in the file, usually retrieved via, e.g., open('path/to/file') as f: return f.readlines()
    :return: dict
        Dictionary with keys the read ID (assumed to be unique), values the associated reads (i.e. DNA sequences)
    '''
    read_dict = OrderedDict()
    read_id = None
    #lines = iter(lines)
    for line in file_handle:
        read = ''
        # iterate through lines to find identifier or set of lines corresponding to reads
        while True:
            # readline() returns '' for final line (dangerous, I shouldn't be using a while loop, yada yada I'll rewrite this later)
            # so catch that here
            if line is None or line == '':
                read_dict[read_id] = read
                break
            elif line.startswith(read_id_delimiter):
                if read_id is not None:
                    # we already have a read ID from last read ID assignment and have been constructing read
                    read_dict[read_id] = read
                # if read_id is None, this is the first line, so just assign the read ID,
                read_id = line.lstrip(read_id_delimiter).strip()
                break
            else:
                # construct the read by iterating through until we hit None (EOF) or next read ID ('>blahblah')
                read += line.strip()
                #line = next(lines, None)
                line = file_handle.readline()
    return dict(read_dict)


def random_genetic_sequence(length=150, alphabet=DNA_ALPHABET):
    return ''.join([choice(alphabet) for _ in range(length)])


def get_reverse_complement(seq):
    return ''.join([RC_DICT[x] for x in seq[::-1]])


# have used in rosalind_hamm and rosalind_corr problems at the least, so putting here
def get_hamming_distance(seq1, seq2):
    '''
    Given two sequences of equal length, return the Hamming distance, i.e. the number of differences between the two
    sequences (computed position by position)

    Raises if sequences are not of equal length
    '''
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must of equal length. '
                         'Instead got\n\tlen(seq1): {}\n\tlen(seq2): {}'.format(len(seq1), len(seq2)))
    return sum(1 for i in range(len(seq1)) if seq1[i] != seq2[i])


def get_reading_frames(seq: str,
                       include_reverse_complement: bool = True,
                       include_offsets: bool = True) -> Generator[str, None, None]:
    '''
    For a given sequence, return the 6 possible reading frames associated to it and its reverse complement,
    without regard for start/stop codons.

    :param seq: str
        Input dna sequence as a string
        Example: 'ATCGCAGATCGA'
    :param include_reverse_complement: bool
        Whether or not to return values that include the reverse complement in the solutions.
        For example, in some instances, seq represents a sense strand, and we do not want to look at
        translations/transcriptions of the reverse complement
    :param include_offsets: bool
        Whether or not to include the three sets of frames that translation could occur at, one for each starting
        location within a possible length-3 codon
        For example, in some instances, seq represents a sense strand where transcription is already occurring-- in
        this case we do not need to account for the other two possibilities where translation can begin

    :return: generator of string reading frames
        a generator that is the (up to) six possible reading frames, as DNA

    These reading frames are translated into AA, so the return result is an array of 6 arrays of AA translations for
    each possible translation
    '''
    if include_reverse_complement:
        # get the reverse complement and include it
        rc_seq = get_reverse_complement(seq)
        # convert to RNA for translation
        seqs = (convert_dna_to_rna(seq), convert_dna_to_rna(rc_seq))
    else:
        seqs = (convert_dna_to_rna(seq), )
    if include_offsets:
        offsets = 3
    else:
        offsets = 1

    rfs = []
    # go through 2x forward and reverse positions
    for tmp_seq in seqs:
        seq_len = len(tmp_seq)
        # go through 3x places to start for unique codons
        for i in range(offsets):
            # iterate through the sequence from the starting point in steps of 3 until you hit last translatable codon
            yield tmp_seq[i: seq_len - (seq_len - i) % 3]


def get_translated_reads_frames(seq: str,
                                include_reverse_complement: bool = True,
                                include_offsets: bool = True) -> Generator[str, None, None]:
    '''
    Given an input sequence, return each of the 6 possible reading frames (offset up to 3 x reverse complement),
    translated into amino acids

    :param seq: str
        Input dna sequence as a string
        Example: 'ATCGCAGATCGA'
    :param include_reverse_complement: bool
        Whether or not to return values that include the reverse complement in the solutions.
        For example, in some instances, seq represents a sense strand, and we do not want to look at
        translations/transcriptions of the reverse complement
    :param include_offsets: bool
        Whether or not to include the three sets of frames that translation could occur at, one for each starting
        location within a possible length-3 codon
        For example, in some instances, seq represents a sense strand where transcription is already occurring-- in
        this case we do not need to account for the other two possibilities where translation can begin

    :return: generator of string translated reading frames
        a generator that is the (up to) six possible reading frames, translated to amino acids
    '''
    rfs = get_reading_frames(seq,
                             include_reverse_complement=include_reverse_complement,
                             include_offsets=include_offsets)
    # _get_reading_frames already shifted the reading from for us, so we just iterate through each rf the same way
    # translate the reading frame with the RNA_CODON_DICT (rcd)
    # yield a generator of a strings that are the translating reading frames
    for rf in rfs:
        yield ''.join(RNA_CODON_DICT[rf[j: j+3]] for j in range(0, len(rf), 3))


def get_open_reading_frames(seq: str,
                            include_reverse_complement: bool = True,
                            include_offsets: bool = True,
                            include_overlapping_solutions: bool = True) -> Generator[str, None, None]:
    '''
    :param seq: str
        Input dna sequence as a string
        Example: 'ATCGCAGATCGA'
    :param include_reverse_complement: bool
        Whether or not to return values that include the reverse complement in the solutions.
        For example, in some instances, seq represents a sense strand, and we do not want to look at
        translations/transcriptions of the reverse complement
    :param include_offsets: bool
        Whether or not to include the three sets of frames that translation could occur at, one for each starting
        location within a possible length-3 codon
        For example, in some instances, seq represents a sense strand where transcription is already occurring-- in
        this case we do not need to account for the other two possibilities where translation can begin
    :param include_overlapping_solutions: bool
        Whether to include multiple solutions if an ORF has multiple start points before a stop point
        Example: If True, if an ORF is MVYIADKQHVASREAYGHMFKVCA, then so is MFKVCA.
                 If False, only includes the longest match

    :return: generator of string ORFs
        a generator that is the unique ORFs associated with the sequence
    '''
    if include_overlapping_solutions:
        orf_regex = re.compile('(?=(?P<orf>M.*?Stop))')
    else:
        orf_regex = re.compile('(?P<orf>M.*?Stop)')
    trfs = get_translated_reads_frames(seq,
                                       include_reverse_complement=include_reverse_complement,
                                       include_offsets=include_offsets)
    orfs = []
    for trf in trfs:
        orfs.extend([m.group('orf').replace('Stop', '') for m in orf_regex.finditer(trf)])
    # casting to set() gives only unique values
    for orf in set(orfs):
        yield orf


def convert_dna_to_rna(seq: str) -> str:
    return seq.replace('T', 'U')


def convert_dna_to_aa(seq: str, 
                      include_reverse_complement: bool = True,
                      include_offsets: bool = True) -> Generator[str, None, None]:
    '''
    Convenience wrapper that is another name for get_translated_reads_frames
    '''
    return get_translated_reads_frames(seq,
                                       include_reverse_complement=include_reverse_complement,
                                       include_offsets=include_offsets)


def convert_dna_to_protein(seq: str,
                           include_reverse_complement: bool = True,
                           include_offsets: bool = True,
                           include_overlapping_solutions: bool = True) -> Generator[str, None, None]:
    '''
    Convenience wrapper that is another name for get_open_reading_frames
    '''
    return get_open_reading_frames(seq,
                                   include_reverse_complement=include_reverse_complement,
                                   include_offsets=include_offsets,
                                   include_overlapping_solutions=include_overlapping_solutions)
