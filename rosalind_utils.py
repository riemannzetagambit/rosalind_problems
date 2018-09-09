from random import choice

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


def process_fasta_file(lines):
    '''
    Given a file in fasta format, return the read IDs and their associated reads (which can span over multiple lines) as
    a dictionary

    This method assumes that the read IDs in the fasta file are unique, that there are read IDs followed immediately but
    an unknown number of lines of sequence data, and no other gaps (i.e. no blank lines). (this is the standard fasta
    format)

    :param lines: list
        List of the lines in the file, usually retrieved via, e.g., open('path/to/file') as f: return f.readlines()

    :return: dict
        Dictionary with keys the read ID (assumed to be unique), values the associated reads (i.e. DNA sequences)
    '''
    read_dict = {}
    read_id = None
    lines = iter(lines)
    for line in lines:
        read = ''
        # iterate through lines to find identifier or set of lines corresponding to reads
        while True:
            if line is None:
                # at end of file, so assign read GC content to last known read ID
                read_dict[read_id] = read
                break
            elif line.startswith('>'):
                if read_id is not None:
                    # we already have a read ID from last read ID assignment and have been constructing read
                    # so just assign GC content of constructed read to last read ID (that read was associated with)
                    read_dict[read_id] = read
                # if read_id is None, this is the first line, so just assign the read ID, 
                # we'll get GC content in next pass
                read_id = line.lstrip('>').strip()
                break
            else:
                # construct the read by iterating through until we hit None (EOF) or next read ID ('>blahblah')
                read += line.strip()
                line = next(lines, None)
    return read_dict


def random_genetic_sequence(length=150, alphabet=DNA_ALPHABET):
    return ''.join([choice(alphabet) for _ in range(length)])


def get_reverse_complement(seq):
    return ''.join([RC_DICT[x] for x in seq[::-1]])
