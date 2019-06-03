from typing import List
import sys

from rosalind_utils import convert_dna_to_protein, get_rosalind_data, process_fasta_file

def _splice_read(read: str, introns: List[str]) -> str:
    '''
    Given a read (gene) and its introns, return the exons of the gene by removing the introns

    :param read: str
        the read to splice
    :param introns: list of strings
        list of strings of the introns to remove from the read

    :return: str
        the spliced read
    '''
    for intron in introns:
        read = read.replace(intron, '')

    return read

def solve_problem(sequence_data: list) -> str:
    '''
    Assumptions: input file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    ...
    '''
    reads_dict = process_fasta_file(sequence_data)
    reads = list(reads_dict.values())
    # since process_fasta_file returns an OrderedDict, we are assured that the first value is the read
    read, introns = reads[0], reads[1:]
    spliced_read = _splice_read(read, introns)
    # because we use the generator twice below, we need to cast this to a list lest we exhaust it
    proteins = list(convert_dna_to_protein(spliced_read))

    print('\n'.join(proteins))
    return '\n'.join(proteins)


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
