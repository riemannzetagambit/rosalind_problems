import re
import sys

import requests

from rosalind_utils import get_rosalind_data, process_fasta_file

_UNIPROT_URL = 'https://www.uniprot.org/uniprot/{uid}.fasta'
# this is the pattern for the N-glycosylation motif; the look-ahead (?=) catches overlapping matches
_N_GLYCOSYLATION_REGEX = re.compile('(?=N[^P][ST][^P])')


def solve_problem(protein_data):
    '''
    Assumptions: input data is of the form
    uniprot_id_0
    uniprot_id_1
    ...
    '''
    n_glycosylation_matches = {}
    for uniprot_id in protein_data:
        resp = requests.get(_UNIPROT_URL.format(uid=uniprot_id))
        # clean up their response files into lines expected for fasta; python 3 returns these as bytes, so we decode
        fasta_lines = resp.content.decode('utf-8').strip().split('\n')
        seq_dict = process_fasta_file(fasta_lines)
        protein_seq = list(seq_dict.values())[0]
        motifs = _N_GLYCOSYLATION_REGEX.finditer(protein_seq)
        motif_loci = []
        for motif in motifs:
            # add 1 since python indexes from 0
            motif_loci.append(motif.start() + 1)
        if len(motif_loci) > 0:
            n_glycosylation_matches[uniprot_id] = motif_loci

    print('\n'.join('{uid}\n{loci}'.format(uid=uid, loci=' '.join(map(str, loci))) 
                    for uid, loci in n_glycosylation_matches.items()))
    return '\n'.join('{uid}\n{loci}'.format(uid=uid, loci=' '.join(map(str, loci))) 
                     for uid, loci in n_glycosylation_matches.items())


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
