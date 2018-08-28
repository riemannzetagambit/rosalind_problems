from collections import Counter
import sys

from rosalind_utils import get_rosalind_data, RC_DICT


def _get_gc_content(dna_seq):
    '''
    Compute GC content using Counter from collections
    Faster than just string parsing in _alt_get_gc_content
    '''
    nt_counts = dict(Counter(dna_seq).most_common())
    return 100*(nt_counts.get('G', 0) + nt_counts.get('C', 0))/sum(nt_counts.values())


def _alt_get_gc_content(dna_seq):
    '''
    Get GC content without any imports, slower than Counter
    '''
    return 100*sum(1 for x in dna_seq if x == 'G' or x == 'C')/len(dna_seq)


def solve_problem(sequence_data):
    read_id_dict = {}
    lines = iter(sequence_data)
    for line in lines:
        if line.startswith('>'):
            # read in pairs, read name, next line is read content
            read_id = line.lstrip('>').strip()
            read_id_dict[read_id] = next(lines).strip() # user iter above to use next() capability

    print(read_id_dict)
    print({k: _get_gc_content(v) for k, v in read_id_dict.items()})
    print({k: _alt_get_gc_content(v) for k, v in read_id_dict.items()})
    # could also compute {k: _get_gc_content(v) for k, v in read_id_dict.items()}, grab max key, and then grab value
    # therein to save an extra _get_gc_content call below
    max_read_id = max(read_id_dict, key=lambda x: _get_gc_content(read_id_dict[x]))
    gc_content = _get_gc_content(read_id_dict[max_read_id])
    print(f'{max_read_id}\n{gc_content}')
    
    return max_read_id, gc_content


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
