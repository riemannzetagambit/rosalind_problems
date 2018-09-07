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
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]

    '''
    gc_content_dict = {}
    read_id = None
    lines = iter(sequence_data)
    for line in lines:
        read = ''
        # iterate through lines to find identifier or set of lines corresponding to reads
        while True:
            if line is None:
                # at end of file, so assign read GC content to last known read ID
                gc_content_dict[read_id] = _get_gc_content(read)
                break
            elif line.startswith('>'):
                if read_id is not None:
                    # we already have a read ID from last read ID assignment and have been constructing read
                    # so just assign GC content of constructed read to last read ID (that read was associated with)
                    gc_content_dict[read_id] = _get_gc_content(read)
                # if read_id is None, this is the first line, so just assign the read ID, 
                # we'll get GC content in next pass
                read_id = line.lstrip('>').strip()
                break
            else:
                # construct the read by iterating through until we hit None (EOF) or next read ID ('>blahblah')
                read += line.strip()
                line = next(lines, None)

    max_read_id = max(gc_content_dict, key=gc_content_dict.get)
    print(f'{max_read_id}\n{gc_content_dict[max_read_id]}')
    
    return max_read_id, gc_content_dict[max_read_id]


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
