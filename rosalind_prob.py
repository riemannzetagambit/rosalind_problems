from functools import reduce
from math import log
import operator
import sys

from rosalind_utils import get_rosalind_data


# this is not a python builtin, so make our own
# see https://stackoverflow.com/questions/595374/whats-the-function-like-sum-but-for-multiplication-product
def prod(iterable):
    return reduce(operator.mul, iterable, 1)

def solve_problem(sequence_data):
    '''
    Assumptions: file is of the form
    >read_id_0
    [arbitrary number of lines of ATCG]
    >read_id_1
    [arbitrary number of lines of ATCG]

    '''
    read, probs = sequence_data[0], sequence_data[1]
    # split space-delimited probabilities into list and cast them to floats with an iterable
    probs = map(float, probs.split())

    log_prob_random_str = []
    for prob in probs:
        prob_dict = {}
        # the given probability is the GC content, so prob(G) = prob(c) = prob/2
        prob_dict['G'], prob_dict['C'] = prob/2, prob/2
        # by conservation of probability
        prob_dict['A'], prob_dict['T'] = (1 - prob)/2, (1 - prob)/2
        # multiply together the probability of each individual nt
        prob_read = prod([prob_dict[x] for x in read])
        log_prob_random_str.append(log(prob_read, 10))

    print(' '.join(map(str, log_prob_random_str)))
    return ' '.join(map(str, log_prob_random_str))


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
