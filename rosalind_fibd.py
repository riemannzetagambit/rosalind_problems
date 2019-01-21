import sys

from rosalind_utils import get_rosalind_data

# memoize
fib_dict = {1: 1, 2: 1}

def count_rabbits(n, m):
    '''

    '''
    if n <= 2:
        return fib_dict[n]
    else:
        # here is the critical factor that is different from usual Fibonacci-- k in front of second sum
        fib_dict[n] = solve_problem(n-1, k) + k*solve_problem(n-2, k)

    return fib_dict[n]


def solve_problem(fib_data):
    n, m = map(int, fib_data.split())
    final_count = solve_problem(n, m)

    print(final_count)
    return final_count


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
