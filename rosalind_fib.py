import sys

from .rosalind_utils import get_rosalind_data

# memoize
fib_dict = {1: 1, 2: 1}

def solve_problem(n, k):
    if n <= 2:
        return fib_dict[n]
    else:
        # here is the critical factor that is different from usual Fibonacci-- k in front of second sum
        fib_dict[n] = solve_problem(n-1, k) + k*solve_problem(n-2, k)

    return fib_dict[n]


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    n, k = map(int, get_rosalind_data(args[-1])[0].split())
    final_count = solve_problem(n, k)
    print(final_count)
