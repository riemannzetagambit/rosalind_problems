from ..rosalind_utils import get_rosalind_data
from ..rosalind_long import solve_problem as long_solve


def test_rosalind_long():
    solution = 'ATTAGACCTGCCGGAATAC'
    test_data = get_rosalind_data('./rosalind_long_test.txt')
    assert long_solve(test_data) == solution
