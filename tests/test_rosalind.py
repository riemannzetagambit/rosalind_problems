import os

from ..rosalind_utils import get_rosalind_data
from ..rosalind_corr import solve_problem as corr_solve
from ..rosalind_long import solve_problem as long_solve


def test_rosalind_long():
    solution = 'ATTAGACCTGCCGGAATAC'
    test_data = get_rosalind_data(os.path.abspath('./rosalind_long_test.txt'))
    assert long_solve(test_data) == solution

def test_rosalind_corr():
    solution = 'TTCAT->TTGAT\nGAGGA->GATGA\nTTTCC->TTTCA'
    test_data = get_rosalind_data(os.path.abspath('./rosalind_corr_test.txt'))
    assert corr_solve(test_data) == solution
