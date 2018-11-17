import os
import sys

# remove the need for relative .s and add directory above tests to path
sys.path.insert(0, '..')

from rosalind_utils import get_rosalind_data
from rosalind_corr import solve_problem as corr_solve
from rosalind_long import solve_problem as long_solve
#from .rosalind_prob import solve_problem as prob_solve

file_path = os.path.dirname(os.path.abspath(__file__))

def test_rosalind_corr():
    solution = 'TTCAT->TTGAT\nGAGGA->GATGA\nTTTCC->TTTCA'
    test_path = os.path.join(file_path, './rosalind_corr_test.txt')
    test_data = get_rosalind_data(test_path)
    assert corr_solve(test_data) == solution

def test_rosalind_long():
    solution = 'ATTAGACCTGCCGGAATAC'
    test_path = os.path.join(file_path, './rosalind_long_test.txt')
    test_data = get_rosalind_data(test_path)
    assert long_solve(test_data) == solution

#def test_rosalind_prob():
#    solution = '-5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009'
#    test_data = get_rosalind_data(os.path.abspath('./rosalind_prob_test.txt'))
#    assert prob_solve(test_data) == solution
