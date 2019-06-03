from pathlib import Path
import sys

# remove the need for relative .s and add directory above tests to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rosalind_utils import get_rosalind_data
from rosalind_corr import solve_problem as corr_solve
from rosalind_long import solve_problem as long_solve
from rosalind_mprt import solve_problem as mprt_solve
from rosalind_orf import solve_problem as orf_solve
from rosalind_prob import solve_problem as prob_solve
from rosalind_splc import solve_problem as splc_solve


def is_close(a, b, rel_tol=1e-4):
    '''
    Function to assert two numbers are 'close enough' to be equal when equality is not guaranteed
    (without using numpy implementation, though there is a math.isclose function)
    '''
    return abs((a - b)/a) <= rel_tol

def test_rosalind_corr():
    solution = 'TTCAT->TTGAT\nGAGGA->GATGA\nTTTCC->TTTCA'
    test_path = Path('{}/../rosalind_corr_test.txt'.format(__file__)).resolve()
    test_data = get_rosalind_data(test_path)
    # account for rearranging
    assert sorted(corr_solve(test_data).split()) == sorted(solution.split())

def test_rosalind_long():
    solution = 'ATTAGACCTGCCGGAATAC'
    test_path = Path('{}/../rosalind_long_test.txt'.format(__file__)).resolve()
    test_data = get_rosalind_data(test_path)
    assert long_solve(test_data) == solution

def test_rosalind_mprt():
    solution = '''B5ZC00\n85 118 142 306 395\nP07204_TRBM_HUMAN\n47 115 116 382 409
                  \nP20840_SAG1_YEAST\n79 109 135 248 306 348 364 402 485 501 614'''
    test_path = Path('{}/../rosalind_mprt_test.txt'.format(__file__)).resolve()
    test_data = get_rosalind_data(test_path)
    assert sorted(mprt_solve(test_data).split()) == sorted(solution.split())

def test_rosalind_prob():
    solution = '-5.737 -5.217 -5.263 -5.360 -5.958 -6.628 -7.009'
    test_path = Path('{}/../rosalind_prob_test.txt'.format(__file__)).resolve()
    test_data = get_rosalind_data(test_path)
    # test every individual value
    prob_solution = prob_solve(test_data)
    for prob_sol, test_sol in zip(prob_solution.split(), solution.split()):
        assert is_close(float(prob_sol), float(test_sol))

def test_rosalind_orf():
    solution = 'MLLGSFRLIPKETLIQVAGSSPCNLS\nM\nMGMTPRLGLESLLE\nMTPRLGLESLLE'
    test_path = Path('{}/../rosalind_orf_test.txt'.format(__file__)).resolve()
    test_data = get_rosalind_data(test_path)
    # account for rearranging
    assert sorted(orf_solve(test_data).split()) == sorted(solution.split())

def test_rosalind_splc():
    solution = 'MVYIADKQHVASREAYGHMFKVCA'
    test_path = Path('{}/../rosalind_splc_test.txt'.format(__file__)).resolve()
    test_data = get_rosalind_data(test_path)
    assert splc_solve(test_data) == solution
