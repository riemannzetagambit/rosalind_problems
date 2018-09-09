import sys

from rosalind_utils import get_rosalind_data

def _get_dominant_allele_probability(k, m, n):
    '''
    Compute the probability that in a population of k homozygous dominant, m heterozygous, and n homozygous recessive
    (i.e. k + m + n total population), two random selections will produce offspring with at least one dominant allele

    Assume k, m, n are all positive semidefinite
    '''
    total_pop = k + m + n
    prob_two_dominant = (k/total_pop)*((k-1)/(total_pop-1))
    # can choose heterozygous individual first or second, then homozygous, so just multiply by factor of two 
    # (and sim. later)
    prob_one_dominant_one_hetero = 2*(k/total_pop)*(m/(total_pop-1))
    prob_one_dominant_one_recessive = 2*(k/total_pop)*(n/(total_pop-1))
    prob_two_hetero = (m/total_pop)*((m-1)/(total_pop-1))
    prob_one_hetero_one_recessive = 2*(m/total_pop)*(n/(total_pop-1))
    prob_two_recessive = (n/total_pop)*((n-1)/(total_pop-1))

    # draw Mendelian squares for this logic
    # 2 homoz. dominant individuals always produce dominant offspring, so multiply by factor 1 for this combo
    # 1 homoz. dominant 1 heteroz. produce dominant offspring every time because of GG in dom. (GG x Gg)
    # 1 homoz. dominant 1 homoz. recessive produce dominant offspring every time because of GG in dom. (GG x gg)
    # two heterozygous individuals have a dominant allele in 3/4 of combinations (Gg x Gg has a G in 3 of the combos)
    # 1 heteroz. 1 homoz. recessive have dominant allele in 2/3 combinations (Gg x gg has G in 2 of combos)
    # 2 homoz. recessive have 0 combos with dominant allele
    return (1*prob_two_dominant + 1*prob_one_dominant_one_hetero + 1*prob_one_dominant_one_recessive
            + (3/4)*prob_two_hetero + (2/4)*prob_one_hetero_one_recessive
            + (0/4)*prob_two_recessive)


def solve_problem(sequence_data):
    '''
    Assumes only one line in file
    '''
    k, m, n = map(int, sequence_data[0].split())
    prob_dominant_allele = _get_dominant_allele_probability(k, m, n)

    return prob_dominant_allele


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) != 2:
        raise ValueError('Must provide a single argument which is the data file path')
    solve_problem(get_rosalind_data(args[-1]))
