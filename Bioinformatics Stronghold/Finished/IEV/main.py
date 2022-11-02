def perm_dom_prom(a: str, b: str):
    combination = [a, b]
    assert len(combination) == 2
    for a in combination:
        assert a in ['AA', 'Aa', 'aa']
    if 'AA' in combination:  # AA-AA, AA-Aa, AA-aa
        return 1
    elif 'Aa' in combination:
        if 'aa' in combination:  # Aa-aa
            return 0.5
        else:  # Aa-Aa
            return 0.75
    else:  # aa-aa
        return 0


def dom_phenotype_prob(genotypes_input):
    a, b, c, d, e, f = (int(i) for i in genotypes_input.split(' '))
    genotype = {'AA-AA': a,
                'AA-Aa': b,
                'AA-aa': c,
                'Aa-Aa': d,
                'Aa-aa': e,
                'aa-aa': f}

    pair_no_offspring = 2
    expected_dom_offspring = 0
    for pair in genotype:
        no_pairs = genotype[pair]
        p_dom = perm_dom_prom(*pair.split('-'))
        expected_dom_offspring += no_pairs * pair_no_offspring * p_dom

    return expected_dom_offspring


genotypes = "17851 16753 17552 18512 17419 17172"

if __name__ == '__main__':
    print(dom_phenotype_prob(genotypes))
