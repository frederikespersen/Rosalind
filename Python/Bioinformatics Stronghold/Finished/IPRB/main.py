def perm_dom_prom(a: str, b: str):
    combination = [a, b]
    assert len(combination) == 2
    for a in combination:
        assert a in ['k', 'm', 'n']
    if 'k' in combination:  # kk, km, kn
        return 1
    elif 'm' in combination:
        if 'n' in combination:  # mn
            return 0.5
        else:  # mm
            return 0.75
    else:  # nn
        return 0


def dom_phenotype_prob(genotypes_input: str):
    k, m, n = (int(i) for i in genotypes_input.split(' '))
    pop = k + m + n
    genotypes = {'k': k,
                 'm': m,
                 'n': n}

    dom_prob = 0
    for a in genotypes:
        for b in genotypes:
            p_a = genotypes[a] / pop
            p_b = genotypes[b] / (pop - 1)
            if a == b:
                p_b -= 1 / (pop - 1)
            p_ab = p_a * p_b

            p_ab_dom = perm_dom_prom(a, b)

            dom_prob += p_ab * p_ab_dom

    return dom_prob


if __name__ == '__main__':
    print(dom_phenotype_prob("28 26 25"))
