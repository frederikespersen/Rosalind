def prob_at_least_N(k: int, N: int):
    # All AaBb-XXXX matings always has a 25% probability of AaBb offspring (Punnett square).
    # This is a binomial experiment;
    # We are looking for the probability of at least N 'successes' (AaBb) with 2^k trials and success rate of 25%.
    # p(N ≥ AaBb)
    # This is equal to 1 minus the probability of less than N AaBb with 2^k and success rate of 25%
    # p(N ≥ AaBb) = 1 - p(N < AaBb)
    # p(N ≥ AaBb) = 1 - p(N - 1 ≤ AaBb)   [Note that N is a natural number]
    # The right hand side can be nicely expressed by the binomial cumulative distribution function (bcdf)
    # p(N ≥ AaBb) = 1 - bcdf(N - 1, 2^k, p)

    p_AaBb = 0.25

    from scipy.stats import binom
    return 1 - binom.cdf(N - 1, 2**k, p_AaBb)


sample_input = '5 9'
if __name__=='__main__':
    k, N = (int(i) for i in sample_input.split(' '))
    print(prob_at_least_N(k, N))

