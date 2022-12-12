def final_pop(n: int, k: int, m: int) -> int:
    pop = []
    a, b, d = 1, 0, 0  # Adult, Born, Died
    a_, b_, d_ = [], [], []

    for i in range(n):

        if i > 0:
            a = pop[i - 1]

        if i > 1:
            b = k * (a_[i - 1] - d_[i - 1])

        if i == m:
            d = a_[0]
        elif i > m:
            d = b_[i - m]

        pop.append(a + b - d)
        a_.append(a)
        b_.append(b)
        d_.append(d)

    return pop[-1]


sample_input = '80 19'
if __name__ == '__main__':
    k = 1
    n, m = (int(i) for i in sample_input.split(' '))
    print(final_pop(n, k, m))
