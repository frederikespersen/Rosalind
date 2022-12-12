def final_pop(n: int, k: int) -> int:
    pop = []
    for m in range(n):
        if m < 2:
            pop.append(1)
        else:
            pop.append(pop[m - 1] + k * pop[m - 2])
    return pop[-1]


sample_input = '29 2'
if __name__ == '__main__':
    n, k = (int(i) for i in sample_input.split(' '))
    print(final_pop(n, k))
