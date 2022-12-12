# For every subset, each element of the set is either included or not.
# This makes the total amount of subsets binarily exponetial
def no_subsets(no_elements: int) -> int:
    return 2**no_elements % 1000000


n = 866
if __name__ == '__main__':
    print(no_subsets(n))