
# The minimum unrooted tree has three leaves l and one internal node i:
# i_3 = 1
# To get one more leaf, the tree must turn one node (leaf) into an internal node; splitting it:
# i_l = i_l-1 + 1
# Thus:
# i_l-i_l-1 = 1
# i_l-i_l-d = d  or  i_l+d-i_l = d
# i_l+d = i_l + d
# Set l = 3:
# i_3+d = 1 + d
# Set n = d + 3 -> d = n - 3
# i_n = 1 + n - 3
# i_n = n - 2
def internal_nodes(n: int) -> int:
    assert n > 2
    return n - 2


n = 8513
if __name__ == '__main__':
    print(internal_nodes(n))