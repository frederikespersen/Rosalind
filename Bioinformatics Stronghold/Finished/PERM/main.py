def permutations(n: int) -> list:
    nums = [str(i+1) for i in range(n)]
    perms = nums
    for i in range(n-1):
        _perms = perms
        perms = []
        for num in nums:
            perms += [perm + num for perm in _perms if num not in perm]
    return perms


n = 7
if __name__ == '__main__':
    perms = permutations(n)
    spaces_perms = [' '.join(list(perm)) for perm in perms]
    with open('perms.txt', 'w') as file:
        file.write(f'{len(perms)}\n')  # Writing number of permutations
        file.writelines([perm + '\n' for perm in spaces_perms])  # Listing permutations