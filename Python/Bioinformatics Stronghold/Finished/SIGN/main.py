def permutations(n: int) -> list:
    nums = [str(i+1) for i in range(n)]
    perms = nums + ['-' + num for num in nums]
    for i in range(n-1):
        _perms = perms
        perms = []
        for num in nums:
            for perm in _perms:
                if num not in perm:
                    perms += [perm + ' ' + num]
                    perms += [perm + ' -' + num]
    return sorted(perms)


n = 6
if __name__ == '__main__':
    perms = permutations(n)
    spaced_perms = [' '.join(list(perm)) for perm in perms]
    with open('perms.txt', 'w') as file:
        file.writelines(f'{len(spaced_perms)}\n')
        file.writelines([perm + '\n' for perm in perms])
