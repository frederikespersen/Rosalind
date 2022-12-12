def alph_strings_lex(alphabet: tuple, n: int) -> list:
    assert n <= len(alphabet), "n may not be greater than the length of the alphabet provided!"

    for triplet in [alphabet[i:i+3] for i in range(len(alphabet)-2)]:
        assert triplet[0] < triplet[1] < triplet[2], "Alphabet must follow the order of the English alphabet!"

    strings = [*alphabet]
    l = 1
    while l < n:
        _strings = strings
        strings = []
        for char in alphabet:
            strings += [string + char for string in _strings]
        l += 1

    return sorted(strings)


n = 4
input = 'A B C D'
if __name__ == '__main__':
    alphabet = tuple(input.split())
    strings = alph_strings_lex(alphabet, n)
    print(*strings, sep='\n')