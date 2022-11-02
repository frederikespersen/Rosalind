def alph_strings_lex(alphabet: list, n: int) -> list:
    assert n <= len(alphabet), "n may not be greater than the length of the alphabet provided!"

    # Assembling strings
    strings = [*alphabet]
    for l in range(n - 1):
        _strings = []
        for char in alphabet:
            _strings += [string + char for string in strings]
        strings += _strings
    strings = list(set(strings))

    # Sorting strings
    filler = 'Ø'
    alphabet = [filler] + alphabet
    order = dict(zip(alphabet, range(len(alphabet))))
    strings = [string + filler * (n - len(string)) for string in strings]
    strings.sort(key=lambda string: [order[char] for char in string])
    strings = [string.replace(filler, '') for string in strings]

    return strings


n = 4
input = 'W F X V U B O Y Z Q E'
if __name__ == '__main__':
    alphabet = input.split()
    strings = alph_strings_lex(alphabet, n)
    print(*strings, sep='\n')