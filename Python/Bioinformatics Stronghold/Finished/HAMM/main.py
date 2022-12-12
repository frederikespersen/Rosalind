def single_substitutions(seq1: str, seq2: str) -> int:
    assert len(seq1) == len(seq2)
    count = 0
    for b in range(len(seq1)):
        if seq1[b] != seq2[b]:
            count += 1
    return count


strings_filename = 'strings.txt'

if __name__ == '__main__':
    with open(strings_filename, 'r') as file:
        strings = [line.strip() for line in file.readlines()]

    print(single_substitutions(strings[0], strings[1]))
