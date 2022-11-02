def load_dna_strings(filename: str) -> dict:
    with open(filename, 'r') as file:
        lines = file.readlines()

    strings = {}
    for line in lines:
        line = line.rstrip('\n')
        if line[0] == '>':
            s = line[1:]
            strings[s] = ''
        else:
            strings[s] += line

    return strings


def gc_content(seq: str) -> float:
    gc = 0
    for p in range(len(seq)):
        if seq[p] in ['G', 'C']:
            gc += 1
    return 100 * gc / len(seq)


if __name__ == '__main__':
    dna_strings = load_dna_strings('sample.fasta')
    gc_strings = {}
    for k in dna_strings:
        gc_strings[k] = gc_content(dna_strings[k])
    k_max = max(gc_strings, key=gc_strings.get)
    print(k_max, gc_strings[k_max], sep='\n')
