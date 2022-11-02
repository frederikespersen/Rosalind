def parse_fasta(filename: str) -> dict:
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    strings = {}
    for line in lines:
        if line[0] == '>':
            k = line[1:]
            strings[k] = ''
        else:
            strings[k] += line
    return strings


def dna_complimentary(sequence: str):
    complimentary = ''
    for p in range(len(sequence)):
        b = sequence[p]
        if b == 'A':
            complimentary += 'T'
        elif b == 'T':
            complimentary += 'A'
        elif b == 'G':
            complimentary += 'C'
        elif b == 'C':
            complimentary += 'G'
    complimentary = complimentary[::-1]
    assert len(sequence) == len(complimentary)
    return complimentary


def load_triplet_codes(filename: str) -> dict:
    with open(filename, 'r') as file:
        lines = ' '.join([line.strip() for line in file.readlines()])
    while '  ' in lines:
        lines = lines.replace('  ', ' ')
    triplets = {}
    for i in lines.split(' '):
        if len(i) == 3:
            triplet = i
        else:
            triplets[triplet] = i
    return triplets


def orfs(strand: str, triplet_code_filename: str) -> list:
    template = strand
    coding = dna_complimentary(strand)
    triplet_codes = load_triplet_codes(triplet_code_filename)
    orfs = []
    for string in [template, coding]:
        for f in [0, 1, 2]:
            frame = string[f:]
            codons = [frame[i:i + 3] for i in range(0, len(frame), 3) if i + 3 <= len(frame)]
            frame_aas = [triplet_codes[triplet.replace('T', 'U')] for triplet in codons]  # Template/Coding strand switch
            translating = False
            aas = ''
            for aa in frame_aas:
                if aa == 'M':
                    translating = True
                if aa == 'Stop':
                    translating = False
                    if aas != '':
                        for i in range(len(aas)):  # In case of orfs in orfs
                            if aas[i] == 'M':
                                orfs.append(aas[i:])
                    aas = ''
                    continue
                if translating:
                    aas += aa
    unique_orfs = list(set(orfs))
    return unique_orfs


fasta_filename = 'string.fasta'
triplet_code_filename = 'codons.txt'
if __name__ == '__main__':
    string = [*parse_fasta(fasta_filename).values()][0]
    orfs = orfs(string, triplet_code_filename)
    print(*orfs, sep='\n')
