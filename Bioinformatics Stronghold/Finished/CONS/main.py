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


def consensus_string(strings: list):
    import numpy as np
    nb = np.array(['A', 'C', 'G', 'T'])
    n = len(strings[0])
    profile_matrix = np.zeros((4, n))
    for string in strings:
        for j in range(n):
            for i in range(4):
                if nb[i] == string[j]:
                    profile_matrix[i, j] += 1
    cons_strings = ['']
    for j in range(n):
        m = max(profile_matrix[:, j])
        cons_strings_ = cons_strings
        cons_strings = []
        m_profile = profile_matrix[:, j] == m
        m_indices = np.transpose(m_profile.nonzero())
        #for i in [i[0] for i in m_indices]:  # Consensus string permutations increases exponentially for multiple hits
        for i in [i[0] for i in [m_indices[0]]]:
            cons_strings += [string + nb[i] for string in cons_strings_]
    return cons_strings, profile_matrix


fasta_filename = 'strings.fasta'
if __name__ == '__main__':
    parsed = parse_fasta(fasta_filename)
    cons_strings, profile_matrix = consensus_string([*parsed.values()])
    print(cons_strings[0])
    print('A:', *[int(i) for i in profile_matrix[0]])
    print('C:', *[int(i) for i in profile_matrix[1]])
    print('G:', *[int(i) for i in profile_matrix[2]])
    print('T:', *[int(i) for i in profile_matrix[3]])