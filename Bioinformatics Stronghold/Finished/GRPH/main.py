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


def overlap_graph_adjlst(nodes: dict, k: int) -> list:
    adjlts = []
    for s in nodes:
        suffix = nodes[s][-k:]
        for p in nodes:
            if p == s:
                continue
            prefix = nodes[p][:k]
            if suffix == prefix:
                adjlts.append([s, p])
    return adjlts


fasta_filename = 'strings.fasta'
k = 3
if __name__ == '__main__':
    parsed = parse_fasta(fasta_filename)
    adjlst = overlap_graph_adjlst(parsed, k)
    for dir_edg in adjlst:
        print(*dir_edg)