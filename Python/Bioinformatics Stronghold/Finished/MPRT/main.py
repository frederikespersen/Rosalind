def get_fasta(uniprot_id: str) -> str:
    import urllib.request as url
    uniprot_link = 'http://www.uniprot.org/uniprot/'
    try:
        fasta = url.urlopen(uniprot_link + uniprot_id + '.fasta')
    except:
        return ''

    pseq = ''.join([str(l) for l in fasta.readlines()[1:]])
    for s in ["'", "b", "\\n"]:
        pseq = pseq.replace(s, "")
    return pseq


def format_motif(motif: str) -> list:
    all_aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    motif_elements = []

    oneof = False
    noneof = False
    e = ''
    for i in motif:

        if i == '{':
            noneof = not noneof

        if i in ['[',']']:
            oneof = not oneof
        else:
            e += i

        if i == '}':
            noneof = not noneof
            e_ = e
            e = ''
            for aa in all_aa:
                if aa not in e_:
                    e += aa

        if not oneof and not noneof:
            motif_elements.append(e)
            e = ''

    motif_variants = ['']
    for e in motif_elements:
        variants = motif_variants
        motif_variants = []
        for aa in e:
            motif_variants += [var + aa for var in variants]

    return motif_variants


def locate_motif(seq, motif):
    motif_variants = format_motif(motif)
    ml = len(motif_variants[0])
    locations = []
    for p in range(len(seq) - ml + 1):
        subseq = seq[p:p + ml]
        if subseq in motif_variants:
            locations.append(p + 1)
    return locations


uniprot_ids_filename = 'uniprot_ids.txt'
motif = 'N{P}[ST]{P}'

if __name__ == '__main__':
    with open(uniprot_ids_filename, 'r') as file:
        ids = [line.strip() for line in file.readlines()]

    for id in ids:
        seq = get_fasta(id.split('_')[0])
        motif_matches = locate_motif(seq, motif)
        if len(motif_matches) > 0:
            print(id)
            print(*motif_matches)
