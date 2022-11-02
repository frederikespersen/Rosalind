def mrna_permutations(seq: str) -> int:
    permutations = 1
    for aa in seq:
        permutations *= aas[aa]
    permutations *= aas['Stop']
    return permutations


codon_file = 'codons.txt'
sample_protein = 'MWRSYVYSCHHGEPNGDYCVWTLGLSRILVRISKMPARQAATMNIINYMCSFDNPYAMSPYYTAHTCKSQVLPYFRQLCLATSEVGTNVKERHIPYVSSCCFQQTYFSGWFCYISRPWLIWVHVFHQHGMEAHLQCNEPLKISRPKDPNETMLQEFMIARSENRYQTQHPMEGQSTSNCQDCEFPCGVEQECWQVDAAVHPTKWDFKFSGEKWFFRVQRHGGWQADCFFECMMSKPINMLWCFYVWYVLYWEMDGHQDKQGVHYYICENIPHFHEHSYRMKLYTDRGMDWDTNFRGQQYYNMSGEQPNHDVLCHGWHNDQSLCMNQMDKMANSRHWDAMSEMPGFYHSNTKKLCHSDVVGNAILFVGLERATEANLIDDNYWNIIMIAVDYGRLYQVKRIPNLMNSNHCTCIDGFNHEISFICNVQAWVIGMRRMFKGQSCFHEPVLIGTMSNMGDQQMKQYYICQAELVRKSLHDYFGNWWCGHIAMYDMSRQSLEPGFDEKYNDQWFMFYNLLIDYARWRELMPQHQLKYHVATQWEFRWGEADAYCLAQYIWNYPGDDWINGASEHEAYKCNGEACECFVVIKGAIIGKGWRPQDEWDKKVHCKSLFRGTEWQRQGQSMSGMIEHDIVDRKCFDYRNNSKTHGIFAHTHDMLDKATFYFKICSFDPRKQMQYECSRFGNVLCTAEPVDVLSRVSMMQNYWADNVIMTVFYPTECEMSGCTGQHWTPMENWGVYKIENTHSCVPHRLENGMCEEMIMVMYQAFWQHMSELKLHLGREGYQGEEVQYIPVAWFTFHPNNCEVCHGLCDPQMWPLKNRYTMRKQIECPLAYEYNHGTDSWYKSIEVMYAYSPKEWATKIHAEYRYWNKIKARYFYTLMGMVPQVDTHKDQTEGIMQMNPQLFTRMLTRKRTIHTRRANQIERYFRAHHCMDLEREQYIASQSIGMMTIHAVWLCTMVVSSFLNIHDNCTPTIYIFAFLLSAITTMASSTFSSLCTMGVDG'

if __name__ == '__main__':

    # Importing codons
    with open(codon_file, 'r') as file:
        lines = ' '.join([line.strip() for line in file.readlines()])
    while '  ' in lines:
        lines = lines.replace('  ', ' ')
    codons = {}
    for i in lines.split(' '):
        if len(i) == 3:
            codon = i
        else:
            codons[codon] = i

    # Mapping codons to amino acids
    aas = {}
    for c in codons:
        aa = codons[c]
        if aa not in aas:
            aas[aa] = []
        aas[aa].append(c)

    # Aggregating codon possibilities to count:
    aas = {aa: len(cs) for aa, cs in aas.items()}

    # Calculating permutations
    mod = 1000000
    print(mrna_permutations(sample_protein) % mod)
