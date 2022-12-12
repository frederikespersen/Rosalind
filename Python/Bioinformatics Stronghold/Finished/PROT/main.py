def translation(seq: str) -> str:
    aas = ''
    for c in range(len(seq) // 3):
        triplet = seq[3*c:3*(c+1)]
        aa = codons[triplet]
        if aa == 'Stop':
            break
        aas += aa
    return aas


codon_file = 'codons.txt'
sample_rna = 'AUGUCACUCCAUGUUGGUCCUCCAUUAAUAUGUGCUAUAUUGACAGUUGGAAACGGAUUAAGUACUUCUUCGAUUGUGGAGGCCAUAAACGGCCAAGAACCUGCGCGCUUACAUACUUCACAAAACACAUCAGCAGGUGAAUCUCAGUUGACACCUUCCGCCCCGUCCCAUGGGGGUCGCUCUCUGGCAGUUUCUCCGCUAGUCAGAUGUGUAUCGUUUCAGGUGCGGCUCUACAUAGUUCCCACCAGCUGCCUAGUACAAUAUCCGAAAGAAUCUGGCCAUAGGGCACGCACUAUCAUAAUAAAAUUCAGAACACGAAUGCGAGCCCAUUCAGACCCGGCUUCAAAGAGGAAAGAAAUAUGUAGGGGGGCGGUAUUCCUCAAGUCUGUUACGGUGUCCAAAACGCCGGAGUUGAAUUACGCAAUGGUGGUUCAGCGUACGCGAGAUCACUCCAUCAUCGAAGUUUCUCCCUUAAACGUGGCCUCCCGCGAUUCCCAAGCUUGCUGCCAUGGUCCUCCGCUUCCUCACGCUGACAACCUGUCAAACUCACUGUUCCAUGGUCUUUCGAUAGAAUAUAGACUCCCCUCAAAUUCGGAGGAUUCAUCACGUUCUUCGCAAACGAAUGUAGGUGGAGAAGAGACCAAGAGUCCAGUCUUGAGGGAUUUUGUCGCGCCUCAUGGGAUUCCAUCGUUUUCGGUAAGUUCCUCAGGUGACGCCGGGCCCAUUACGAGAACUCAACGCUACUCUCAACCUUCGCCACACCAGGAUGCUUCCCCGGUUGAAGACGAACUACCACGCGAUGUCAAACUAGACAGUAGCCUCAAUUUUCGGGGUGUGCGUCCUAAGCGGAUCGAUGUCGCUAGAAUUAUGGGGAAUUUGGACAAGCGACGUAAAUGGUGCAUCCACCCCUCGAAACUAACUGAUGGCGCCCAUUCAGCAAUUCCCCCUUCCCCAAGUUUACGAGCCCCCAUGACUUGCGAUACAAAGUUCAGAUCUUUUAUACCUGCUAUCUGUCUAGGCACGGAGGCGGGUUCGCACACUUCGACGCUUGCCGACCGCUACCAAAAAAUUGGAGGAAUAUCUUGCCCGAAAGUACAUCACUUACGGACAGUGGUUCUCACAAGGUGCGCCGAUUAUUUGCGCCCUCUCGAUCAGCUUCCGCGAUCCCGCGCUUUAGACCAAAAACUCAAGAUUUGUAUUACCCGCGAUUCGUUGGCGAUUUUAGACGGAUUGUACUCUGGGACGAUAAGGUACAAAAGCUUACACACCGCAGGGCCUGAUUCAGCUUCAGAGAACGCUAGCUAUCAGGCCGUGAUUUGGUCUAACUAUCAAUGGCUCAGUGUAACCGAAAUUACAGUCAAAAACCGCUGGACCCAUAAUCGCUAUGUUCCUCGCUUUUUCUAUUUGGCCACUACCAUAAUACUUGCCUAUGCUUCACUACAAGGGGACGUUCUCGGUCUCCUAGCCCCGCCUUCCGCUCUAAGGACGACUCAAGAUACGUUGCUGAGACGACUGGUUUUAAAGGCCUCACUCACGGCUAAAAACGGAACUCUAAGCUUGAUAUCAGCCACGUGUCAGAAUUUACGCGGCGGGCAUUUCCGCAACGGCGGUCCUGACGUCCACGGACGGUCAGUUAUACACAGCGUACCUAGUGCGUCCUGUGGCACCCUCGUGCAUAAGCCGCGACAUAAGUACCUGCCCUCAUGUAGGCGAGUCCGUAAGAUGGCUGCAGGAUACCUCACCGAAAUGGCCUUACGCAGUAAUUCCCUCUUGGUGAUCGUCACGCUCCGUGGAAGGCACGAUUACAUGCAUCGAGGCAUGGAACUUAUCAAAAGAACAGUGCCGAAUAUGUGGAAAGUGAAAGGUAGGCCUGUUGUCAUCCUGUUGCCUUGGCUGGCCGUCAAGGGGCCACAGUCUAGCAGGCCUUGUCGAGAACCGCGCUUCAUUGCGGGUUUCAGGUCUCGUAACCAUGAAAUGGGUCCUGAUCUGCGUUCGUGGGAGCCAGAGGGAAGAGAGUCCCCAAGAACGUGGGGUUAUUUUAACUGCAUUCAAGCAGCGUUUUGCAGAACUUCCCGCACGAAAAAUCGAGACCAACACAGUUUGACUGCCUCCACGCCGGCCCAAAGCGUGCCCGCGCAGACAAGUGUUGAGUGGUCAAGAACUCGCUAUCAACGUCUACUUUCGUUAGCGGAUGUGUACGAGGACAGGCCACUGUGCAAUCAGGUAAGCGUUGCAAAUUCCAUUCUGAAUACCCGCGAGAGCGUAUUACGCAUCCACCGAAUAAUUAGGCAGAAAGUGCCCGCUUAUUCCAACGAGGGGGGUGACGUCGGCCUUCGGACAAUAGCAAAUGGCUCCGAUUCGACUCGGUGGAAACAGGUAAAUGCAUCAUAUAGAUAUGUCACAUCCCGCGCUGUCGGCCUUAGCCGCCCGCCAGCUUUUACCCUGAAUAUCCGGCAGUGUUGUGCGAGGGAGCUAUGCUCAUUUGGUACGACUAUCGCCCGCAAUCAGUGGAAGAUCCCCCGCUCGGAUCCUUAUACAACGCGAAACGACACCGGGCGGCGAUAUCGGUUUCAUCACGGUUGUCUGGAUAGUGGGUCGGCUCAGUCAUACAACCUUACCGGUAUUCAAUGGCUCAUGGAGCUCCCCGGGGGGUCAGCACCAUUCAUAACCCGGUAUGUUCAUUUGAGACUCAUCCGCUUGAAUCCUUAUAGUCUCGUGCGCAUGGAUACCAACUUGGAAGCCCGUUUAUAUGGAUGGAUGAACUUCGGCAAAACAACGUGGUCACUUACGAAGGUGGGGAUCGGUGGUGCUGCAAAUCCGGAGGGGACAAAUCCCGUAGCUACCGCCCACUAUCUUACGGUCCGCGAUGUCGUUUUUAAUCAAUUUGAGGAUCUCGAUACAAAUGCGUGGUAUAGUCGAUCUACCGAUGCUCUGCGUCUGUACCGGUCCCCAUUUGGGUCACAUUCAAUUCCUCGACCGGAAAUCUCCCCAAACAAUAAGGGGUUUGAGCAUGGAGACCACUCCUGCAACAUCCUGCAUCCCGUCUUCGGUUUAAGACCUAUACCGUCGGGCGGCAUGAACAGUGUUCAUAGGCGAGGAUCUGUGGGGACUAGGUGGGUCCGGCUAGAUCGGCGGAGUGGGUUCACACUUUGCUGGGCGUCUGGGCCUGCAGGCCGGUUCAACACGGGUAGGAAUCCUGGACUAUCACUGAAGACGGCCGUAAGGAUUUUCGAUGACCUGAAACAAAAGUUGAAUUUGGGGCAAUCUCAUGGCUGUUACGCCAAUCGAUGGGGGCCCCGUGUGCCUCAUCCAGGUCCCAUAAUAGGACUAAUGUUAGCGAAUGGAGAGAAUCGCUUACUAAUCAAGUUCGCACCGCGGUGUCCGCCUAUCUGGACGACGCGAUAUUUGCUGCUGCGUGAACCAGUAGAGUACAAAAGUACUGCGUGGAAACGAAUGGUGCCACCUGACGAUCCUGUAGUUCUCACCAUAAAAUCCACAAGAGGACUACAUCUUGAUAUAGGCCUGACUCACUGCUAUAUUCACGAAGUACAACAUAUGAUCGCGAACAGUCGUCAUUCGCGUCUAGAAUUCAACAAGCGACUAUCUCGUCUACCAGGUCUAGGGUGGCGACAUUAUACCGAGAUAGCGUAUAACCAACACUUGUGGAUACAUUCGCGUCUAAGGAGUCUUUGCGUCCAAAGGAUCUCGACUGCGGGGACCGGAGGCGAAAGCAUCCCGCGAGAAAUAGCCUGGGUUCCACAGCCAAGGUGUGUAGCCUUCGAAUUUUCGUACGGUUCGAAAACAUGUCACUAUCAGCGUUUCUACCAGCUUUCUAACCCAGCGUGGGUCGCGUAUGCCGCAUGUUCUAUAAGACGCCUACAUGUCAACGCACCUGCAUGUAAUGUAAAACGCAAUACUCAAUAUCAACAGAUCCCCAUUUCCGUUCUAGUAAAGGUGAAGUAUUACUUAUUCUCGUUUAUCCACAAGGUACGAGCUGACGUAAGGUGUAACUUAUGUCCCUAUCUGAGUUCGGCAAUUACGGUUUGUUACAUAUUCUGGCCGAUCAUAAAACAUAGCGACAGGGGGGACAGUAAAGGUCUAAGUGUGGACCGCGCUUUGACUUACAACGACUCAUCGCAGCCCCCUAGUAACCUAUCACACUCUCAGAUUGUGCGCAUAAGGUUAGGCUUAUGCCACAUUAGGGCCGGGAGGUUCACAGACGGGCAAACGCCAGAAGCUACCGGAGGCUCGGGUCGCCCGUGGCGUGCUCUAAGUGGAGUCCCCUACAGAACGACUAUAAGCCUCUUGCUCAGAACUAUUCCGAUCCAUCUAAAAAUUUCCAAUGUGGGAUCUGCAUUCCCUGACGAUUGGGACCGACCACUGAUAUUUGCUGCCGGUUUAUGUUUUUUGUCUGAUGGGUAUAUCGUCAGGCACGUAACCUACCACAUCGUAGGAGCUAGCCGCUUGUCACUCUUAAUGGUCACAGAAAUAUGUUUUAGAUCUAACCGCGAACAGUGUAAAAGCACCAACCCAAAUAUACCACAACGUCUUAGUGAUGCGAGCGGACGUACGAUAAUGCAGUAUCAGAGCUUAGGAGUUAUAAUGGAAUACCUUCUAAUCCAGUCGCCAUCCCAGUUGUUACCUGCGGGAACAACAACGCUGGUGCCUCGCGACGUCGGCUUUGAUGCAGAAACGUUGUAUCGGAGGGGGAUUGUCCGCCUCGUUGCCUCGGAGGAGAUUCAUGCGGGAACUACAUCGAUUCCUCGGCCGUCAUUAAACAUUAGCUCCCACCACGAGACAAAAUUCAAAAUGACUUGUUUGACUUACCAUCGCGAAAAGGACAGUCGCCAACACGAGGUGACCAUGUAUUUUCAUCACCGAUGUUCCCACCUCCCCACUAGACACUCCGCUGAUGUAGUAACCUUUGUACCAGUCUUAAACUUUGUUCUCUUAAAUCGUCGUGAGAUCAGCGAUAUCCUUAACGAAAAGUCUCCCUUGGAGAGAUGUCAACCCACGACCUCAUCUCCCUUUCAGUUCGGGCGUGCCCAGAGUAUAUUGCAUGCUCGGAAGGCUCAACGUACGCACAGAUCCUGGCGAGGCCAUUAUAUAAGGCGGGGCACUCGGUGCAAUACUCAACUUAAGCGCAAUUUAAUCCGUUCCUUUGAAGCCCACGAAAGCCAGAGAGUACAUUUGGAUCGUCGUUGCUCUGUAAUUGAUGGCGAUUCGAACACGCCAACCGACCUCCAUGGGGUGAAUAGCUUACCGUCGAUCAGAACCUUUAAAUACCUCAAGUCCAAAUAUCGCGUUUGGCUCCUUGUUGUACAAUGGGCCAGCGUAUGGGGUGCCAAGAUCUCUCUGGAUGGUGGAACCAAAAGCUCCAUGAUCUCCCGGCGUCUUUGCCGUCCUUACUAUUUGAGCACAAUGUUCGCUCUCAGCACGCAACGCUCCCACCUUCCGCUUGCUUGGCUUUAUCAAUAUGUCCGUUACAAGCAGAUGGUCAUAAUCGCCGAGAGAACGGGUGUCACAGUUCCGUUACGAAUUUUAUCCUGCGAAUCACUUUGGUGGAGACCGUCUUCCAUUGUAAUAAACGACACUAAGACCGGCGCCACGUAUGCCGGAAGGGAACCUAACAGUCUGGAAGGCCUUUCACCUGUCCUGGAGCGCCUAUUCGCACUGCUCACUUAUAUCGAAAGUACCCCUAGUGUAUUUGCCGGGCCCAUGUCACUCGAGGUCAAUAUUGGAUGCGACUCUGCAAAUCGUUGGCCGGACUGCGAGGGUGCACAAAUGCUGAUCGCUAAAAGCCUUGCGACUGCCUGUAAAGCGACGGGCACUCCUCUAGGACGGUUCCAGGUUGCUUGGGUCAAUAGGACAAGUGUGUCGGAAAGGACAGUUCCAUGGAGCGAAGAAUCAAGGAAGUCUAGUCCGUUACAUGCUAGUCGGUAUACCGCCUCCAAAGCGCUGACCCUGGGGCCGGCGAUAACGCAACGUAGACGAGAGCGCAUAAUGCGCGCCUGCUACGGACUACCUAUGUCGUCGCGGGUCGGUAGCCCACAACUAUCUCUGCUUAGUGCCCUAGUAUGCAGCCACUACUGCCUACAUGCUCAAAACUUCAUCCCCGCCGACACAUCGAGUGGCUGGAGGGCCGAAACCGUAGCGAGCUGUUUUAGAAGAGCCUUCUCUUGUCUCUUUCACACGGACAAGAGUUCGUACUGUCUGAGUGACAUAUUUUGCAUGAAAUUCAAGUUGAGAUCUUACAGACACUUUUUCGCCAAGGGCCGUGGGGACAUAGUCAAGGCCCGCUUAGCAACGCGAUGGCGCAUUCCGUGCAUCGUGCGGUUCGCGGCUGACUGCAUUCUCGCCGUCGUUAGGCUAUGCGGUGAAAUCUACAUCGAACCUGGGCAACUUCACACGGCGCCUCGUGAGGUUACACCACUUGAGGAGUCACCCAUGUCGUCGCAUUCGCCGGCAUGGGAGAGAUAUAGAAUUGACUGUCCAACUCGAACGAGCUUUGAAAAGCUACGGUGGGAAUUUAACACUCGACACAGUUAUAAGCGCUUCCUGUUUUCCCCUUGGACGCAGGGCAUCCGACGCGAGAAGGUGAGACCAUGUCUAAAUACCACGGGGGAGUCAGGAAUUCCAGACCGAGGGAGUGGAGCGACGCUGUCCCGGCAAAAGCUGCCUAAGAACUGUCGAUUUAGUGCCUUGCACAGGAGUCGCUGGUUAAUAGUCGGCGCCACUAUGUGGCGAAUGGCGGGAAUCAAUCGAAGGAAUACGAGCCGCUUCCGAGAAGGACACGCGCCACGAUCGAUCCCAGGAGCCGAACGCUCCGUUAGCUCAGAAUGUUGUACUAUUACGAACUCAUUAUGCUUGGGCCUUAUGUCCAGCUCCGUGCGACCGAUGCGGAGCAAUGUCUGUUGUGCUCUAUGGCGACAGUUAUGUUCGGCCGGGCUGUUGGCCCGUACUUUGGGUUUAGUUGUUGAAGUUGGGUCGCCCGGUGUUCAUCGCCCGAUCAGGUUUCCUCUCACUUCCACCCAGAAGGACUGCUCUGAAUGUUGUUGUGCCGACUGCGGAAAAGGCAGCACCGUAGGUGCAGCAAGAAGCACGCGCCUCAAGAGACAUCUGGCUGAAGUACACCGGGUGCACCAGCCCCCUGGAGUACAACUGCGCUCGGCCCUCACAUCAGUAAUCGCUCCAGUAACACCUAUACACCGCAAACGGGCCCGGGUAGCCCGCACCACGAAUAGUUUGAACACUGAGAUCUAUGUUCAAUUACGCAGUCUUAAACACUUGCAGUUGAUUUGGCGAGUCACGGACCAUCUCAGGACAUGCCAAUGGCUUGCUAUAUUUCCCGCGAGCCUAGGCUUUAUAUCAUCCCGCCGGCAGAGUCAGCAUAAUUCUGGAGAGGUGCCUUGCGCAUACAAUAGUCGAUGCCAGAGCUUGGGCUGUGAGCACUUUAACGUGAAUAACCGGGAAUUACAAUCCGAAACUACUGCCAGAGCCAUUCAUAUUUUCCUGGCGUUAUACUAUUUACCCCUAUUUAGUCGGACGAAUCGCUAUCAAGUAUGUAGAGCGCCGAGUGUUUGCAAUCUUACUCGAGCCAUGCGUGUAGGUCCCUCUAUUACGGCCGUCUUCGAUGUGAAGAAACACGGUCACCGUAGGUGGGAAAGAAAUGCAUCGAUGAACCGACGGUGCACCCAACUGGAACUCUCCUGUACAGCUCUCCGUUACGUGCUAUCCUUACUUGAUCGCACAUUGCCGCGCAGUUCAUUUUCCCCAACUCGAUGGCGGAGCGUCGUGACUGCACUUGCGCGGGUAGGGCAAGUGGUAAGUGGACGAGCCGGACGUGCCCAUUUAUGCGGGGUCAUCCGGUCAAUCAGAUACUGCCAGCGUGGCUAUAGCACACCAAUCUAUCCCAAGGUAGCGCUCCGUAUGGCCGGUCACAGGACACUCAGCCCAACAGGUUCUACAUGUUGCUGUGCGCGCAACUCCCCUCCUCCGGUAUCAAAAGGUAACCUGCAUUAUAUCCCGUCACACUUUGGUUUCGAUCUAGUUUUCGAAUCGAAGCCAUCCUUCGACAGGUCAGCUCACUUACACCAUGCAAAUAGGCCCAACUCUCUCAGAAAUGGGGACCAUAAAAACCCCGAGAGUGAUUGCAUCGGCCACGCGUUCCAUCUAGCUAAUUAUGUGCACCAAAGAGGUCUAAUGAGCAAGAGGAACGUAAGUCAUUAUUGGACAGUUCCAUUGAGAGGCCUUGGUACAGUAGUGUCUUUGCAUGUGAGGGACCAGAACUUUCGCAGUCGAGCUAUUCAGAGCGUGGUUGUAGUGCGUACGUAUUCCCACGCCUUGAUGCCGAUAGAGAUCAAGAAACCCCCAUACGCUAAUGGCAGAGGUAAUAAGGAUACAAUUUAUAGAAGAUAUGCGCAGACGCCUUGUCCUCCCCGCUAUACGAGACUGCCAUUGCAUACCAGUAGACGAAAAGCCAGAGGUCUUGACCUUUCCACCUACCGGUACUCUCCAUUAUAUUUUGGAUUGCCGGGCCAGUGGGAAGCCUGGCGUAAGCUCAAUUGUCCGGUUAGUCCGCGCGCUGCCGGCUACCUACGUGCUCCAUUUCGGUGUAUUGAAGCAGACCCUCUUGGACUAAGGCUGCUAAUGGGGAGGUAUCGGUACGUACGCGAUUAUAAAUCAUUCAGGUCCGCUAAACUGUGUGUAAUCAGUACGGUUCACCGAUGCUCAAAGAUUAACCAGCAGAGAACUAUUACGAGCAUCAUUGCUCAGACGCCGGUAUCCGAUAAAUCGAAUACGGUGAGUCUACUUGAGCUCCUGCGAGCAUCUGCUCUGCCAUUAUUUUGCUCCUGCAGUCCAGCUGUAAACAAAGGCCUAUUCGCGACUAAAUUGGUUGUUAAAUAUUCCAAAGCAUGCCUGUGGCGCGCAGCCAGCGACAGUAAGGGGAGGGGGCUAUCCCCAUUUGAACCGAGCCAAAGAGUUGAUAAAGAACAAGACUUUGUUGUUAAGCCGAAAUUGAAUUUUGGGUGUGCAUGCUGGCAUUACUGCGCCUCGUCAGUGCUGCAGCCUUGGAAUGAGGGUGCAGGACUGGUCAGAGACCCGGUGAUUUUUCGCUAUUGCACUAGUGAAACUGCCGACAGGGCUUCGUACAUUAAAUCAGUGACAGACUUACUUCUGAAAGUAUGCUCCAGGGGAUCACAGCCCCUAACUAACAAAUCGGGACUAGUCGCCGCGGGAUUUUCUGCGGAUGGGGCCGCCGGCGUCUUCAACGCUAACCGGCCGACGCCGCGGAACCCGCCGGCGGCGAUCCAGUCUGCGGGUAUGCCCUCAUCAGGUGGUUACCGUGUUCCUUAG'

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

    print(translation(sample_rna))