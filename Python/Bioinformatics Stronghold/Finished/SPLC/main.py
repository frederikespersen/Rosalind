import sys
sys.path.append('../..')
import molbio as mb


def splice_introns(strand: str, introns: list) -> str:
    splice = strand
    for intron in introns:
        splice = splice.replace(intron, '')
    return splice


fasta_filename = 'exons_introns.fasta'
if __name__ == '__main__':
    strings = mb.parse_fasta(fasta_filename)
    coding_strand = [*strings.values()][0]
    introns = [*strings.values()][1:]
    exons = splice_introns(coding_strand, introns)
    mrna = mb.transcribe(exons, template=False)
    orfs = mb.orfs(mrna)
    peptide = mb.translate(orfs[0])
    print(peptide)