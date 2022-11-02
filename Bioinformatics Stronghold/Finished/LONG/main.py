import sys
sys.path.append('../..')
import molbio as mb

def contig(strands: list) -> str:
    superstrand = strands.pop(0)
    while len(strands) > 0:
        for strand in strands:
            for l in [len(strand) - i for i in range(len(strand)//2+1)]:
                str_left = strand[:l]
                sup_left = superstrand[:l]
                str_right = strand[len(strand)-l:]
                sup_right = superstrand[len(superstrand)-l:]
                if str_left == sup_right:  # Strand-right = Superstrand-left
                    superstrand = superstrand + strand[l:]
                    strands.remove(strand)
                    break
                if sup_left == str_right: # Strand-left = Superstrand-right
                    superstrand = strand[:len(strand)-l] + superstrand
                    strands.remove(strand)
                    break

    return superstrand


if __name__ == '__main__':
    strands = [*mb.parse_fasta('strands.fasta').values()]
    print(contig(strands))
