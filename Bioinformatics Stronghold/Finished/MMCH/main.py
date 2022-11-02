import sys
sys.path.append('../..')
import molbio as mb


def maximum_matchings(strand: str) -> int:
    a = sum([1 for nb in strand if nb == 'A'])
    c = sum([1 for nb in strand if nb == 'C'])
    g = sum([1 for nb in strand if nb == 'G'])
    u = sum([1 for nb in strand if nb == 'U'])

    au_max = max(a, u)
    cg_max = max(c, g)
    au_min = min(a, u)
    cg_min = min(c, g)

    # Permutation: Pick from xx_max for a total of xx_min times without repetition
    from math import factorial as f
    matchings_au = f(au_max) // f(au_max - au_min)
    matchings_gc = f(cg_max) // f(cg_max - cg_min)
    max_matchs = matchings_au * matchings_gc

    return max_matchs


if __name__ == '__main__':
    strand = [*mb.parse_fasta('strand.fasta').values()][0]
    n = maximum_matchings(strand)
    print(n)
