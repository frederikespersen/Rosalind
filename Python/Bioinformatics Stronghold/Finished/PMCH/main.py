import sys
sys.path.append('../..')
import molbio as mb


def perfect_matchings(strand: str) -> int:
    a = sum([1 for nb in strand if nb == 'A'])
    c = sum([1 for nb in strand if nb == 'C'])
    g = sum([1 for nb in strand if nb == 'G'])
    u = sum([1 for nb in strand if nb == 'U'])

    assert a == u
    assert c == g

    from math import factorial as f
    perf_matchs = f(g) * f(a)
    return perf_matchs


if __name__ == '__main__':
    strand = [*mb.parse_fasta('PMCH/strand.fasta').values()][0]
    n = perfect_matchings(strand)
    print(n)

