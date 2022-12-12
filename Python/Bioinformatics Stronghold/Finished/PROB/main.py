def log_rsp(strand:str, gc: float) -> float:
    """
    Logarithmic random sampling probability
    """
    pGC = gc
    pAT = 1 - pGC
    ps = {'A': pAT / 2,
          'C': pGC / 2,
          'G': pGC / 2,
          'T': pAT / 2}

    import math
    p = 1
    for nb in strand:
        p *= ps[nb]

    pp = math.log(p, 10)
    return pp


strand_input = 'TCCAGTTGCGAGCTCGACGCCAAACAAATCCGCCATCGGTTGCACCCTAAGATATAATGAGCAACTCCCGAATTGACGCTTAAATAAACTGGACAGTAGT'
gcs_input = '0.069 0.168 0.199 0.261 0.342 0.402 0.442 0.485 0.562 0.628 0.660 0.758 0.786 0.871 0.925'
if __name__ == '__main__':
    pps = [log_rsp(strand_input, float(gc)) for gc in gcs_input.split()]
    print(*pps)