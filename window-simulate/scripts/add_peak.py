#!/usr/bin/env python3

import sys

def generate_gts(alt):
    alleles = [0 for _ in range(40)]
    for i in range(alt):
        alleles[i] = 1
    gts = [f"{i}|{j}" for i,j in zip(alleles[::2], alleles[1::2])]
    return ":".join(gts)

# Makes a little smiley face :)
PEAKS = [
    f"{generate_gts(19)}:{generate_gts(18)}",
    f"{generate_gts(18)}:{generate_gts(19)}",
    f"{generate_gts(21)}:{generate_gts(19)}",
    f"{generate_gts(18)}:{generate_gts(20)}",
    f"{generate_gts(18)}:{generate_gts(21)}",
    f"{generate_gts(21)}:{generate_gts(22)}",
    f"{generate_gts(19)}:{generate_gts(23)}",
]

if __name__ == "__main__":
    chromosome_length = int(sys.argv[1])
    peak_height = int(sys.argv[2])

    i = 0
    peak = 0
    n = 0

    for line in sys.stdin:
        chrom, pos, gts = line.split()
        if peak < len(PEAKS):
            while True:
                i += 1
                if i >= chromosome_length:
                    i = 0            
                    chrom = str(int(chrom) + 1)
                if i >= int(pos):
                    break
                if n >= peak_height:
                    n = 0
                    peak += 1
                    if peak >= len(PEAKS):
                        break
                n += 1
                print(chrom, i, PEAKS[peak], sep="\t", file=sys.stdout)
        print(chrom, pos, gts, sep="\t", file=sys.stdout)
          
