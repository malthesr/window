#!/usr/bin/env python3

import numpy as np


def read_counts(path):
    """ Read counts from plain text file for a single population """
    counts = np.loadtxt(path, dtype = "uint8")
    return counts


def calc_sfs(c1, c2, n, sites_considered=None):
    """ Create 2d SFS from counts for two populations """
    sfs = np.zeros([2 * n + 1 for _ in range(2)], dtype="uint64")
    for i1, i2 in zip(c1, c2):
        sfs[i1, i2] += 1
    if sites_considered:
        sfs[0, 0] += sites_considered - np.sum(sfs)
    return sfs


def format_sfs(sfs):
    """
    Write SFS in same format as winsfs (i.e. ANGSD with header): 
    first line is header giving the shape, next line is SFS in flat, C-major order
    """
    fmt_dim = "/".join(str(x) for x in sfs.shape)
    header = f"#SHAPE=<{fmt_dim}>"
    fmt_sfs = [f"{x:.0f}" for x in sfs.flatten(order="C")]
    flat_sfs = ' '.join(fmt_sfs)
    return f"{header}\n{flat_sfs}"
    

if __name__ == "__main__":
    paths = snakemake.input["counts"]
    assert len(paths) == 2, "Must provide exactly two paths"
    n = snakemake.params["n"]
    sites_considered = snakemake.params["sites_considered"]

    c1, c2 = map(read_counts, paths)
    sfs = calc_sfs(c1, c2, n, sites_considered)
    assert np.sum(sfs) == sites_considered, "Number of sites doesn't fit"
    
    with open(snakemake.output["sfs"], "w") as f:
        print(format_sfs(sfs), file=f)
