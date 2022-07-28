#!/usr/bin/env python3

import numpy as np
import os
import tskit


if __name__ == "__main__":
    ts = tskit.load(snakemake.input["ts"])
    vcf_path = snakemake.output["vcf"]
    with open(vcf_path, "w") as vcf:
        ts.write_vcf(vcf, contig_id=snakemake.wildcards["chr"], position_transform="legacy")

    # Get last record of newly created VCF
    with open(vcf_path, "rb") as vcf:
        try:
            vcf.seek(-2, os.SEEK_END)
            while vcf.read(1) != b'\n':
                vcf.seek(-2, os.SEEK_CUR)
        except OSError:
            vcf.seek(0)
        last_record = vcf.readline().decode()

    # If last record is not at end of chromosome, add a fake monomorphic record at the end.
    # This makes simulation downstream easier to deal with.
    fields = last_record.split("\t")
    pos = int(fields[1])
    n = len(fields) - 9
    m = snakemake.params["chromosome_length"]
    if pos < m:
        fields[1] = str(m)
        fields[3] = fields[4] = "N"
        fields[9:] = ["0|0" for _ in range(n)]
        last_record = "\t".join(fields)            
        with open(vcf_path, "a") as vcf:
            print(last_record, file=vcf)
