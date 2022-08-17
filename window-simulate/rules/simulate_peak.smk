"""
Simulate sequencing errors in genotypes with spike-in peaks
"""

rule add_peak:
    """
    Add peak to simulated genotypes in plain text format.
    """
    input:
        simgl="bin/simgl",
        bcf=expand(rules.concat_vcfs.output.bcf, n=20, allow_missing=True),
    output:
        txt=GLF_DIR / "n20.peak.txt",
    params:
        query=r"'%CHROM\t%POS\t[%GT:]\n'",
        chromosome_length=config["chromosome_length"],
        peak_height=10_000,
    shell:
        """
        bcftools query --format {params.query} {input.bcf} | \
        sed 's/:$//' | \
        scripts/add_peak.py {params.chromosome_length} {params.peak_height} \
            > {output.txt}
        """
        

rule peak_txt_to_counts:
    """ Simplify plain text counts to format for making true SFS """
    input:
        txt=rules.add_peak.output.txt,
    output:
        counts=TRUTH_DIR / "{id}_n20.peak.counts.txt",
    params:
        cut_fields=lambda wc: "1-20" if wc.id == "A" else "21-40"
    shell:
        """
        cat {input.txt} | \
            cut -f3 | \
            cut -d ':' -f {params.cut_fields} | \
            sed 's/|\|:/+/g' | \
            bc > {output.counts}
        """


rule simulate_peak_gl:
    """
    Simulate genotype likelihoods from genotypes with added peaks
    """
    input:
        simgl="bin/simgl",
        txt=rules.add_peak.output.txt,
    output:
        glf=GLF_DIR / "n20_d{d}_e{e}ppt.peak.glf.gz",
    params:
        error_rate=lambda wildcards: int(wildcards.e) / 1_000,
        mean_depth=lambda wildcards: wildcards.d,
    shell:
        """
        cat {input.txt} | \
        {input.simgl} \
            --error-rates {params.error_rate} \
            --interleave \
            --mean-depths {params.mean_depth} \
            --seed 1 \
            > {output.glf}
        """