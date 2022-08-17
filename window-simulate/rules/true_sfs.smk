"""
Process demographic simulations to get true counts in text format.
"""

TRUTH_DIR = RESULTS_DIR / "truth"


rule create_samples_file:
    """
    Creates a file containing the names of samples in the VCF files corresponding to a particular
    population ID from the original model.
    """
    input:
        yaml=rules.create_model.output.yaml,
    output:
        samples_file=TRUTH_DIR / "{id}_n{n}.samples_file",
    run:
        start, end = get_pop_range(wildcards)
        with open(output["samples_file"], "w") as samples_file:
            for i in range(start, end):
                print(f"tsk_{i}", file=samples_file)


rule vcf_to_counts:
    """
    Get counts from true genotypes and output in text format.
    """
    input:
        bcf=rules.concat_vcfs.output.bcf,
        samples_file=rules.create_samples_file.output.samples_file,
    output:
        counts=TRUTH_DIR / "{id}_n{n}.counts.txt",
    params:
        bcftools_query=r"'%INFO/AC\n'",
    threads: 4
    shell:
        """
        bcftools view \
            -O b \
            -S {input.samples_file} \
            --threads {threads} \
            {input.bcf} | \
        bcftools query -f {params.bcftools_query} > {output.counts}
        """


def get_counts(wc):
    """Helper to get plain counts for 2d SFS creation"""
    fst, snd = wc.pair.split("-")
    counts = expand(TRUTH_DIR / "{id}_n{n}{peak}.counts.txt", id=[fst, snd], **wc)
    return counts


rule true_sfs:
    """
    Calculates global 2D SFS from true counts
    """
    input:
        counts=get_counts,
    output:
        sfs=TRUTH_DIR / "{pair}_n{n}{peak}.sfs",
    params:
        n=lambda wc: int(wc.n),
        sites_considered=config["chromosomes"] * config["chromosome_length"],
    script:
        f"{workflow.basedir}/scripts/true_sfs.py"
