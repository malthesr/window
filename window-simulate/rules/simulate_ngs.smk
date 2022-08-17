"""
Simulate sequencing errors and output as joint GLF file.
"""

GLF_DIR = RESULTS_DIR / "glf"
SAF_DIR = RESULTS_DIR / "saf"


rule simulate_gl:
    """
    Simulate genotype likelihoods from genotypes; depths and error rates cycle to the number of
    populations as necessary.
    """
    input:
        simgl="bin/simgl",
        bcf=rules.concat_vcfs.output.bcf,
    output:
        glf=temp(GLF_DIR / f"{BASE}.glf.gz"),
    params:
        error_rate=lambda wildcards: int(wildcards.e) / 1_000,
        mean_depth=lambda wildcards: wildcards.d,
        query=r"'%CHROM\t%POS\t[%GT:]\n'",
    shell:
        """
        bcftools query --format {params.query} {input.bcf} | \
        sed 's/:$//' | \
        {input.simgl} \
            --error-rates {params.error_rate} \
            --interleave \
            --mean-depths {params.mean_depth} \
            --seed 1 \
            > {output.glf}
        """


rule split_glf:
    """
    Split the full GLF file by population.
    """
    input:
        subsetglf="bin/subsetglf",
        glf=GLF_DIR / (str(BASE) + "{peak}.glf.gz"),
    output:
        glf=temp(GLF_DIR / (str(ID_BASE) + "{peak}.glf.gz")),
    params:
        subset=lambda wc: ",".join([str(n) for n in range(*get_pop_range(wc))]),
        total_n=lambda wc: get_total_n(wc),
    shell:
        """
        {input.subsetglf} \
            {params.subset} \
            {params.total_n} \
            {input.glf} \
            > {output.glf}
        """


rule create_reference:
    """
    Create fake reference for ANGSD analysis. Note that despite having simulated different
    chromosomes, we create a reference with only one chromosome, since ANGSD treats all sites in
    GLF as lying on the same chromosome.
    """
    output:
        ref=SAF_DIR / "ref.fa.gz",
    params:
        uncompressed=lambda wc, output: output.ref.removesuffix(".gz"),
        sequence_length=config["chromosome_length"] * config["chromosomes"],
    shell:
        """
        echo "> chr1" > {params.uncompressed};
        # Naive bash solutions slow, see https://stackoverflow.com/a/30288267
        perl -E 'say "A" x {params.sequence_length}' >> {params.uncompressed};
        bgzip {params.uncompressed}
        """


rule index_reference:
    """
    Create fai index reference. Note that ANGSD is very finicky about the time stamps on the index
    file, so its best to keep this in a separate rule from the index creation.
    """
    input:
        ref=rules.create_reference.output.ref,
    output:
        fai=SAF_DIR / "ref.fa.gz.fai",
        gzi=SAF_DIR / "ref.fa.gz.gzi",
    shell:
        """
        sleep 1; # ensure different timestamp from reference (enforced by ANGSD...)
        samtools faidx {input.ref}
        """


rule create_saf:
    """ Calculate SAF likelihoods from genotype likelihoods """
    input:
        glf=rules.split_glf.output.glf,
        ref=rules.create_reference.output.ref,
        fai=rules.index_reference.output.fai,
    output:
        saf_files=multiext(
            str(SAF_DIR / (str(ID_BASE) + "{peak}")),
            ".saf.idx",
            ".saf.pos.gz",
            ".saf.gz",
            ".arg",
        ),
    log:
        SAF_DIR / (ID_BASE + "{peak}.log"),
    params:
        prefix=lambda wc, output: output.saf_files[0].removesuffix(".saf.idx"),
    shell:
        """
        angsd \
            -anc {input.ref} \
            -dosaf 1 \
            -fai {input.fai} \
            -glf {input.glf} \
            -nind {wildcards.n} \
            -out {params.prefix} \
            2> {log};
        # ANGSD returns exit code 0 even when unsuccesful (...),
        # so we grep the log to produce a better exit code
        tail -n 1 {log} | grep -Piq '^[\s]+\[all done\]';
        """
