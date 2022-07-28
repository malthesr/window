"""
Simulate demography using msprime for independent chromosomes using demes input,
and concatenate independent simulations into BCF.
"""

MODEL_DIR = RESULTS_DIR / "model"
TS_DIR = RESULTS_DIR / "ts"
BCF_DIR = RESULTS_DIR / "bcf"


wildcard_constraints:
    chr="[0-9]+",


rule create_model:
    """
    Creates a serial bottleneck model with a migration and outputs it in demes YAML format.
    The model is plotted as a "tubes" plot and a "sizes" plot for overview and debugging.
    """
    output:
        yaml=MODEL_DIR / "model.yaml",
        model_fig=MODEL_DIR / "model.pdf",
        sizes_fig=MODEL_DIR / "sizes.pdf",
    params:
        ids=IDS,
        migrations=config["migrations"],
    script:
        f"{workflow.basedir}/scripts/create_model.py"


rule simulate_demography:
    """
    Simulate a single chromosome using msprime, overlaying mutations.
    """
    input:
        yaml=rules.create_model.output.yaml,
    output:
        ts=temp(TS_DIR / "n{n}.chr{chr}.ts"),
    params:
        mutation_rate=config["mutation_rate"],
        recombination_rate=config["recombination_rate"],
        chromosome_length=config["chromosome_length"],
        sample_sizes=lambda wildcards: get_all_pop_n(wildcards),
        seed=lambda wildcards: int(wildcards.chr),
    script:
        f"{workflow.basedir}/scripts/simulate_demography.py"


rule ts_to_vcf:
    """
    Convert chromosome tree-sequence file to VCF with appropriate chromosome ID.
    """
    input:
        ts=rules.simulate_demography.output.ts,
    output:
        vcf=BCF_DIR / "n{n}.chr{chr}.vcf",
    params:
        chromosome_length=config["chromosome_length"],
    script:
        f"{workflow.basedir}/scripts/ts_to_vcf.py"


rule concat_vcfs:
    """
    Concats VCF files for all chromosomes and removes all multiallelic sites.
    """
    input:
        vcfs=expand(
            rules.ts_to_vcf.output.vcf,
            chr=range(1, config["chromosomes"] + 1),
            allow_missing=True,
        ),
    output:
        bcf=BCF_DIR / "n{n}.bcf.gz",
    threads: 4
    shell:
        """
        ( \
        bcftools concat -O u --threads {threads} {input.vcfs} | \
            bcftools view -M 2 -O b --threads {threads} > {output.bcf} \
        ) 2> /dev/null
        """
