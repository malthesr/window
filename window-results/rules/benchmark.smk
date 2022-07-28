rule bench_realsfs:
    input:
        realsfs=rules.install_realsfs.output.realsfs,
        safs=get_safs,
    output:
        sfs=BENCH_DIR / "{pops}_i{i}.realsfs.raw.sfs",
    params:
        tolerance=1e-100,
    benchmark:
        BENCH_DIR / "{pops}_i{i}.realsfs.benchmark"
    threads: 20
    shell:
        """
        {input.realsfs} \
            -seed 1 \
            -m 0 \
            -P {threads} \
            -tole {params.tolerance} \
            -maxIter {wildcards.i} \
            {input.safs} \
            2> /dev/null \
            > {output.sfs}
        """


rule bench_winsfs:
    input:
        winsfs=rules.install_winsfs.output.winsfs,
        safs=get_safs,
    output:
        sfs=BENCH_DIR / "{pops}_i{i}.winsfs.sfs",
    benchmark:
        BENCH_DIR / "{pops}_i{i}.winsfs.benchmark"
    threads: 20
    shell:
        """
        {input.winsfs} \
            --seed 1 \
            --threads {threads} \
            --max-epochs {wildcards.i} \
            --blocks 500 \
            --window-size 100 \
            {input.safs} \
            > {output.sfs}
        """


rule bench_stream_winsfs:
    input:
        winsfs=rules.install_winsfs.output.winsfs,
        saf_shuf=rules.shuffle.output.saf_shuf,
    output:
        sfs=BENCH_DIR / "{pops}_i{i}.stream_winsfs.sfs",
    benchmark:
        BENCH_DIR / "{pops}_i{i}.stream_winsfs.benchmark"
    threads: 1
    shell:
        """
        {input.winsfs} \
            --seed 1 \
            --max-epochs {wildcards.i} \
            --blocks 500 \
            --window-size 100 \
            {input.saf_shuf} \
            > {output.sfs} \
            2> /dev/null
        """


rule gather_bench:
    """
    Gather benchmarking results into tidy CSV.
    """
    input:
        expand(rules.bench_realsfs.benchmark, i=[1] + list(range(10, 101, 10)), allow_missing=True),
        expand(rules.bench_winsfs.benchmark, i=range(1, 11), allow_missing=True),
        expand(rules.bench_stream_winsfs.benchmark, i=range(1, 11), allow_missing=True),
        expand(rules.shuffle.benchmark, allow_missing=True)
    output:
        csv=BENCH_DIR / "{pops}.bench.csv"
    script:
        f"{workflow.basedir}/scripts/gather_bench.R"


rule all_bench:
    """ Phony rule to expand benchmarking results """
    input:
        expand(rules.gather_bench.output.csv, pops=config["run"]["pops"])
