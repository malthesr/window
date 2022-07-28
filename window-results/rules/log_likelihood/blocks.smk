"""
Rules for getting the log-likelihood of the individual winsfs blocks
"""

BLOCK_DIR = WINSFS_DIR / "blocks"


rule winsfs_tolerance:
    """ Run winsfs in memory with block log-likelihood loggin """
    input:
        winsfs=rules.install_winsfs.output.winsfs,
        safs=get_safs,
    output:
        sfs=BLOCK_DIR / "{pops}_b{b}_w{w}.sfs",
    params:
        max_epochs=50,
        tolerance=1e-100, # Trigger logging, not stopping
    log:
        BLOCK_DIR / "{pops}_b{b}_w{w}.log",
    threads: 20
    shell:
        """
        {input.winsfs} \
            -vvv \
            --seed 1 \
            --threads {threads} \
            --max-epochs {params.max_epochs} \
            --tolerance {params.tolerance} \
            --blocks {wildcards.b} \
            --window-size {wildcards.w} \
            {input.safs} \
            > {output.sfs} \
            2> {log}
        """


rule tidy_block_log_likelihoods:
    """ Collect epoch sum of block log-likelihoods from log """
    input:
        sfs=rules.winsfs_tolerance.output.sfs, # Ensure finished log
        log=rules.winsfs_tolerance.log,
    output:
        csv=BLOCK_DIR / "{pops}_b{b}_w{w}.block_loglik.csv",
    params:
        grep="^DEBUG \[stop\] Current log-likelihood -[0-9\.e]+, Δ=[0-9-e\.]+",
        awk="""BEGIN{OFS=","} {print NR+1,$0}""",
    shell:
        """
        ( \
            echo "epoch,loglik,change" &&
            grep -Po '{params.grep}' {input.log} | \
                grep -Po '[0-9].*$' | \
                sed 's/, Δ=/,/' | \
                awk '{params.awk}' \
        ) > {output.csv}
        """


def get_gather_block_log_likelihoods_input(wc):
    """Helper function to get tidy block log-likelihood CSVs"""
    name = rules.tidy_block_log_likelihoods.output.csv
    lst = expand(name, pops=wc.pops, b=config["run"]["b"], w=config["run"]["w"])
    return lst


rule gather_block_log_likelihoods:
    """ 
    Gather block log-likelihood CSVs across methods and parameters for populations 
    """
    input:
        get_gather_block_log_likelihoods_input,
    output:
        csv=RESULTS_DIR / "{pops}.block_loglik.csv",
    params:
        conditions=expand("winsfs_b{b}_w{w}", b=config["run"]["b"], w=config["run"]["w"]),
        fmt="idd",
    script:
        f"{workflow.basedir}/scripts/gather.R"
