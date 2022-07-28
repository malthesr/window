"""
Rules for getting the log-likelihood of the individual SFS iterations
"""


if "simulation" in config:
    rule nth_log_likelihood:
        """ Estimate log-likelihood after n full epochs based on true SFS """
        input:
            sfs=Path("{dir}") / "{pops}{maybe}.after{n}.sfs",
        output:
            log_likelihood=Path("{dir}") / "{pops}{maybe}.after{n}.{split}.loglik",
        params:
            true_sfs=lambda wc: config["simulation"]["truth"][wc.pops]
        script:
            f"{workflow.basedir}/scripts/sim_log_likelihood.R"
else:
    rule nth_log_likelihood:
        """ Estimate log-likelihood after n full epochs based on data set """
        input:
            winsfs=rules.install_winsfs.output.winsfs,
            safs=lambda wc: get_safs(wc, split=wc.split),
            sfs=Path("{dir}") / "{pops}{maybe}.after{n}.sfs",
        output:
            log_likelihood=Path("{dir}") / "{pops}{maybe}.after{n}.{split}.loglik",
        shell:
            """
            {input.winsfs} log-likelihood --sfs {input.sfs} {input.safs} > {output.log_likelihood}
            """


def get_log_likelihoods(dir, pops, maybe, epochs, split):
    """Helper to get log-likelihoods for set of conditions"""
    return expand(
        rules.nth_log_likelihood.output.log_likelihood,
        dir=dir,
        pops=pops,
        maybe=maybe,
        n=epochs,
        split=split,
    )


def get_all_log_likelihoods(wc):
    """Find test and train log-likelihood after all epochs depending on wildcards"""
    epochs = get_epochs(wc)
    loglik_dir = Path(wc.dir) / "iterations"
    train_loglik = get_log_likelihoods(loglik_dir, wc.pops, wc.maybe, epochs, "train")
    if "simulation" in config:
        return {"train_loglik": train_loglik}
    else:
        test_loglik = get_log_likelihoods(loglik_dir, wc.pops, wc.maybe, epochs, "test")
        return {"train_loglik": train_loglik, "test_loglik": test_loglik}


rule tidy_log_likelihoods:
    """ 
    Tidy (test, train) log-likelihoods for single (pair of) 
    populations into tidy CSV 
    """
    input:
        unpack(get_all_log_likelihoods),
    output:
        csv=Path("{dir}") / "{pops}{maybe}.loglik.csv",
    params:
        pops=lambda wc: wc.pops.split("-"),
        epochs=lambda wc: get_epochs(wc),
    wildcard_constraints:
        dir=f"({WINSFS_DIR}|{WINSFS_STREAM_DIR}|{REALSFS_DIR})",
    script:
        f"{workflow.basedir}/scripts/tidy_log_likelihoods.R"


rule gather_log_likelihoods:
    """
    Gather log likelihood CSVs across methods and parameters for populations
    """
    input:
        get_gather_input,
    output:
        csv=RESULTS_DIR / "{pops}.{csv}.csv",
    params:
        conditions=get_conditions(),
        fmt="idd",
    wildcard_constraints:
        csv="loglik",  # Input rule requires wildcard be present
    script:
        f"{workflow.basedir}/scripts/gather.R"
