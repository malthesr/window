""" Rules for tidying and gathering SFS results """

TRUTH_DIR = RESULTS_DIR / "truth"


def get_all_sfs(wc):
    """Find SFS after all epochs depending on wildcards"""
    epochs = get_epochs(wc)
    if str(wc.dir) == str(TRUTH_DIR):
        sfs = config["simulation"]["truth"][wc.pops]
    elif str(wc.dir) == str(REALSFS_DIR):
        sfs = expand(rules.nth_realsfs.output.sfs, **wc, n=epochs)
    else:
        maybe_pattern = re.compile("_b(?P<b>[0-9]+)_w(?P<w>[0-9]+)")
        params = maybe_pattern.match(wc.maybe).groupdict()
        if str(wc.dir) == str(WINSFS_DIR):
            sfs = expand(rules.nth_winsfs.output.sfs, **wc, **params, n=epochs)
        elif str(wc.dir) == str(WINSFS_STREAM_DIR):
            sfs = expand(rules.nth_stream_winsfs.output.sfs, **wc, **params, n=epochs)
    return sfs


rule tidy_sfs:
    """ 
    Tidy spectra and stats for single (pair of) populations into tidy CSV 
    """
    input:
        sfs=get_all_sfs,
    output:
        sfs_csv=Path("{dir}") / "{pops}{maybe}.sfs.csv",
        stat_csv=Path("{dir}") / "{pops}{maybe}.stat.csv",
    params:
        pops=lambda wc: wc.pops.split("-"),
        epochs=lambda wc: get_epochs(wc),
    wildcard_constraints:
        dir=f"({WINSFS_DIR}|{WINSFS_STREAM_DIR}|{REALSFS_DIR}|{TRUTH_DIR})",
    script:
        f"{workflow.basedir}/scripts/tidy_sfs.R"


def get_gather_sfs_input(wc):
    """
    Wrapper around get_gather_input, but including the truth if simulated.
    """
    input = get_gather_input(wc)
    if "simulation" in config:
        truth = expand(TRUTH_DIR / "{pops}.{csv}.csv", **wc)
        input += truth
    return input


def get_sfs_conditions():
    """
    Wrapper around get_conditions, but including the truth if simulated.
    """
    conditions = get_conditions()
    if "simulation" in config:
        conditions += ["truth"]
    return conditions


rule gather_sfs:
    """ Gather SFS/stat CSVs across methods and parameters for populations """
    input:
        get_gather_sfs_input,
    output:
        csv=RESULTS_DIR / "{pops}.{csv}.csv",
    params:
        conditions=get_sfs_conditions(),
        fmt=lambda wc: "id" if 1 == len(wc.pops.split("-")) else "idd",
    wildcard_constraints:
        csv="(sfs|stat)",
    script:
        f"{workflow.basedir}/scripts/gather.R"
