"""
General helper functions and input functions used in several places
"""

def get_safs(wc, split="train"):
    """ Get SAF index paths from config based on wildcards """
    match = re.search(POP_PATTERN, wc.pops)
    pops = match.group(0).split(POP_SEP)
    safs = [config["input"][split][pop] for pop in pops]
    return safs


def get_only(lst):
    """ Get only element in list, or error """
    if len(lst) == 1:
        return lst[0]
    else:
        raise ValueError("list does not have exactly one element")


def get_realsfs_epochs(wc):
    """
    Get all epochs for realSFS.
    
    This really is a bit of a hack, since it relies on realSFS having already
    run to produce a log file. The realSFS epochs can vary due to the way 
    acceleration is implementing, giving jumps of two or three epochs at a time.
    The alternative is snakemake checkpoints, which is bothersome.
    The obvious implication is that realSFS must be run before any rule that
    requires inspecting which epochs were run.
    """
    log = get_only(expand(rules.realsfs.log, **wc))
    pattern = re.compile("(?<=Iteration )[0-9]+")
    lines = [pattern.search(line) for line in open(log)]
    epochs = [int(line[0]) for line in filter(None, lines)]
    return epochs


def get_winsfs_epochs():
    """ Get all epochs for window SFS """
    max_epochs = config["winsfs"]["max_epochs"]
    epochs = list(range(1, max_epochs + 1))
    return epochs
    

def get_winsfs_stream_epochs():
    """ Get all epochs for streaming window SFS """
    max_epochs = config["winsfs_stream"]["max_epochs"]
    epochs = list(range(1, max_epochs + 1))
    return epochs
    
    
def get_epochs(wc):
    """ Get all relevant epochs depending on methods """
    if str(wc.dir) == str(REALSFS_DIR):
        return get_realsfs_epochs(wc)
    elif str(wc.dir) == str(WINSFS_DIR):
        return get_winsfs_epochs()
    elif str(wc.dir) == str(WINSFS_STREAM_DIR):
        return get_winsfs_stream_epochs()


def get_conditions():
    """ Get names of different conditions to put in final CSV """
    stream = [""]
    if "winsfs_stream" in config:
        stream += ["stream_"]
    winsfs_conditions = expand(
        "{stream}winsfs_b{b}_w{w}", stream=stream, **config["run"]
    )
    # Note that the order here may be significant, 
    # if used in condition with get_gather_input()
    return ["realsfs"] + winsfs_conditions


def get_gather_input(wc):
    """
    Helper rule to gather everything for population(s) based on wildcards 
    """
    name = Path("{dir}") / "{pops}{maybe}.{csv}.csv"
    realsfs = expand(
        name,
        dir=REALSFS_DIR,
        pops=wc.pops,
        maybe="",
        csv=wc.csv,
    )
    winsfs_dirs = [WINSFS_DIR]
    if "winsfs_stream" in config:
        winsfs_dirs += [WINSFS_STREAM_DIR]
    winsfs = expand(
        name,
        dir=winsfs_dirs,
        pops=wc.pops,
        maybe=expand("_b{b}_w{w}", **config["run"]),
        csv=wc.csv,
    )
    # Note that the order here may be significant, 
    # if used in condition with get_conditions()
    lst = realsfs + winsfs
    return lst
