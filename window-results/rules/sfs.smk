""" Rules for creating SFS in various ways """


include: "sfs/realsfs.smk"
include: "sfs/winsfs.smk"


if "winsfs_stream" in config:
    include: "sfs/stream_winsfs.smk"


include: "sfs/tidy.smk"


rule all_sfs:
    """ 
    Phony rule to produce all SFS based on config and gather into tidy CSVs
    """
    input:
        expand(
            rules.gather_sfs.output.csv,
            pops=config["run"]["pops"],
            csv=["sfs", "stat"],
        ),
