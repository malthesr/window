""" Rules for running winsfs in streaming mode """


rule shuffle:
    """
    Perform block pseudo-shuffle of SAFs on disk
    in preparation for running winsfs out of memory
    """
    input:
        winsfs=rules.install_winsfs.output.winsfs,
        safs=get_safs,
    output:
        saf_shuf=WINSFS_STREAM_DIR / "{pops}.saf.shuf",
    params:
        blocks=config["winsfs_stream"]["shuffle_blocks"],
    benchmark:
        WINSFS_STREAM_DIR / "{pops}.shuffle.benchmark"
    log:
        WINSFS_STREAM_DIR / "{pops}.log",
    threads: 1
    shell:
        """
        {input.winsfs} shuffle \
            -vvv \
            --blocks {params.blocks} \
            --output {output.saf_shuf} \
            {input.safs} \
            2> {log}
        """


rule stream_winsfs:
    """ Run winsfs out of memory """
    input:
        winsfs=rules.install_winsfs.output.winsfs,
        saf_shuf=rules.shuffle.output.saf_shuf,
    output:
        sfs=WINSFS_STREAM_DIR / "{pops}_b{b}_w{w}.sfs",
    params:
        max_epochs=config["winsfs_stream"]["max_epochs"],
    log:
        WINSFS_STREAM_DIR / "{pops}_b{b}_w{w}.log",
    threads: 1
    shell:
        """
        {input.winsfs} \
            -vvv \
            --seed 1 \
            --max-epochs {params.max_epochs} \
            --blocks {wildcards.b} \
            --window-size {wildcards.w} \
            {input.saf_shuf} \
            > {output.sfs} \
            2> {log}
        """


rule nth_stream_winsfs:
    """
    Extract SFS after n full epochs from log.

    The shape information occurs at the top of the shuffle log.
    """
    input:
        log=rules.stream_winsfs.log,
        shuffle_log=rules.shuffle.log,
    output:
        sfs=WINSFS_STREAM_DIR / "iterations" / "{pops}_b{b}_w{w}.after{n}.sfs",
    params:
        header_pattern="^INFO  \[init\] Pre-allocating [0-9]+ bytes on disk for [0-9]+ for populations with shape [0-9/]+$",
        sfs_pattern=rules.nth_winsfs.params.sfs_pattern,
    shell:
        """
        (
            grep -Pm 1 '{params.header_pattern}' {input.shuffle_log} | 
                grep -Po '[0-9/]+$' | 
                sed 's/^/#SHAPE=</;s/$/>/' &&
            grep -P '{params.sfs_pattern}' {input.log} |
                sed '{wildcards.n}q;d' |
                grep -Po '([0-9\.]+( )?)+$'
        ) > {output.sfs}
        """
