""" Rules for running winsfs """


rule install_winsfs:
    """
    Installs winsfs to ./bin
    """
    output:
        winsfs="bin/winsfs",
    params:
        winsfs_git="https://github.com/malthesr/winsfs",
        rev="13c4dcd",
        root=workflow.basedir,
    threads: 8
    shell:
        """
        cargo install -qj {threads} --git {params.winsfs_git} --rev {params.rev} --root {params.root}
        """


rule winsfs:
    """ Run winsfs in memory """
    input:
        winsfs=rules.install_winsfs.output.winsfs,
        safs=get_safs,
    output:
        sfs=WINSFS_DIR / "{pops}_b{b}_w{w}.sfs",
    params:
        max_epochs=config["winsfs"]["max_epochs"],
    log:
        WINSFS_DIR / "{pops}_b{b}_w{w}.log",
    threads: 20
    shell:
        """
        {input.winsfs} \
            -vvv \
            --seed 1 \
            --threads {threads} \
            --max-epochs {params.max_epochs} \
            --blocks {wildcards.b} \
            --window-size {wildcards.w} \
            {input.safs} \
            > {output.sfs} \
            2> {log}
        """


rule nth_winsfs:
    """
    Extract SFS after n full epochs from log.

    The shape information occurs once at the top and has to be combined with the
    appropriate epoch SFS to feed into e.g. log-likelihood calculations.
    """
    input:
        log=rules.winsfs.log,
    output:
        sfs=WINSFS_DIR / "iterations" / "{pops}_b{b}_w{w}.after{n}.sfs",
    params:
        header_pattern="^DEBUG \[init\] Found [0-9]+ \(intersecting\) sites in SAF files with shape [0-9/]+$",
        sfs_pattern="^DEBUG \[windowem\] Current SFS: [0-9\. ]+$",
    shell:
        # The log files are big since winsfs is run with '-vvv', so we jump
        # some regex hoops to avoid look-aheads and such.
        """
        ( \
            grep -Pm 1 '{params.header_pattern}' {input.log} | \
                grep -Po '[0-9/]+$' | \
                sed 's/^/#SHAPE=</;s/$/>/' && \
            grep -P '{params.sfs_pattern}' {input.log} | \
                sed '{wildcards.n}q;d' | \
                grep -Po '([0-9\.]+( )?)+$' \
        ) > {output.sfs}
        """
