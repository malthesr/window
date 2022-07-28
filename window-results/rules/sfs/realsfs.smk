""" Rules for running realSFS """


rule install_realsfs:
    """
    Clones angsd, htslib, and replaces angsd realSFS with modified version;
    builds and installs to ./bin
    """
    input:
        mod="realsfs_mod/realSFS_optim.cpp",
    output:
        realsfs="bin/realSFS",
    params:
        angsd="angsd",
        angsd_git="https://github.com/ANGSD/angsd",
        angsd_rev="9bd93ee",
        htslib="htslib",
        htslib_git="https://github.com/samtools/htslib",
        htslib_rev="bd133ac",
        old_realsfs_optim="angsd/misc/realSFS_optim.cpp",
        old_realsfs="angsd/misc/realSFS",
    threads: 8
    shell:
        """
        (
            # Clone and build fixed revision of htslib
            git clone --single-branch --branch 'master' {params.htslib_git} {params.htslib};
            cd {params.htslib};
            git reset --hard {params.htslib_rev};
            git submodule update --init --recursive;
            make -j {threads};
            cd ..;
            # Clone and build fixed revision of angsd, using local htslib and modified realSFS
            git clone {params.angsd_git} {params.angsd};
            cd {params.angsd};
            git reset --hard {params.angsd_rev};
            cd ..;
            cp {input.mod} {params.old_realsfs_optim};
            cd {params.angsd};
            make -j {threads} HTSSRC=../{params.htslib};
            cd ..;
            # Move realSFS to output location
            mv {params.old_realsfs} {output.realsfs};
            # Clean-up, since temp(directory(..)) doesn't work
            rm -rf {params.angsd};
            rm -rf {params.htslib};
        ) 2> /dev/null > /dev/null
        """


rule realsfs:
    """ Run realSFS (with modified logging) """
    input:
        realsfs=rules.install_realsfs.output.realsfs,
        safs=get_safs,
    output:
        sfs=REALSFS_DIR / "{pops}.raw.sfs",
    params:
        tolerance=1e-10,
        max_iterations=config["realsfs"]["max_epochs"],
    log:
        REALSFS_DIR / "{pops}.log",
    threads: 20
    shell:
        """
        {input.realsfs} \
            -seed 1 \
            -v 1 \
            -P {threads} \
            -tole {params.tolerance} \
            -maxIter {params.max_iterations} \
            {input.safs} \
            2> {log} \
            > {output.sfs}
        """


rule format_realsfs:
    """ Add header containing dimensions to realSFS output """
    input:
        sfs=rules.realsfs.output.sfs,
        log=rules.realsfs.log,
    output:
        sfs=REALSFS_DIR / "{pops}.sfs",
    shell:
        """
        ( \
            grep -P '^\s\-\> dim\(.*\):[0-9]+' {input.log} \
                | grep -Po '[0-9]+$' \
                | paste -sd "/" - \
                | sed 's/^/#SHAPE</;s/$/>/' \
            && cat {input.sfs}
        ) > {output.sfs}
        """


rule nth_realsfs:
    """ Extract SFS after n full epochs from log """
    input:
        log=rules.realsfs.log,
        sfs=rules.format_realsfs.output.sfs,
    output:
        sfs=REALSFS_DIR / "iterations" / "{pops}.after{n}.sfs",
    params:
        pattern=lambda wc: f"'(?<=Iteration {wc.n}, current SFS: )[0-9\. ]+'",
    shell:
        """
        ( \
            head -n1 {input.sfs} && \
            grep -Pom1 {params.pattern} {input.log} \
        ) > {output.sfs}
        """


rule all_realsfs:
    """
    Phony rule to run realSFS
    
    This is especially handy since realSFS must run first to avoid checkpoint. 
    See also doc comment on get_realsfs_epochs in helpers.
    """
    input:
        expand(rules.format_realsfs.output, pops=config["run"]["pops"])
