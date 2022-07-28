#!/usr/bin/env fish

# Call snakemake for the various outputs required.
# Takes a single argument setting the number of cores to use;
# if unset, will try to get the available number of cores (on linux).

if set -q argv[1]
    set jobs $argv[1]
else
    set -g jobs (grep processor /proc/cpuinfo | wc -l)
end
echo "Using $jobs cores"

alias snk 'python3 -m snakemake -p'

# Run simulations
set -l sim_config configs/sim.yaml
snk all_realsfs --configfile $sim_config -j $jobs
snk all --configfile $sim_config -j $jobs

# Run impala
set -l impala_config configs/impala.yaml
snk all_realsfs --configfile $impala_config -j $jobs
snk all --configfile $impala_config -j $jobs
snk all_bench --configfile $impala_config -j $jobs

# Run 1000g
set -l 1000g_config configs/1000g.yaml
snk all_realsfs --configfile $1000g_config -j $jobs
snk all --configfile $1000g_config -j $jobs

# Copy final CSVs to shared folder
set -l out data
mkdir -p $out
Rscript scripts/gather_sim.R results/sim/*.{{block_,}loglik,sfs,stat}.csv
cp results/{impala,1000g}/*.{{block_,}loglik,sfs,stat}.csv $out
cp results/impala/bench/*.bench.csv $out
