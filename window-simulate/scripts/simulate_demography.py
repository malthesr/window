#!/usr/bin/env python3

import demes
import msprime as ms


if __name__ == "__main__":
    par = snakemake.params

    mutation_rate = par["mutation_rate"]
    recombination_rate = par["recombination_rate"]
    chromosome_length = par["chromosome_length"]
    sample_sizes = par["sample_sizes"]
    seed = par["seed"]

    graph = demes.load(snakemake.input["yaml"])
    demography = ms.Demography.from_demes(graph)

    ts = ms.sim_ancestry(samples=sample_sizes,
                         demography=demography,
                         sequence_length=chromosome_length,
                         recombination_rate=recombination_rate,
                         discrete_genome=True,
                         ploidy=2,
                         random_seed=seed)
    mts = ms.sim_mutations(ts, rate=mutation_rate, random_seed=seed)
    mts.dump(snakemake.output["ts"])