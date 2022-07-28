#!/usr/bin/env python3

from collections import namedtuple

import demes
import demesdraw


def serial_bottleneck_model(
        ids,
        first_branch_time=10_000,
        branch_interval=1_000,
        normal_size=10_000,
        bottleneck_size=1_000, 
        bottleneck_duration=1_000,
        pulses=None,
):
    builder = demes.Builder(
        description="A simple serial bottlenecks demographic model.",
        time_units="generations",
    )
    populations = len(ids)
    builder.add_deme(ids[0], epochs=[dict(start_size=normal_size)])
    for i in range(1, populations):
        name = ids[i]
        ancestor = ids[i - 1]
        start_time = first_branch_time - branch_interval * (i - 1)
        bottleneck_epoch = dict(end_time=start_time - bottleneck_duration,
                                start_size=bottleneck_size,
                                end_size=normal_size)
        normal_epoch = dict(end_size=normal_size)
        builder.add_deme(name,
                         ancestors=[ancestor],
                         start_time=start_time,
                         epochs=[bottleneck_epoch, normal_epoch])
    if pulses is not None:
        for pulse in pulses:
            builder.add_pulse(sources=[pulse.source],
                              dest=pulse.dest,
                              time=pulse.time,
                              proportions=[pulse.proportion])
    return builder


Pulse = namedtuple("Pulse", "source dest time proportion")


if __name__ == "__main__":
    par = snakemake.params
    ids = par["ids"]

    if migrations := par["migrations"]:
        pulses = [Pulse(**d) for d in migrations]
    else:
        pulses = None

    builder = serial_bottleneck_model(ids=ids, pulses=pulses)
    graph = builder.resolve()
    
    out = snakemake.output
    demes.dump(graph, out["yaml"])

    ax = demesdraw.tubes(graph)
    ax.figure.savefig(out["model_fig"])

    ax = demesdraw.size_history(graph)
    ax.figure.savefig(out["sizes_fig"])
