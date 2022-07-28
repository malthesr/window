# Simulations

This subproject contains code for creating the simulations used in the paper.

## Running

All the required output can be created by running:

```shell
python3 -m snakemake all -j $threads
```

using whatever amount of threads desired.

By default, the output will be created in the `results/` directory in the placed expected by the config files for `window-results`.