# Analysing and tidying data

This subproject contains code for creating the main results for the paper and tidying these long CSVs files to be consumed by plotting. Input files are taken from the `window-input` and `window-simulate` subprojects via the YAML config files in `configs/` dir and are assumed to have been already created. See the corresponding subdirectories for how to do so.

## Running

The different data sets (impala, human 1000g, and simulations) run through the snakemake pipeline separately. A helper `run.fish` script exists to run all analyses and copy the CSV files required for plotting into the `data/` directory.

```shell
fish run.fish $threads
```

For the plotting code itself, see the `window-paper` sub-directory.