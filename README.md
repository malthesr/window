# Code and data for winsfs paper

This repository contains the code used for the results shown in the paper "Estimation of site frequency spectra from low-coverage sequencing data using stochastic EM reduces overfitting, runtime, and memory usage". The manuscript is available on [bioRxiv](https://doi.org/10.1101/2022.05.24.493190). The winsfs software itself can be found in a separate repository [here](https://github.com/malthesr/winsfs). This repository represents an attempt to make the results more reproducible, as well as to make available the final, tidied data for anyone who wants to play around with it.

## Overview

The repository contains a number of sub-projects. Each sub-folder should contain a separate `README` to give a description of its contents, but the following should give an overview:

- `window-input`: A snakemake pipeline for processing the real-world data used in the article (human 1000G, impala) from BAM to SAF files split into train/test sets.
- `window-simulate`: A snakemake pipeline for creating the simulated data used in the article and processing the true spectra as well as producing SAF files.
- `window-results`: A snakemake pipeline for almost all processing from SAF files onwards, including SFS creation, data tidying, benchmarking, and more.
- `window-paper`: Tidied final results data in CSV format, plotting code, and LaTeX code for the article itself.

## Accessing the final results data

If you wish to play around with the final data used for the paper, this can be found in tidy long-format CSVs in the `window-paper/data` sub-folder.

## Reproducing the results

While this repository is an attempt at increased reproducibility, this repository unfortunately does not represent a fully self-contained, one-click reproducible workflow. In particular, for the real-world data, it assumes the raw data has already been downloaded and processed to BAMs. More information is available in the `window-input` sub-directory. Moreover, dependency management is sadly lacking (see below for a few pointers to remedy this). Apart from these caveats, however, it should be possible to fairly closely reproduce the results shown in the paper. To do so, start in the `window-input` and `window-simulate` sub-folders, which are independent of each other, to produce SAF files, and then move on to `window-results`. See the respective READMEs for detail.

Feel free open an issue if you're trying to reproduce the results and run into problems.

### Dependencies

Running the main workflows in the various sub-directories will require at least the following:

- `snakemake`
- `R` (>4.1) with `dplyr`, `tidyr`, and `readr` installed
- `python` (>3.9) with `numpy`, `msprime`, `tskit`, `demes`, and `demesdraw` installed
- A Rust toolchain, e.g. `cargo`
- `bcftools`/`samtools`
- `angsd`
- Optionally, a `fish` shell to run some of the helper scripts
- Common utilities like `sed`, `awk`, `perl`, `cut`, etc. are also assumed to be present.

The plotting scripts in the `window-paper` folder have significantly more R dependencies, but these can be easily grepped in the `scripts/` folder and installed if required.
