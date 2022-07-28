#!/usr/bin/env Rscript

source("scripts/sfs.R")

sfs_path <- snakemake@input[["sfs"]]
true_sfs_path <- snakemake@params[["true_sfs"]]

sfs <- read_sfs(sfs_path)
true_sfs <- read_sfs(true_sfs_path)

stopifnot(dim(sfs) == dim(true_sfs))
stopifnot(abs(sum(sfs) - sum(true_sfs)) < 1)

# Add epsilon where 0 (if any) for log below
sfs[sfs < 1e-4] <- 1e-4

v <- sum(true_sfs * (sfs |> normalise() |> log()))

cat(v, file = snakemake@output[["log_likelihood"]], sep = "\n")