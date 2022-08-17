#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)   
    library(readr)   
    library(tidyr)   
})

paths <- commandArgs(trailingOnly = TRUE) |> sort()

subset_paths <- \(p, pattern) grep(pattern, p, value = TRUE)
get_names <- \(p) gsub("\\..*\\.csv$", "", basename(p))

sfs_paths <- paths |> subset_paths("\\.sfs\\.csv$")
stat_paths <- paths |> subset_paths("\\.stat\\.csv$")
loglik_paths <- paths |> subset_paths("\\.loglik\\.csv$")
block_loglik_paths <- paths |> subset_paths("\\.block_loglik\\.csv$")

re <- "^A(_peak)?_n([0-9]+)_d([0-9]+)_e[0-9]+ppt-B(_peak)?_n[0-9]+_d[0-9]+_e[0-9]+.*"
read <- \(p, fmt) lapply(p, read_csv, col_types = fmt) |> 
    setNames(get_names(p)) |>
    bind_rows(.id = "name") |>
    extract(name, into = c("peak", "n", "d", NA), regex = re, convert = TRUE) |>
    mutate(peak = peak == "_peak")
dedup_truth <- \(df) df |> 
    filter(condition == "truth") |> 
    mutate(d = NA, epoch = NA) |> 
    distinct() |>
    arrange(n) |> 
    bind_rows(filter(df, condition != "truth"))
write <- \(df, p) write_csv(df, p)

sfs_df <- sfs_paths |> 
    read(fmt = "ciiid") |>
    dedup_truth() |>
    filter(condition %in% c("realsfs", "truth") | epoch < 20) # reduce csv size
write(sfs_df, "data/sim.sfs.csv")

stat_df <- stat_paths |> 
    read(fmt = "cidd") |> 
    dedup_truth()
write(stat_df, "data/sim.stat.csv")

loglik_df <- loglik_paths |> 
    read(fmt = "cid") |>
    rename("loglik" = train_loglik)
write(loglik_df, "data/sim.loglik.csv")

block_loglik_df <- block_loglik_paths |>
    read(fmt = "cidd")
write(block_loglik_df, "data/sim.block_loglik.csv")

