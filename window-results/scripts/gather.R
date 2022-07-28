#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

paths <- snakemake@input
conditions <- snakemake@params[["conditions"]]
fmt <- snakemake@params[["fmt"]]

df <- paths |>
    lapply(read_csv, col_types = fmt) |>
    setNames(conditions) |>
    bind_rows(.id = "condition")

csv <- snakemake@output[["csv"]]
write_csv(df, csv)
