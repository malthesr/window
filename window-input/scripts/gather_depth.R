#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

paths <- commandArgs(trailingOnly = TRUE)
names <- gsub("\\.depth\\.csv", "", basename(paths))

df <- paths |>
    lapply(read_csv, col_types = "ciiid") |>
    setNames(names) |>
    bind_rows(.id = "name") |>
    separate(name, into = c(NA, "split"), sep = "_") |>
    relocate(population)

write_csv(df, "data/depth.csv")
