#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

path <- snakemake@input[["path"]]

df <- path |> 
    read_table(col_names = FALSE, col_types = cols(.default = col_integer())) |> 
    select(-102) |> # ANGSD leaves a trailing whitespace, leading to an empty col
    mutate(individual = row_number()) |>
    pivot_longer(-individual, names_to = "depth", values_to = "count") |>
    mutate(
        population =  gsub("_.*", "", basename(path)),
        depth = as.integer(gsub("^X", "", depth)) - 1,
    ) |>
    group_by(individual) |>
    mutate(density = count / sum(count)) |>
    ungroup() |>
    relocate(population)

csv <- snakemake@output[["csv"]]
write_csv(df, csv)