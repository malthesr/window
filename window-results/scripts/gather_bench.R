#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

paths <- unlist(snakemake@input)
run_paths <- grep("sfs\\.benchmark$", paths, value = TRUE)
shuffle_paths <- grep("shuffle\\.benchmark$", paths, value = TRUE)

name <- \(p) gsub("\\.benchmark$", "", basename(p))
read <- \(p) read_tsv(p, col_types = "d_d________")
prep <- \(p) p |> lapply(read) |> setNames(name(p)) |> bind_rows(.id = "name")

run_df <- run_paths |>
    prep() |>
    extract(
        name,
        into = c("pops", "epochs", "condition"),
        regex = "^([A-Za-z-]+)_i([0-9]+).([A-Za-z_-]+)$",
        convert = TRUE
    ) 
    
shuffle_df <- shuffle_paths |>
    prep() |>
    separate(name, into = c("pops", NA), sep = "\\.") |>
    mutate(epochs = NA_integer_, condition = "shuffle")
    
df <- bind_rows(run_df, shuffle_df) |>
    rename(
        "seconds" = s,
        "mb" = max_rss,
    )

csv <- snakemake@output[["csv"]]
write_csv(df, csv)
