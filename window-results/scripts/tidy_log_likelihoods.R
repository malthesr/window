#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

read_log_likelihoods <- function(ps, epochs, split) {
    tibble(
        epoch = as.character(epochs),
        "{split}_loglik" := sapply(ps, scan, quiet = TRUE)
    )
}

epochs <- snakemake@params[["epochs"]]

train_paths <- snakemake@input[["train_loglik"]]
train_df <- read_log_likelihoods(train_paths, epochs, "train")

test_paths <- snakemake@input[["test_loglik"]]
if (!is.null(test_paths)) {
    test_df <- read_log_likelihoods(test_paths, epochs, "test")
    df <- full_join(train_df, test_df, by = "epoch")
} else {
    df <- train_df
}

csv <- snakemake@output[["csv"]]
write_csv(df, csv)
