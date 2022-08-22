#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/depth.R")

df <- read_csv("data/depth.csv", col_types = "cciiid")

plot <- df |>
    filter(population %in% c("MasaiMara", "Shangani") & split == "train") |>
    mutate(population = if_else(population == "MasaiMara", "Maasai Mara", population)) |>
    plot_depth()

ggsave(
    "figures/impala_depth.png",
    plot,
    device = agg_png,
    width = 8,
    height = 4
)

