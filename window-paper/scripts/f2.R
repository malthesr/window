#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/2d.R")
source("scripts/plotting/stats.R")

# Human f2

h_stat_path <- "data/CEU-YRI.stat.csv"
h_stat_df <- h_stat_path |>
    read_stat_csv() |>
    filter(!startsWith(condition, "stream"))

h_plot <- plot_f2(
    h_stat_df,
    y_breaks = c(0, 5e-7, 1e-6)
) +
    ggtitle("Human")

# Impala f2

i_stat_path <- "data/Shangani-MasaiMara.stat.csv"
i_stat_df <- i_stat_path |>
    read_stat_csv() |>
    filter(!startsWith(condition, "stream"))

i_plot <- plot_f2(
    i_stat_df,
    y_breaks = c(0, 5e-6, 1e-5)
) +
    ggtitle("Impala")

# Compose

plot <- wrap_plots(
    h_plot,
    plot_spacer(),
    i_plot,
    nrow = 1,
    widths = c(1, 0.05, 1)
) +
    plot_annotation(
        tag_levels = "a",
    ) +
    plot_layout(guides = "collect") &
    theme(
        legend.position = "bottom",
    )

ggsave(
    "figures/f2.png",
    plot,
    device = agg_png,
    width = 8,
    height = 4
)
