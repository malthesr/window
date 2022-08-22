#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/stats.R")
source("scripts/plotting/bench.R")

stop_epoch <- get_stop_epoch("data/Shangani-MasaiMara.block_loglik.csv")
df <- read_csv("data/Shangani-MasaiMara.bench.csv", col_types = "cicdd") |>
    prepare_bench() |>
    filter(condition == "realsfs" | epochs <= stop_epoch)

time_plot <- plot_time(df, y_breaks = c(0, 2, 4))
memory_plot <- plot_memory(df)

# Compose

plot <- wrap_plots(
    time_plot,
    plot_spacer(),
    memory_plot,
    widths = c(0.6, 0.03, 0.45)
) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(margin = margin(r = -5)))

ggsave(
    "figures/bench_impala.png",
    plot,
    device = agg_png,
    width = 8,
    height = 3
)
