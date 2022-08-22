#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/2d.R")
source("scripts/plotting/stats.R")
source("scripts/plotting/sim.R")

# Setup

stop_df <- read_sim_stop_df("data/sim.block_loglik.csv", peak = TRUE)

# SFS plots

sfs_df <- read_sim_sfs_csv("data/sim.sfs.csv", peak = TRUE)

truth_df <- sfs_df |>
    filter(condition == "truth") |>
    prepare_sfs_df() |>
    ungroup()

realsfs_df <- sfs_df |>
    filter(condition == "realsfs") |>
    group_by(d) |>
    filter(epoch == max(epoch)) |>
    prepare_sfs_df() |>
    ungroup() |>
    mutate(normvar = pmax(normvar, 1e-14)) |>
    left_join(
        truth_df,
        by = c("i", "j"),
        suffix = c("", "_truth")
    ) |>
    mutate(residual = normvar - normvar_truth) |>
    select(-ends_with("_truth"))

winsfs_df <- sfs_df |>
    filter(condition == "winsfs_b500_w100") |>
    right_join(
        filter(stop_df, condition == "winsfs_b500_w100"),
        by = c("n", "d", "condition", "epoch")
    ) |>
    select(-c(loglik, change)) |>
    group_by(d) |>
    prepare_sfs_df() |>
    ungroup() |>
    left_join(
        truth_df,
        by = c("i", "j"),
        suffix = c("", "_truth")
    ) |>
    mutate(residual = normvar - normvar_truth) |>
    select(-ends_with("_truth"))

get_values <- \(df) df |> filter(keep) |> pull(residual)
y_limits <- range(
    get_values(realsfs_df),
    get_values(winsfs_df),
    na.rm = TRUE
)

realsfs_plot_lst <- realsfs_df |>
    group_by(d) |>
    group_map(
        ~plot_realsfs2d_residual(
            .,
            names = SIM_NAMES,
            y_limits = y_limits,
            n.breaks = 10
        )
    )

winsfs_plot_lst <- winsfs_df |>
    group_by(d) |>
    group_map(
        ~plot_winsfs2d_residual(
            .,
            names = SIM_NAMES,
            y_limits = y_limits,
            n.breaks = 10
        )
    )

truth_plot <- plot_spacer() + plot_truesfs2d(
    truth_df,
    names = SIM_NAMES,
    y_limits = truth_df |> filter(keep) |> pull(normvar) |> range(),
    n.breaks = 6
) + plot_spacer()

ggsave(
    "figures/peak_truth.png",
    truth_plot,
    device = agg_png,
    width = 8,
    height = 3
)

winsfs_plot <- winsfs_plot_lst |>
    remove_axis_labels(nrow = 1) |>
    pseudo_facet_1x3_d() +
    plot_layout(guides = "collect")

ggsave(
    "figures/peak_winsfs.png",
    winsfs_plot,
    device = agg_png,
    width = 8,
    height = 3
)

realsfs_plot <- realsfs_plot_lst |>
    remove_axis_labels(nrow = 1) |>
    pseudo_facet_1x3_d() +
    plot_layout(guides = "collect")

ggsave(
    "figures/peak_realsfs.png",
    realsfs_plot,
    device = agg_png,
    width = 8,
    height = 3
)
