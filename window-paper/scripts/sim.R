#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/2d.R")
source("scripts/plotting/stats.R")
source("scripts/plotting/sim.R")

# Setup

stop_df <- read_sim_stop_df("data/sim.block_loglik.csv", peak = FALSE)

# SFS plots

sfs_df <- read_sim_sfs_csv("data/sim.sfs.csv", peak = FALSE)

truth_df <- sfs_df |>
    filter(condition == "truth") |>
    group_by(n) |>
    prepare_sfs_df() |>
    ungroup()

realsfs_df <- sfs_df |>
    filter(condition == "realsfs") |>
    group_by(n, d) |>
    filter(epoch == max(epoch)) |>
    prepare_sfs_df() |>
    ungroup()

winsfs_df <- sfs_df |>
    filter(condition == "winsfs_b500_w100") |>
    right_join(
        filter(stop_df, condition == "winsfs_b500_w100"),
        by = c("n", "d", "condition", "epoch")
    ) |>
    select(-c(loglik, change)) |>
    group_by(n, d) |>
    prepare_sfs_df() |>
    ungroup()

rm(sfs_df)

get_values <- \(df) df |> filter(keep) |> pull(normvar)
y_limits <- range(
    get_values(truth_df),
    get_values(realsfs_df),
    get_values(winsfs_df),
    na.rm = TRUE
)

truth_plot_lst <- truth_df |>
    group_by(n) |>
    group_map(
        ~plot_truesfs2d(
            .,
            names = SIM_NAMES,
            y_limits = y_limits,
            n.breaks = 10
        )
    ) |>
    remove_axis_labels(nrow = 1)

truth_plot <- pseudo_facet_1x3_n(truth_plot_lst) +
    plot_layout(guides = "collect") &
    theme(legend.key.height = unit(1, "cm"))

realsfs_plot_lst <- realsfs_df |>
    group_by(d, n) |>
    group_map(
        ~plot_realsfs2d(
            .,
            names = SIM_NAMES,
            y_limits = y_limits,
            n.breaks = 10
        )
    ) |>
    remove_axis_labels(nrow = 3)

realsfs_plot <- pseudo_facet_3x3(realsfs_plot_lst) +
    plot_layout(guides = "collect") &
    theme(legend.key.height = unit(2, "cm"))

winsfs_plot_lst <- winsfs_df |>
    group_by(d, n) |>
    group_map(
        ~plot_winsfs2d(
            .,
            names = SIM_NAMES,
            y_limits = y_limits,
            n.breaks = 10
        )
    ) |>
    remove_axis_labels(nrow = 3)

winsfs_plot <- pseudo_facet_3x3(winsfs_plot_lst) +
    plot_layout(guides = "collect") &
    theme(legend.key.height = unit(2, "cm"))

ggsave(
    "figures/sim_truth.png",
    truth_plot,
    device = agg_png,
    width = 8,
    height = 3
)

ggsave(
    "figures/sim_realsfs.png",
    realsfs_plot,
    device = agg_png,
    width = 8,
    height = 8
)

ggsave(
    "figures/sim_winsfs.png",
    winsfs_plot,
    device = agg_png,
    width = 8,
    height = 8
)

# Log-likelihood plots

loglik_df <- read_csv("data/sim.loglik.csv", col_types = "liicid") |>
    filter(peak == FALSE) |>
    select(-peak)

loglik_cfg <- list(
    "main_y_breaks" = list(
        "d2" = c(-4.5e8, -3e8, -1.5e8, 0),
        "d4" = c(  -9e7, -6e7,   -3e7, 0),
        "d8" = c(-7.5e6, -5e6, -2.5e6, 0)
    ),
    "inset_y_limits" = list(
        "n5"  = c( -25000,  100),
        "n10" = c(-250000, 1000),
        "n20" = c(-600000, 2500)
    ),
    "inset_y_breaks" = list(
        "n5"  = c(-30000,   -20000,  -10000, 0),
        "n10" = c(-300000, -200000, -100000, 0),
        "n20" = c(-750000, -500000, -250000, 0)
    )
)

loglik_p <- function(n, d) {
    plot_sim_loglik(
        df = loglik_df |> filter(n == !!n & d == !!d),
        stop_df = stop_df |> filter(n == !!n & d == !!d),
        main_y_breaks = loglik_cfg$main_y_breaks[[paste0("d", d)]],
        inset_y_limits = loglik_cfg$inset_y_limits[[paste0("n", n)]],
        inset_y_breaks = loglik_cfg$inset_y_breaks[[paste0("n", n)]]
    )
}

loglik_plot_lst <- expand_grid(d = c(2, 4, 8), n = c(5, 10, 20)) |>
    apply(1, \(x) loglik_p(x[2], x[1])) |>
    remove_axis_labels(n = 3)

loglik_plot <- pseudo_facet_3x3(loglik_plot_lst) +
    plot_layout(guides = "collect") &
    theme(
        legend.position = "bottom",
        plot.margin = margin(r = 5, t = 5),
        legend.margin = margin(t = -10),
    )

ggsave(
    "figures/sim_loglik.png",
    loglik_plot,
    device = agg_png,
    width = 8,
    height = 7
)
