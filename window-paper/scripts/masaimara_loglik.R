#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/stats.R")
source("scripts/plotting/2d.R")

loglik_path <- "data/MasaiMara.loglik.csv"
loglik_df <- read_loglik_csv(loglik_path) |>
    filter(!startsWith(condition, "stream"))

# Train log-lik plot

train_loglik_plot <- plot_loglik(
    loglik_df,
    split = "train",
    first_realsfs_epoch = 7,
    main_y_breaks = c(0, -5e5, -1e6),
    inset_y_limits = c(-20, 2),
    inset_y_breaks = c(0, -10, -20)
)

# Test log-lik plot

test_loglik_plot <- plot_loglik(
    loglik_df,
    split = "test",
    first_realsfs_epoch = 7,
    main_y_breaks = c(0, -5e5, -1e6),
    inset_y_limits = c(-50, 5),
    inset_y_breaks = c(0, -25, -50)
)

# Compose

plot <- wrap_plots(
    train_loglik_plot,
    plot_spacer(),
    test_loglik_plot,
    nrow = 1,
    widths = c(1, 0.05, 1)
) +
    plot_layout(guides = "collect") +
    plot_annotation(
        title = "Impala",
        tag_levels = "a"
    ) &
    theme(
        legend.position = "bottom",
    )

ggsave(
    "figures/masaimara_loglik.png",
    plot,
    device = agg_png,
    width = 8,
    height = 4
)
