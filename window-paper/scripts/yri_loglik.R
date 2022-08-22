#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/stats.R")
source("scripts/plotting/2d.R")

loglik_path <- "data/YRI.loglik.csv"
loglik_df <- read_loglik_csv(loglik_path) |>
    filter(!startsWith(condition, "stream"))

# Train log-lik plot

train_loglik_plot <- plot_loglik(
    loglik_df,
    split = "train",
    first_realsfs_epoch = 6,
    main_y_breaks = c(0, -2.5e4, -5e4),
    inset_y_limits = c(-50, 2),
    inset_y_breaks = c(0, -25, -50)
)

# Test log-lik plot

test_loglik_plot <- plot_loglik(
    loglik_df,
    split = "test",
    first_realsfs_epoch = 6,
    main_y_breaks = c(0, -2.5e4, -5e4),
    inset_y_limits = c(-250, 10),
    inset_y_breaks = c(0, -100, -200)
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
        title = "Human",
        tag_levels = "a"
    ) &
    theme(
        legend.position = "bottom",
    )

ggsave(
    "figures/yri_loglik.png",
    plot,
    device = agg_png,
    width = 8,
    height = 4
)
