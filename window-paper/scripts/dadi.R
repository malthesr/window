#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/dadi.R")

# Setup

residual_path <- "data/dadi_residuals.csv"
residual_df <- read_csv(residual_path, col_types = "cciid") |>
    filter(!(i == 0 & j == 0))

names <- c("Shangani", "Maasai Mara")
limits <- range(residual_df$residual)
breaks <- pretty(limits) - 25
breaks <- c(breaks, breaks[length(breaks)] + 50)

remove_x_axis <- theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
)
remove_y_axis <-
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
    )
remove_axes <- list(remove_x_axis, remove_y_axis)

# First model residuals

first_df <- residual_df |> filter(model == "first")

first_realsfs_plot <- plot_residuals(
    first_df |> filter(method == "realsfs"),
    condition = "realSFS,<br>epoch 100",
    names = names,
    limits = limits,
    breaks = breaks
) +
    remove_x_axis

first_winsfs_plot <- plot_residuals(
    first_df |> filter(method == "winsfs"),
    condition = "winsfs,<br>epoch 2",
    names = names,
    limits = limits,
    breaks = breaks
) +
    remove_axes

# Second model residuals

second_df <- residual_df |> filter(model == "second")

second_realsfs_plot <- plot_residuals(
    second_df |> filter(method == "realsfs"),
    condition = "realSFS,<br>epoch 100",
    names = names,
    limits = limits,
    breaks = breaks
) +
    remove_x_axis

second_winsfs_plot <- plot_residuals(
    second_df |> filter(method == "winsfs"),
    condition = "winsfs,<br>epoch 2",
    names = names,
    limits = limits,
    breaks = breaks
) +
    remove_axes

# Third model residuals

third_df <- residual_df |> filter(model == "third")

third_realsfs_plot <- plot_residuals(
    third_df |> filter(method == "realsfs"),
    condition = "realSFS,<br>epoch 100",
    names = names,
    limits = limits,
    breaks = breaks
)

third_winsfs_plot <- plot_residuals(
    third_df |> filter(method == "winsfs"),
    condition = "winsfs,<br>epoch 2",
    names = names,
    limits = limits,
    breaks = breaks
) +
    remove_y_axis

# Models

model_colours <- c("#CD8D5A", "#72B173")

hack_legend_df <- tibble(
    x = -Inf,
    y = -Inf,
    fill = rep(c("realSFS", "winsfs"), 6)
)

base_model <- plot_basic_model() +
    geom_tile(
        aes(x, y, fill = fill),
        hack_legend_df,
        inherit.aes = FALSE,
        colour = "grey30"
    ) +
    scale_fill_manual(
        name = NULL,
        values = model_colours,
        guide = guide_legend(
            label.position = "bottom"
        )
    )

first_par <- read_par(
    "data/dadi_first_parameters.csv",
    col_types = "cdddddd",
    col_sci = "m",
    colours = model_colours
)

first_model <- base_model |>
    add_shared_size(label = first_par$n) |>
    add_split_sizes(labels = c(first_par$n1, first_par$n2)) |>
    add_constant_time(label = first_par$t) |>
    add_symmetric_migration(label = first_par$m) |>
    add_log_likelihood(label = first_par$log_likelihood) +
    remove_axes

second_par <- read_par(
    "data/dadi_second_parameters.csv",
    col_types = "cddddddd",
    col_sci = c("m12", "m21"),
    colours = model_colours
)

second_model <- base_model |>
    add_shared_size(label = second_par$n) |>
    add_split_sizes(labels = c(second_par$n1, second_par$n2)) |>
    add_constant_time(label = second_par$t) |>
    add_asymmetric_migration(labels = c(second_par$m1, second_par$m2)) |>
    add_log_likelihood(label = second_par$log_likelihood) +
    remove_axes

third_par <- read_par(
    "data/dadi_third_parameters.csv",
    col_types = "cdddddddddd",
    col_sci = c("m12", "m21"),
    colours = model_colours
)

third_model <- base_model |>
    add_shared_size(label = third_par$n) |>
    add_split_sizes(labels = c(third_par$n1a, third_par$n2a), y = 0.8) |>
    add_split_sizes(labels = c(third_par$n1b, third_par$n2b), y = 0.2) |>
    add_split_time(labels = c(third_par$t1, third_par$t2)) |>
    add_asymmetric_migration(labels = c(third_par$m1, third_par$m2), y = 0.55) |>
    add_log_likelihood(label = third_par$log_likelihood) +
    scale_x_continuous(
        name = "Shangani⠀⠀⠀⠀Maasai Mara⠀⠀",
        breaks = NULL,
    ) +
    remove_y_axis

# Compose

residual_plot <- wrap_plots(
    first_realsfs_plot, first_winsfs_plot,
    second_realsfs_plot, second_winsfs_plot,
    third_realsfs_plot, third_winsfs_plot,
    nrow = 3
) +
    plot_layout(
        guides = "collect"
    )

model_plot <- wrap_plots(
    first_model,
    second_model,
    third_model,
    nrow = 3
) +
    plot_layout(
        guides = "collect"
    )

plot <- wrap_plots(
    model_plot,
    plot_spacer(),
    residual_plot,
    nrow = 1,
    widths = c(0.3, 0.05, 0.6)
) +
    plot_annotation(
        tag_levels = list(c("a", "c", "e", "b", "", "d", "", "f"))
    ) &
    theme(
        plot.tag = element_text(margin = margin(b = -7.5, t = 10)),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        legend.title = element_markdown(vjust = 0.8),
        legend.margin = margin(t = 8, r = 25)
    )

ggsave(
    "figures/dadi.png",
    plot,
    device = agg_png,
    width = 8,
    height = 8
)

