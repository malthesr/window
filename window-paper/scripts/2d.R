#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/2d.R")
source("scripts/plotting/stats.R")

# Human SFS plot

h_sfs_df <- "data/CEU-YRI.sfs.csv" |> read_sfs_csv()
h_names <- c("CEU", "YRI")

h_realsfs_epoch <- 93
h_realsfs_df <- h_sfs_df |>
    filter(condition == "realsfs" & epoch == h_realsfs_epoch) |>
    prepare_sfs_df()

h_winsfs_epoch <- 1
h_winsfs_df <- h_sfs_df |>
    filter(condition == "winsfs_b500_w100" & epoch == h_winsfs_epoch) |>
    prepare_sfs_df()

h_y_limits <- range(
    h_realsfs_df |> filter(keep) |> pull(normvar),
    h_winsfs_df |> filter(keep) |> pull(normvar),
    na.rm = TRUE
)

h_realsfs_plot <- plot_realsfs2d(
    h_realsfs_df,
    names = h_names,
    y_limits = h_y_limits
) + ggtitle("Human")

h_winsfs_plot <- plot_winsfs2d(
    h_winsfs_df,
    names = h_names,
    y_limits = h_y_limits
) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
    )

h_sfs_plot <- h_realsfs_plot +
    h_winsfs_plot +
    plot_layout(guides = "collect")

# Human stat plots

h_loglik_df <- "data/CEU-YRI.loglik.csv" |>
    read_loglik_csv() |>
    filter(!startsWith(condition, "stream"))
h_stat_df <- "data/CEU-YRI.stat.csv" |>
    read_stat_csv() |>
    filter(!startsWith(condition, "stream"))

h_train_loglik_plot <- plot_loglik(
    h_loglik_df,
    split = "train",
    first_realsfs_epoch = 6,
    main_y_breaks = c(0, -2.5e4, -5e4),
    inset_y_limits = c(-200, 10),
    inset_y_breaks = c(-200, -100, 0)
)

h_test_loglik_plot <- plot_loglik(
    h_loglik_df,
    split = "test",
    first_realsfs_epoch = 6,
    main_y_breaks = c(0, -2.5e4, -5e4, -7.5e5),
    inset_y_limits = c(-1000, 50),
    inset_y_breaks = c(-800, -400, 0)
)

h_fst_plot <- plot_fst(
    h_stat_df,
    first_realsfs_epoch = 6,
    y_breaks = c(-1e-3, -5e-4, 0)
)

h_stat_plot <- wrap_plots(
    h_train_loglik_plot,
    h_test_loglik_plot,
    h_fst_plot,
    nrow = 1,
    widths = c(0.4, 0.4, 0.2)
) +
    plot_layout(guides = "collect") &
    theme(
        legend.position = "bottom",
        plot.margin = margin(r = 5),
        legend.margin = margin(t = -10),
    )

layout <- paste(
    c(rep("A", 99), "C", "\n", rep("B", 100)),
    collapse = ""
)

h_plot <- wrap_plots(
    A = h_sfs_plot,
    B = h_stat_plot,
    C = guide_area(),
    design = layout,
    heights = c(1, 0.8)
)

# Impala SFS plot

i_sfs_df <- "data/Shangani-MasaiMara.sfs.csv" |> read_sfs_csv()
i_names <- c("Shangani", "Maasai Mara")

i_realsfs_epoch <- 100
i_realsfs_df <- i_sfs_df |>
    filter(condition == "realsfs" & epoch == i_realsfs_epoch) |>
    prepare_sfs_df()

i_winsfs_epoch <- 2
i_winsfs_df <- i_sfs_df |>
    filter(condition == "winsfs_b500_w100" & epoch == i_winsfs_epoch) |>
    prepare_sfs_df()

i_y_limits <- range(
    i_realsfs_df |> filter(keep) |> pull(normvar),
    i_winsfs_df |> filter(keep) |> pull(normvar),
    na.rm = TRUE
)

i_realsfs_plot <- plot_realsfs2d(
    i_realsfs_df,
    names = i_names,
    y_limits = i_y_limits
) + ggtitle("Impala")

i_winsfs_plot <- plot_winsfs2d(
    i_winsfs_df,
    names = i_names,
    y_limits = i_y_limits
) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
    )

i_sfs_plot <- i_realsfs_plot +
    i_winsfs_plot +
    plot_layout(guides = "collect")

# Impala stat plots

i_loglik_df <- "data/Shangani-MasaiMara.loglik.csv" |>
    read_loglik_csv() |>
    filter(!startsWith(condition, "stream"))
i_stat_df <- "data/Shangani-MasaiMara.stat.csv" |>
    read_stat_csv() |>
    filter(!startsWith(condition, "stream"))

i_train_loglik_plot <- plot_loglik(
    i_loglik_df,
    split = "train",
    first_realsfs_epoch = 7,
    main_y_breaks = c(0, -7.5e5, -1.5e6),
    inset_y_limits = c(-200, 10),
    inset_y_breaks = c(0, -100, -200)
)

i_test_loglik_plot <- plot_loglik(
    i_loglik_df,
    split = "test",
    first_realsfs_epoch = 7,
    main_y_breaks = c(0, -7.5e5, -1.5e6),
    inset_y_limits = c(-300, 10),
    inset_y_breaks = c(0, -150, -300)
)

i_fst_plot <- plot_fst(
    i_stat_df,
    first_realsfs_epoch = 10,
    y_breaks = c(-1e-3, -5e-4, 0)
)

i_stat_plot <- wrap_plots(
    i_train_loglik_plot,
    i_test_loglik_plot,
    i_fst_plot,
    nrow = 1,
    widths = c(0.4, 0.4, 0.2)
) +
    plot_layout(guides = "collect") &
    theme(
        legend.position = "bottom",
        plot.margin = margin(r = 5),
        legend.margin = margin(t = -10),
    )

layout <- paste(
    c(rep("A", 99), "C", "\n", rep("B", 100)),
    collapse = ""
)

i_plot <- wrap_plots(
    A = i_sfs_plot,
    B = i_stat_plot,
    C = guide_area(),
    design = layout,
    heights = c(1, 0.8)
)

# Compose

plot <- h_plot / i_plot +
    plot_annotation(
        tag_levels = list(c(
            "a", "", "b", "c", "d",
            "e", "", "f", "g", "h"
        ))
    ) &
    theme(
        plot.tag = element_text(margin = margin(b = -7.5, r = 5)),
        plot.title = element_markdown(margin = margin(b = 5))
    )

ggsave(
    "figures/2d.png",
    plot,
    device = agg_png,
    width = 8,
    height = 9
)
