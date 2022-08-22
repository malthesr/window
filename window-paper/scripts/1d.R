#!/usr/bin/env R

source("scripts/plotting/plotting.R")
source("scripts/plotting/1d.R")
source("scripts/plotting/stats.R")

# YRI plots

yri_sfs_df <- "data/YRI.sfs.csv" |> read_sfs_csv()
yri_stat_df <- "data/YRI.stat.csv" |>
    read_stat_csv() |>
    filter(!startsWith(condition, "stream"))

yri_realsfs_df <- yri_sfs_df |>
    filter(condition == "realsfs" & epoch %in% c(6, 15, 39)) |>
    prepare_sfs_df()

yri_winsfs_df <- yri_sfs_df |>
    filter(condition == "winsfs_b500_w100" & epoch %in% c(1, 2, 3)) |>
    prepare_sfs_df()

yri_sfs_y_max <- max(yri_realsfs_df$normvar, yri_winsfs_df$normvar)

yri_realsfs_plot <- plot_realsfs(
    yri_realsfs_df,
    name = "YRI",
    y_max = yri_sfs_y_max
) + ggtitle("Human")

yri_winsfs_plot <- plot_winsfs(
    yri_winsfs_df,
    name = "YRI",
    realsfs_df = yri_realsfs_df,
    y_max = yri_sfs_y_max
)

yri_theta_plot <- plot_theta(
    yri_stat_df,
    first_realsfs_epoch = 6,
    y_breaks = c(0, 5e-6, 1e-5)
) + ggtitle("Human")

yri_sfs_plot <- yri_realsfs_plot + yri_winsfs_plot

# Maasai Mara plot

mm_sfs_df <- "data/MasaiMara.sfs.csv" |> read_sfs_csv()
mm_stat_df <- "data/MasaiMara.stat.csv" |>
    read_stat_csv() |>
    filter(!startsWith(condition, "stream"))

mm_realsfs_df <- mm_sfs_df |>
    filter(condition == "realsfs" & epoch %in% c(10, 40, 100)) |>
    prepare_sfs_df()

mm_winsfs_df <- mm_sfs_df |>
    filter(condition == "winsfs_b500_w100" & epoch %in% c(1, 2, 3)) |>
    prepare_sfs_df()

mm_sfs_y_max <- max(mm_realsfs_df$normvar, mm_winsfs_df$normvar)

mm_realsfs_plot <- plot_realsfs(
    mm_realsfs_df,
    name = "Maasai Mara",
    y_max = mm_sfs_y_max
) + ggtitle("Impala")

mm_winsfs_plot <- plot_winsfs(
    mm_winsfs_df,
    name = "Maasai Mara",
    realsfs_df = mm_realsfs_df,
    y_max = mm_sfs_y_max
)

mm_sfs_plot <- mm_realsfs_plot + mm_winsfs_plot

mm_theta_plot <- plot_theta(
    mm_stat_df,
    first_realsfs_epoch = 7,
    y_breaks = c(0, 2e-4, 4e-4)
) + ggtitle("Impala")

# Compose

sfs_plot <- wrap_plots(
    yri_sfs_plot,
    mm_sfs_plot,
    ncol = 1
)

theta_plot <- wrap_plots(
    yri_theta_plot,
    mm_theta_plot,
    ncol = 1
) +
    plot_layout(guides = "collect")

plot <- wrap_plots(
    sfs_plot,
    plot_spacer(),
    theta_plot,
    widths = c(0.8, 0.01, 0.2)
) +
    plot_annotation(
        tag_levels = list(c("a", "", "c", "", "b", "d"))
    ) &
    theme(
        plot.tag = element_text(margin = margin(b = -7.5)),
        plot.title = element_markdown(margin = margin(b = 5))
    )

ggsave(
    "figures/1d.png",
    plot,
    device = agg_png,
    width = 8,
    height = 6
)


