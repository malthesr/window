#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/stats.R")

# Prepare

loglik_paths <- list.files(
    "data",
    "(MasaiMara|YRI).loglik.csv",
    full.names = TRUE
)
names <- gsub("\\.loglik.csv$", "", basename(loglik_paths))

loglik_df <- loglik_paths |>
    lapply(
        read_csv,
        col_types = "cid",
        col_select = c("condition", "epoch", "test_loglik")
    ) |>
    setNames(names) |>
    bind_rows(.id = "population") |>
    filter(condition != "realsfs") |>
    mutate(
        window_size = gsub(".*w(?=[0-9]{3})", "", condition, perl = TRUE),
        window_size = as.integer(window_size),
        population = fmt_populations(population)
    ) |>
    group_by(population, window_size) |>
    top_n(1, test_loglik) |>
    select(population, window_size, epoch, test_loglik)


# Clean

block_loglik_paths <- list.files(
    "data",
    "(MasaiMara|YRI)\\.block_loglik\\.csv",
    full.names = TRUE
)

df <- block_loglik_paths |>
    lapply(read_csv, col_types = "cidd") |>
    setNames(names) |>
    bind_rows(.id = "population") |>
    extract(
        condition,
        into = "window_size",
        regex = "winsfs_b[0-9]+_w([0-9]+)",
        convert = TRUE,
    ) |>
    mutate(population = fmt_populations(population))

cutoff <- 1e-4
cutoff_df <- df |>
    filter(change < cutoff) |>
    group_by(population, window_size) |>
    top_n(1, change)
max_epoch <- 10

plot <- ggplot(df, aes(epoch, change, colour = factor(window_size), group = window_size)) +
    geom_line(key_glyph = "path") +
    geom_point(data = cutoff_df) +
    geom_hline(
        aes(yintercept = cutoff),
        data = tibble(cutoff),
        linetype = "dashed",
        colour = "grey30",
        key_glyph = "blank"
    ) +
    geom_vline(
        aes(xintercept = epoch, colour = factor(window_size)),
        loglik_df,
        linetype = "dashed",
        size = 0.75,
        key_glyph = "blank"
    ) +
    geom_label_repel(
        aes(label = epoch),
        data = cutoff_df,
        fill = alpha("white", alpha = 0.8),
        key_glyph = "blank"
    ) +
    facet_wrap(~population) +
    scale_x_continuous(
        name = "Epoch",
        limits = c(-0.5, max(df$epoch)),
    ) +
    scale_y_continuous(
        name = "Block log-likelihood sum difference",
        breaks = c(cutoff * 10, cutoff * 5, 0),
        labels = label_sci
    ) +
    scale_color_manual(
        name = NULL,
        values = viridis_pal(option = "plasma", begin = 0.3, end = 0.7)(3),
        labels = function(x) paste("W =", x)
    ) +
    coord_cartesian(
        expand = FALSE,
        ylim = c(0, cutoff * 11)
    ) +
    theme(
        panel.spacing = unit(20, "pt"),
        legend.position = "bottom"
    )

ggsave(
    "figures/stop.png",
    plot,
    device = agg_png,
    width = 8,
    height = 6
)
