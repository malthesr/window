#!/usr/bin/env Rscript

source("scripts/plotting/plotting.R")
source("scripts/plotting/stats.R")

prep <- function(p) {
    full_df <- read_csv(p, col_types = "cidd") |>
        filter(condition != "realsfs") |>
        extract(
            condition,
            into = c("method", "window_size"),
            regex = "([a-z_]+)_b[0-9]+_w([0-9]+)"
        )

    stream_df <- filter(full_df, method == "stream_winsfs")

    full_df |>
        filter(method == "winsfs") |>
        right_join(
            stream_df,
            by = c("window_size", "epoch"),
            suffix = c("_in_ram", "_stream")
        ) |>
        select(-starts_with("method")) |>
        pivot_longer(
            -c(epoch, window_size),
            names_to = c("split", ".value"),
            names_pattern = "(.*)_loglik_(.*)"
        ) |>
        mutate(diff = in_ram - stream)
}

capitalise <- \(x) gsub("([\\w])([\\w]+)", "\\U\\1\\L\\2", x, perl = TRUE)

pops <- c("YRI", "CEU-YRI", "MasaiMara", "Shangani-MasaiMara")

paths <- paste0("data/", pops, ".loglik.csv")
df <- paths |> lapply(prep) |>
    setNames(pops) |>
    bind_rows(.id = "population") |>
    arrange(desc(split)) |>
    mutate(
        fmt_population = fmt_populations(population, sep = "<br>"),
        fmt_split = capitalise(split)
    ) |>
    filter(window_size == 100)

plot <- ggplot(df, aes(epoch, diff, fill = window_size)) +
    geom_col(width = 0.9, position = position_dodge(0.9)) +
    facet_grid2(fmt_split ~ fmt_population, strip = strip_vanilla(clip = "off")) +
    scale_x_continuous(
        name = "Epoch",
        breaks = 1:max(df$epoch)
    ) +
    scale_y_continuous(
        name = "Log-likelihood difference,<br>in-memory minus streaming"
    ) +
    scale_fill_manual(
        guide = "none",
        name = NULL,
        values = WINSFS_COLOURS,
        labels = \(x) paste("W =", x)
    ) +
    coord_cartesian(clip = "off") +
    theme(
        legend.position = "bottom",
        plot.margin = margin(r = 25),
        legend.margin = margin(t = -5),
        strip.text.y = element_text(angle = 270, hjust = 0.5),
    )

ggsave(
    "figures/stream.png",
    plot,
    device = agg_png,
    width = 8,
    height = 5
)
