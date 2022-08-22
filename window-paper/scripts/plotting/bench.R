BENCH_CONDITIONS <- c(
    "realsfs" = "realSFS",
    "winsfs" = "winsfs",
    "stream_winsfs" = "winsfs, stream",
    "stream_winsfs_with_shuffle" = "winsfs, stream\n(incl. shuffle)"
)

BENCH_COLOURS <- c(
    "realsfs" = STAT_COLOURS[1],
    "winsfs" = STAT_COLOURS[2],
    "stream_winsfs" = "steelblue",
    "stream_winsfs_with_shuffle" = "steelblue"
)

get_stop_epoch <- function(p,  tolerance = 1e-4) {
    read_csv(p, col_types = "cidd") |>
        filter(change <= tolerance) |>
        slice_head(n = 1) |>
        pull(epoch)
}

prepare_bench <- function(df) {
    run_df <- df |> filter(condition != "shuffle") |>
        mutate(fmt_condition = BENCH_CONDITIONS[condition])
    shuffle_df <- df |> filter(condition == "shuffle")
    stopifnot(nrow(shuffle_df) == 1)

    run_df |>
        filter(condition == "stream_winsfs") |>
        mutate(
            condition = "stream_winsfs_with_shuffle",
            seconds = seconds + shuffle_df$seconds,
            mb = max(mb, shuffle_df$mb),
            fmt_condition = BENCH_CONDITIONS[condition]
        ) |>
        bind_rows(run_df)
}

plot_time <- function(df, y_breaks) {
    ggplot(
        df,
        aes(
            epochs,
            seconds,
            colour = fmt_condition,
            shape = fmt_condition,
            group = fmt_condition
        )
    ) +
        geom_line() +
        geom_point(size = 2) +
        scale_x_continuous(
            name = "Epochs",
            breaks = c(1, seq(10, max(df$epochs), 10)),
        ) +
        scale_y_continuous(
            name = "Runtime",
            labels = label_comma(accurate = 1, scale = 1/60^2, suffix = "h"),
            breaks = y_breaks * 60^2
        ) +
        expand_limits(y = 0) +
        scale_colour_manual(
            name = NULL,
            values = unname(BENCH_COLOURS),
            guide = guide_legend(ncol = 2)
        ) +
        scale_shape_manual(
            name = NULL,
            values = c(16, 16, 16, 17),
            guide = guide_legend(ncol = 2)
        ) +
        scale_linetype(guide = "none") +
        coord_cartesian(clip = "off", expand = FALSE) +
        theme(
            plot.margin = margin(r = 5, t = 5),
            legend.position = c(1, 0),
            legend.justification = c(1, -0.25),
            legend.background = element_rect(
                fill = "white",
                colour = "white",
            ),
            legend.key.width = unit(1, "cm")
        )
}


plot_memory <- function(df) {
    conditions <- c("realsfs", "winsfs", "stream_winsfs_with_shuffle")
    first_df <- df |>
        filter(epochs == 1 & condition %in% conditions) |>
        mutate(
            bytes = mb * 2^20,
            fmt_condition = ordered(fmt_condition, BENCH_CONDITIONS)
        )

    ggplot(first_df, aes(fmt_condition, bytes, fill = condition)) +
        geom_col(
            colour = "grey30",
            alpha = 0.5
        ) +
        geom_text(
            aes(label = label_bytes()(bytes)),
            vjust = -0.2
        ) +
        scale_x_discrete(
            name = NULL
        ) +
        scale_y_continuous(
            name = "Memory usage",
            labels = label_bytes()
        ) +
        scale_fill_manual(
            values = BENCH_COLOURS,
            guide = "none"
        ) +
        coord_cartesian(clip = "off", expand = FALSE)
}
