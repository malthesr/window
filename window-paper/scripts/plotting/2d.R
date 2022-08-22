read_sfs_csv <- function(p) {
    read_csv(p, col_types = "ciiid")
}

read_loglik_csv <- function(p) {
    read_csv(p, col_types = "cidd")
}

read_stat_csv <- function(p) {
    read_csv(p, col_types = "cidd")
}


guide_sfs <- function(...) {
    guide_coloursteps(
        frame.colour = "grey30",
        ticks.colour = "grey30",
        ticks = TRUE,
        direction = "vertical",
        title.position = "top",
        label.position = "right",
        ...
    )
}

scale_fill_sfs2d_residual <- function(palette = "RdBu",
                                  guide = guide_sfs(),
                                  ...) {
    scale_fill_distiller(
        palette = palette,
        guide = guide,
        ...
    )
}

scale_fill_sfs2d <- function(trans = "log10",
                             labels = trans_format("log10", math_format(10^.x)),
                             guide = guide_sfs(),
                             ...) {
    scale_fill_gradientn(
        colours = cet_pal(20, name = "r1"),
        trans = trans,
        labels = labels,
        guide = guide,
        ...
    )
}

prepare_sfs_df <- function(df) {
    df |>
        mutate(
            norm = normalise(x),
            var = !(i == 0 & j == 0 | i == max(i) & j == max(j)),
            keep = var | x <= max(x * var),
            normvar = normalise(x * keep),
        )
}

annotate_arrow <- function(label, x, y, xend, yend, curvature = 0.3, ...) {
    list(
        annotate(
            x = x,
            y = y,
            xend = xend,
            yend = yend,
            geom = "curve",
            colour = "black",
            curvature = curvature,
            arrow = arrow(length = unit(2, "mm")),
        ),
        annotate(
            x = x,
            y = y,
            label = label,
            colour = "black",
            geom = "richtext",
            fill = NA,
            label.colour = NA,
            size = 4,
            ...
        )
    )
}

plot_sfs2d <- function(df, condition, names, y_limits, ...) {
    df |>
        filter(keep) |>
        ggplot(aes(i, j, fill = normvar)) +
        geom_tile() +
        scale_x_continuous(
            name = paste("Allele count,", names[1]),
            breaks = seq(1, max(df$i), 4)
        ) +
        scale_y_continuous(
            name = paste("Allele count,", names[2]),
            breaks = seq(1, max(df$j), 4)
        ) +
        scale_fill_sfs2d(
            name = "Frequency",
            limits = y_limits,
            ...
        ) +
        annotate(
            geom = "richtext",
            x = 0.5,
            y = max(df$j),
            label = condition,
            colour = "black",
            size = 5,
            vjust = 1,
            hjust = 0,
            fill = NA,
            label.colour = NA,
        ) +
        coord_cartesian(expand = FALSE) +
        theme(
            legend.key.height = unit(0.8575, "cm"),
            legend.margin = margin(b = 16)
        )
}

plot_sfs2d_residual <- function(df, condition, names, y_limits, ...) {
    df |>
        filter(keep) |>
        ggplot(aes(i, j, fill = residual)) +
        geom_tile() +
        scale_x_continuous(
            name = paste("Allele count,", names[1]),
            breaks = seq(1, max(df$i), 4)
        ) +
        scale_y_continuous(
            name = paste("Allele count,", names[2]),
            breaks = seq(1, max(df$j), 4)
        ) +
        scale_fill_sfs2d_residual(
            name = "Frequency<br>residual",
            limits = y_limits,
        ) +
        annotate(
            geom = "richtext",
            x = 0.5,
            y = max(df$j),
            label = condition,
            colour = "black",
            size = 5,
            vjust = 1,
            hjust = 0,
            fill = NA,
            label.colour = NA,
        ) +
        coord_cartesian(expand = FALSE) +
        theme(
            legend.key.height = unit(0.8575, "cm"),
            legend.margin = margin(b = 16)
        )
}

annotate_fixed <- function(p, shape, fixed) {
    p <- p +
        annotate_arrow(
            percent(fixed[1], 0.01),
            x = 1,
            y = 3,
            xend = 0,
            yend = 0,
            hjust = 0
        )

    if (length(fixed) > 1) {
        p <- p + annotate_arrow(
            percent(fixed[2], 0.01),
            x = shape[1] - 1,
            y = shape[2] - 3,
            xend = shape[1],
            yend = shape[2],
            hjust = 1
        )
    }

    p
}

plot_truesfs2d <- function(df, names, y_limits, ...) {
    fixed <- df |>
        filter(j == 0 & i == 0 | j == max(j) & i == max(i)) |>
        filter(!keep) |>
        pull(norm)
    shape <- c(max(df$i), max(df$j))

    plot_sfs2d(
        df,
        condition = "True SFS",
        names = names,
        y_limits = y_limits,
        ...
    ) |>
        annotate_fixed(shape, fixed)
}

plot_realsfs2d <- function(df, names, y_limits, ...) {
    epoch <- df |> pull(epoch) |> unique()
    fixed <- df |>
        filter(j == 0 & i == 0 | j == max(j) & i == max(i)) |>
        filter(!keep) |>
        pull(norm)
    shape <- c(max(df$i), max(df$j))

    plot_sfs2d(
        df,
        condition = paste("realSFS, epoch", epoch),
        names = names,
        y_limits = y_limits,
        ...
    ) |>
        annotate_fixed(shape, fixed)
}

plot_realsfs2d_residual <- function(df, names, y_limits, ...) {
    epoch <- df |> pull(epoch) |> unique()
    shape <- c(max(df$i), max(df$j))
    plot_sfs2d_residual(
        df,
        condition = paste("realSFS, epoch", epoch),
        names = names,
        y_limits = y_limits,
        ...
    )
}

plot_winsfs2d <- function(df, names, y_limits, window_size = 100, ...) {
    epoch <- df |> pull(epoch) |> unique()
    fixed <- df |>
        filter(j == 0 & i == 0 | j == max(j) & i == max(i)) |>
        filter(!keep) |>
        pull(norm)
    shape <- c(max(df$i), max(df$j))

    plot_sfs2d(
        df,
        condition = paste("winsfs, epoch", epoch),
        names = names,
        y_limits = y_limits,
        ...
    ) |>
        annotate_fixed(shape, fixed) +
        annotate(
            geom = "richtext",
            x = 1,
            y = shape[2],
            label = paste0(
                "<span style='font-size:10pt; color:grey30'>",
                "W = ", window_size,
                "</span>"
            ),
            size = 6,
            vjust = 1.5,
            hjust = 0,
            fill = NA,
            label.colour = NA,
        )
}

plot_winsfs2d_residual <- function(df, names, y_limits, window_size = 100, ...) {
    epoch <- df |> pull(epoch) |> unique()
    shape <- c(max(df$i), max(df$j))
    plot_sfs2d_residual(
        df,
        condition = paste("winSFS, epoch", epoch),
        names = names,
        y_limits = y_limits,
        ...
    ) + annotate(
        geom = "richtext",
        x = 1,
        y = shape[2],
        label = paste0(
            "<span style='font-size:10pt; color:grey30'>",
            "W = ", window_size,
            "</span>"
        ),
        size = 6,
        vjust = 1.5,
        hjust = 0,
        fill = NA,
        label.colour = NA,
    )
}
