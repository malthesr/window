WINSFS_COLOURS <- viridis_pal(option = "plasma", begin = 0.3, end = 0.7)(3)
STAT_COLOURS <- c(
    "black",
    WINSFS_COLOURS
)

better_comma <- function(x, ...) {
    if_else(x == 0, "0", comma(x, ...))
}

plot_sim_loglik <- function(df,
                            stop_df,
                            main_y_breaks,
                            inset_y_limits,
                            inset_y_breaks) {
    loglik_df <- df |>
        rename("stat" = loglik) |>
        mutate(stat = stat - max(stat))

    y_min <- loglik_df |>
        filter(condition == "realsfs") |>
        slice_min(epoch) |>
        pull(stat)

    joined_stop_df <- stop_df |>
        left_join(loglik_df, by = c("n", "d", "condition", "epoch")) |>
        ungroup() |>
        select(epoch, stat, condition) |>
        bind_rows(
            loglik_df |>
                filter(condition == "realsfs") |>
                slice_max(epoch, n = 1) |>
                select(epoch, stat, condition)
        )

    stop_marker <- geom_point(
        aes(epoch, stat, fill = condition, shape = "Default stop"),
        joined_stop_df,
        size = 3,
        colour = "grey30",
        alpha = 0.75
    )
    stop_guide <- scale_shape_manual(
        name = NULL,
        values = 23,
        guide = guide_legend(override.aes = list(size = 3)),
    )

    plot_stat(
        loglik_df,
        main_x_scale = scale_x_continuous(
            name = "Epoch",
            trans = "log10",
            breaks = c(1, 10, 100)
        ),
        main_y_scale = scale_y_continuous(
            name = "Log-likelihood,<br>relative to highest",
            labels = \(x) better_comma(x, scale = 1/1000000, suffix = "M"),
            breaks = main_y_breaks,
        ),
        inset_x_scale = scale_x_continuous(
            name = "Epoch",
            trans = "log10",
            breaks = c(1, 5, 50)
        ),
        inset_y_scale = scale_y_continuous(
            labels = \(x) better_comma(x, scale = 1/1000, suffix = "K"),
            breaks = inset_y_breaks
        ),
        main_y_min = y_min,
        inset_y_min = inset_y_limits[1],
        inset_y_max = inset_y_limits[2],
        inset_x_max = 100,
        box_x_min = 8,
        box_x_max = 100,
        inset_extra = list(stop_marker, stop_guide, guides(shape = "none"))
    ) +
        stop_marker +
        stop_guide
}

plot_loglik <- function(df,
                        split,
                        first_realsfs_epoch,
                        main_y_breaks,
                        inset_y_limits,
                        inset_y_breaks) {

    loglik_df <- df |>
        rename("stat" = paste(split, "loglik", sep = "_")) |>
        mutate(stat = stat - max(stat))

    y_min <- loglik_df |>
        filter(condition == "realsfs" & epoch == first_realsfs_epoch) |>
        pull(stat)

    name <- paste(
        gsub("\\b([a-z])", "\\U\\1", split, perl = TRUE),
        "log-likelihood,<br><i>relative to highest</i>"
    )

    plot_stat(
        loglik_df,
        main_y_scale = scale_y_continuous(
            name = name,
            labels = \(x) better_comma(x, scale = 1/1000, suffix = "K"),
            breaks = main_y_breaks,
        ),
        inset_y_scale = scale_y_continuous(
            labels = comma,
            breaks = inset_y_breaks
        ),
        main_y_min = y_min,
        inset_y_min = inset_y_limits[1],
        inset_y_max = inset_y_limits[2]
    )
}

plot_f2 <- function(df, y_breaks) {
    last_realsfs_f2 <- df |>
        filter(condition == "realsfs") |>
        top_n(1, epoch) |>
        pull(f2)

    f2_df <- df |>
        rename("stat" = f2) |>
        mutate(stat = stat - last_realsfs_f2)

    y_max <- f2_df |>
        filter(condition != "realsfs") |>
        pull(stat) |>
        max()

    y_min <- f2_df |> pull(stat) |> min()

    plot_stat_main(
        f2_df,
        y_min = y_min,
        y_max = y_max
    ) +
        scale_x_continuous(
            name = "Epoch",
            trans = "log10",
            breaks = c(1, 5, 50, 500)
        ) +
        scale_y_continuous(
            name = "f<sub>2</sub>, <i>relative to last realSFS",
            labels = label_sci,
            breaks = y_breaks
        )
}

plot_fst <- function(df, first_realsfs_epoch, y_breaks) {
    last_realsfs_fst <- df |>
        filter(condition == "realsfs") |>
        top_n(1, epoch) |>
        pull(fst)

    fst_df <- df |>
        rename("stat" = fst) |>
        mutate(stat = stat - last_realsfs_fst)

    y_min <- fst_df |>
        filter(condition == "realsfs" & epoch == first_realsfs_epoch) |>
        pull(stat)

    y_max <- fst_df |> pull(stat) |> max()

    labels <- function(x) {
        x <- if_else(x == 0, "0",  comma(x))
        gsub("(?!^)[0]+$", "", x, perl = TRUE)
    }

    plot_stat_main(
        fst_df,
        y_min = y_min,
        y_max = y_max
    ) +
        scale_x_continuous(
            name = "Epoch",
            trans = "log10",
            breaks = c(1, 5, 50, 500)
        ) +
        scale_y_continuous(
            name = "Hudson's F<sub>st</sub>,<br><i>relative to last realSFS",
            labels = label_sci,
            breaks = y_breaks,
            sec.axis = sec_axis(
                \(x) x + last_realsfs_fst,
                name = "Hudson's F<sub>st</sub>",
                labels = label_comma(0.001),
                breaks = last_realsfs_fst + range(y_breaks)
            )
        ) +
        theme(
            axis.text.y.right = element_text(
                angle = 270,
                margin = margin(l = 0)
            )
        )
}

plot_theta <- function(df, first_realsfs_epoch, y_breaks) {
    last_realsfs_theta <- df |>
        filter(condition == "realsfs") |>
        top_n(1, epoch) |>
        pull(theta)

    theta_df <- df |>
        rename("stat" = theta) |>
        mutate(stat = stat - last_realsfs_theta)

    y_max <- theta_df |>
        filter(condition == "realsfs" & epoch == first_realsfs_epoch) |>
        pull(stat)

    plot_stat_main(
        theta_df,
        y_min = 0,
        y_max = y_max,
        sep = "\n"
    ) +
        scale_x_continuous(
            name = "Epoch",
            trans = "log10",
            breaks = c(1, 5, 50, 500)
        ) +
        scale_y_continuous(
            name = "Tajima's θ, <i>relative to last realSFS",
            labels = label_sci,
            breaks = y_breaks,
            sec.axis = sec_axis(
                \(x) x + last_realsfs_theta,
                name = "Tajima's θ",
                labels = function(x) label_sci(x, 3),
                breaks = last_realsfs_theta + range(y_breaks)
            )
        ) +
        guides(
            color = guide_legend(
                title.position = "top",
                title.hjust = 0.5,
                label.position = "bottom"
            )
        ) +
        theme(
            axis.text.y.right = element_text(
                angle = 270,
                margin = margin(l = 5)
            )
        )
}

plot_stat <- function(df,
                      main_x_scale = NULL,
                      main_y_scale,
                      inset_x_scale = NULL,
                      inset_y_scale,
                      main_y_min,
                      main_y_max = 0,
                      inset_y_min,
                      inset_y_max = 0,
                      inset_x_min = 1,
                      inset_x_max = 500,
                      box_x_min = 15,
                      box_x_max = 500,
                      box_y_min = main_y_min,
                      box_y_max = main_y_min / 5,
                      label_zoom = TRUE,
                      inset_extra = list()) {
    if (is.null(main_x_scale)) {
        main_x_scale <- scale_x_continuous(
            name = "Epoch",
            trans = "log10",
            breaks = c(1, 5, 50, 500)
        )
    }

    main <- plot_stat_main(
        df,
        y_min = main_y_min,
        y_max = main_y_max
    ) +
        main_x_scale +
        main_y_scale

    inset <- plot_stat_inset(
        df,
        scale_x = inset_x_scale,
        x_min = inset_x_min,
        x_max = inset_x_max,
        y_min = inset_y_min,
        y_max = inset_y_max
    ) +
        inset_y_scale +
        inset_extra

    p <- main +
        annotation_custom(
            grob = ggplotGrob(inset),
            xmin = log10(box_x_min),
            xmax = log10(box_x_max),
            ymin = box_y_min,
            ymax = box_y_max
        )

    if (label_zoom) {
        p <- p +
            annotate(
                "text",
                label = "Zoom",
                x = 10^((log10(box_x_min) + log10(box_x_max)) / 2),
                y = box_y_max,
                colour = "grey30",
                vjust = -0.1,
                hjust = 0.5,
                size = 3.5
            )
    }

    p
}

plot_stat_base <- function(df) {
    ggplot(df, aes(epoch, stat, fill = condition, colour = condition)) +
        geom_line(size = 1, alpha = 0.75) +
        geom_point(size = 1.5, alpha = 0.5, pch = 21) +
        scale_fill_manual(
            guide = "none",
            values = STAT_COLOURS,
        ) +
        theme(
            axis.ticks.x = element_line(
                colour = "grey30",
                size = 0.2
            ),
            axis.ticks.length.x = unit(2, "pt")
        )
}

plot_stat_main <- function(df, y_min, y_max = 0, sep = " ") {
    plot_stat_base(df) +
        scale_color_manual(
            name = NULL,
            values = STAT_COLOURS,
            labels = function(x) condition_labeller(x, sep = sep)
        ) +
        coord_cartesian(ylim = c(1.01 * y_min, y_max))
}

plot_stat_inset <- function(df, scale_x = NULL, x_min, x_max, y_min, y_max = 0) {
    if (is.null(scale_x)) {
        scale_x <- scale_x_continuous(
            name = "Epoch",
            trans = "log10",
            breaks = c(1, 10, 100)
        )
    }

    plot_stat_base(df) +
        coord_cartesian(
            xlim = c(x_min, x_max),
            ylim = c(y_min, y_max)
        ) +
        scale_x +
        scale_color_manual(
            guide = "none",
            values = STAT_COLOURS,
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.background = element_rect(
                fill = "white",
                colour = "grey30",
                size = 1
            ),
        )
}

