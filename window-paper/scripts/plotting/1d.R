read_sfs_csv <- function(p) {
    read_csv(p, col_types = "ciid")
}

prepare_sfs_df <- function(df) {
    df |>
        group_by(epoch) |>
        mutate(
            norm = normalise(x),
            var = !(i == 0 | i == max(i)),
            normvar = normalise(x * var)
        ) |>
        ungroup()
}

read_stat_csv <- function(p) {
    read_csv(p, col_types = "ciddd")
}

plot_sfs <- function(df, name, condition, colours, y_max) {
    var_df <- df |> filter(var)
    labels <- df |>
        group_by(epoch) |>
        summarise(fixed = sum(norm * !var)) |>
        mutate(label = glue("{epoch} ({percent(fixed, 0.01)})")) |>
        pull(label)

    ggplot(var_df, aes(i, normvar, fill = factor(epoch))) +
        geom_col(
            colour = "grey30",
            position = position_dodge(0.8),
            width = 0.8
        ) +
        annotate(
            geom = "richtext",
            x = 2,
            y = y_max,
            label = condition,
            size = 5,
            vjust = 1,
            hjust = 0,
            fill = NA,
            label.colour = NA,
        ) +
        scale_x_continuous(
            name = paste("Allele count,", name),
            breaks = seq(2, max(df$i), 2),
        ) +
        scale_y_continuous(
            name = "Frequency",
            limits = c(0, y_max),
            breaks = seq(0.1, 1, 0.1),
        ) +
        scale_fill_manual(
            name = "Epoch (fixed)",
            labels = labels,
            values = colours,
        ) +
        coord_cartesian(expand = FALSE) +
        theme(
            legend.background = element_rect(
                fill = "white",
                colour = "white",
            ),
            legend.position = c(0.98, 0.98),
            legend.justification = c(1, 1)
        )
}

plot_realsfs <- function(df, name, y_max) {
    colours <- brewer_pal(palette = "YlGnBu")(9)[c(3, 5, 7)]

    plot_sfs(
        df,
        name = name,
        condition = "realSFS",
        colours = colours,
        y_max = y_max
    )
}

plot_winsfs <- function(df, realsfs_df, name, y_max) {
    colours <- brewer_pal(palette = "YlOrRd")(9)[c(3, 5, 7)]

    last_realsfs_epoch <- realsfs_df |> pull(epoch) |> max()
    overlay_df <- realsfs_df |>
        filter(epoch == last_realsfs_epoch & var)

    plot_sfs(
        df,
        name = name,
        condition = "winsfs",
        colours = colours,
        y_max = y_max
    ) +
        geom_point(
            mapping = aes(i, normvar, size = 1),
            data = overlay_df,
            inherit.aes = FALSE,
            pch = 21,
            colour = "black",
            fill = NA,
        ) +
        scale_size_continuous(
            name = NULL,
            labels = paste("realSFS,\nepoch", last_realsfs_epoch),
            range = c(2)
        ) +
        guides(
            fill = guide_legend(order = 1),
            size = guide_legend(order = 2)
        ) +
        annotate(
            geom = "richtext",
            x = 2,
            y = y_max,
            label = "<span style='font-size:10pt; color:grey30'>W=100</span>",
            size = 6,
            vjust = 1.5,
            hjust = 0,
            fill = NA,
            label.colour = NA,
        ) +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.margin = margin(b = -10)
        )
}
