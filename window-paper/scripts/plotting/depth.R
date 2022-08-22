plot_depth <- function(df) {
    ggplot(df, aes(depth, density, group = individual)) +
        geom_line(alpha = 0.5) +
        scale_x_continuous(
            name = "Depth",
        ) +
        scale_y_continuous(
            name = "Density",
            labels = comma,
            breaks = seq(0.1, 1, 0.1)
        ) +
        coord_cartesian(xlim = c(0, 30), expand = FALSE) +
        facet_wrap(~population) +
        theme(
            panel.spacing = unit(20, "pt"),
            plot.margin = margin(r = 10)
        )
}
