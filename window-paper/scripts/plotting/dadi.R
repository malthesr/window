plot_residuals <- function(df, condition, names, limits, breaks) {
    ggplot(df, aes(i, j, fill = residual)) +
        geom_tile() +
        scale_x_continuous(
            name = paste("Allele count,", names[1]),
            breaks = seq(1, max(df$i), 4)
        ) +
        scale_y_continuous(
            name = paste("Allele count,", names[2]),
            breaks = seq(1, max(df$j), 4)
        ) +
        scale_fill_continuous_diverging(
            name = "Residual",
            palette = "Blue-Red 3",
            limits = limits,
            breaks = breaks,
            guide = guide_coloursteps(
                frame.colour = "grey30",
                ticks.colour = "grey30",
                ticks = TRUE
            )
        ) +
        annotate(
            geom = "richtext",
            x = Inf,
            y = max(df$j),
            label = condition,
            colour = "black",
            size = 5,
            vjust = 1,
            hjust = 1,
            fill = NA,
            label.colour = NA,
        ) +
        coord_cartesian(expand = FALSE) +
        theme(
            panel.grid.major.y = element_blank()
        )
}

read_par <- function(path, col_types, col_sci, col_k, colours) {
    par <- read_csv(path, col_types = col_types) |>
        arrange(method) |>
        mutate(
            across(all_of(col_sci), label_sci, digits = 2),
            across(where(is.numeric) & !all_of(col_sci), label_k)
        ) |>
        lapply(function(col) label(col[1], col[2], colours = colours))

    sapply(
        names(par),
        function(name) {
            ifelse(
                name == "log_likelihood",
                gsub("<br>", " ", par[name]),
                par[name]
            )

        }
    )
}

plot_basic_model <- function(names = NULL, ...) {
    df <- tribble(
        ~x,  ~y, ~xend, ~yend,
        0.0, 0.0,   0.0,   1.1,
        0.0, 1.1,   0.2,   1.1,
        0.2, 1.1,   0.2,   1.4,
        0.2, 1.4,   0.4,   1.4,
        0.4, 1.4,   0.4,   1.1,
        0.4, 1.1,   0.6,   1.1,
        0.6, 1.1,   0.6,   0.0,
        0.6, 0.0,   0.4,   0.0,
        0.4, 0.0,   0.4,   1.0,
        0.4, 1.0,   0.2,   1.0,
        0.2, 1.0,   0.2,   0.0,
        0.2, 0.0,   0.0,   0.0,
    )

    p <- ggplot(df, aes(x, y, xend = xend, yend = yend)) +
        geom_segment() +
        geom_polygon(colour = NA, fill = "grey95") +
        coord_cartesian(clip = "off", expand = FALSE) +
        theme(
            panel.grid.major.y = element_blank()
        )

    if (!is.null(names)) {
        p <- p +
            text(x = 0.1, y = 0, label = names[1], vjust = 1.25, ...) +
            text(x = 0.5, y = 0, label = names[2], vjust = 1.25, ...)
    }

    p
}

text <- function(..., size = 4.5) {
    annotate(
        "richtext", size = size, fill = NA, label.color = NA, ...
    )
}

textbox <- function(..., size = 4.5) {
    annotate(
        "richtext", size = size, label.size = 0.5, ...
    )
}

label <- function(realsfs, winsfs, colours, sep = "<br>") {
    paste0(
        "<span style='color:", colours[1], "'>", realsfs, "</span>",
        sep,
        "<span style='color:", colours[2], "'>", winsfs, "</span>"
    )
}

time_arrow <- function(y, yend, ...) {
    annotate(
        "segment",
        x = 0.65,
        y = y,
        xend = 0.65,
        yend = yend,
        arrow = arrow(ends = "both", length = unit(2, "mm")),
        size = 0.5,
        colour = "grey30",
        ...
    )
}

add_constant_time <- function(p, label = NULL, text_x = 0.725, ...) {
    p <- p + time_arrow(y = 0.0, yend = 1.1)

    if (!is.null(label)) {
        p <- p + text(x = text_x, y = 0.55, label = label, ...)
    }

    p
}

add_split_time <- function(p, labels = NULL, ...) {
    p <- p +
        time_arrow(y = 0.0, yend = 0.5) +
        time_arrow(y = 0.6, yend = 1.1)

    if (!is.null(labels)) {
        p <- p +
            text(x = 0.725, y = 0.85, label = labels[1], ...) +
            text(x = 0.725, y = 0.25, label = labels[2], ...)
    }

    p
}

add_shared_size <- function(p, label, ...) {
    p + text(x = 0.3, y = 1.2, label = label, ...)
}

add_split_sizes <- function(p, labels, y = 0.5, ...) {
    p +
        text(x = 0.1, y = y, label = labels[1], ...) +
        text(x = 0.5, y = y, label = labels[2], ...)

}

migration_arrow <- function(y, ends, ...) {
    annotate(
        "segment",
        x = 0.25, y = y, xend = 0.35, yend = y,
        arrow = arrow(ends = ends, length = unit(2, "mm")),
        colour = "grey30",
        size = 0.5,
        ...
    )
}

add_symmetric_migration <- function(p, label = NULL, y = 0.5, ...) {
    p <- p + migration_arrow(y = y, ends = "both")

    if (!is.null(label)) {
        p <- p + text(x = 0.3, y = y + 0.05, label = label, vjust = -0.2, ...)
    }

    p
}

add_asymmetric_migration <- function(p, labels = NULL, y = 0.5, ...) {
    p <- p +
        migration_arrow(y = y - 0.025, ends = "first") +
        migration_arrow(y = y + 0.025, ends = "last")

    if (!is.null(labels)) {
        p <- p +
            text(x = 0.3, y = y - 0.05, label = labels[1], vjust = 1.2, ...) +
            text(x = 0.3, y = y + 0.05, label = labels[2], vjust = -0.2, ...)
    }

    p
}

add_log_likelihood <- function(p, label, ...) {
    p +
        text(
            x = 0.425,
            y = 1.25,
            hjust = 0,
            vjust = 0,
            label = "Log-likelihood:",
            size = 4,
            ...
        ) +
        text(
            x = 0.425,
            y = 1.25,
            label = label,
            hjust = -0.2,
            vjust = 0.8,
            ...
        )
}
