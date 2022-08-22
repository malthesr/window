SIM_NAMES <- paste("Population", c("A", "B"))

read_sim_sfs_csv <- function(p, peak) {
    read_csv(p, col_types = "liiciiid") |>
        filter(peak == !!peak) |>
        select(-peak)
}

read_sim_stop_df <- function(p, peak, tolerance = 1e-4) {
    block_loglik_df <- read_csv(p, col_types = "liicidd") |>
        filter(peak == !!peak) |>
        select(-peak)
    block_loglik_df |>
        group_by(n, d, condition) |>
        filter(change <= tolerance) |>
        slice_head(n = 1)
}

plot_text <- function(s, angle = 0) {
    ggplot(data.frame(s = s, x = 1, y = 1)) +
        geom_text(aes(x, y, label = s), angle = angle) +
        theme_void() +
        coord_cartesian(clip = "off")
}

plot_d <- function(d, angle = 270) {
    plot_text(paste("Mean depth", d), angle = angle)
}

plot_n <- function(n) {
    plot_text(paste(n, "individuals"))
}

remove_axis_labels <- function(lst, nrow) {
    n <- length(lst)
    stopifnot(n %% nrow == 0)
    ncol <- n %/% nrow

    for (i in 1:n) {
        is_first_col <- (i - 1) %% ncol == 0
        if (!is_first_col) {
            lst[[i]] <- lst[[i]] +
                theme(
                    axis.title.y = element_blank(),
                )
        }

        is_last_row <- i > n - ncol
        if (!is_last_row) {
            lst[[i]] <- lst[[i]] +
                theme(
                    axis.title.x = element_blank(),
                )
        }
    }

    lst
}

pseudo_facet_1x3_n <- function(lst, n = c(5, 10, 20)) {
    wrap_plots(
        plot_n(n[1]), plot_n(n[2]), plot_n(n[3]),
        lst[[1]], lst[[2]], lst[[3]],
        ncol = 3,
        heights = c(0.10, 1),
        widths = c(1, 1, 1)
    )
}

pseudo_facet_1x3_d <- function(lst, d = c(2, 4, 8)) {
    wrap_plots(
        plot_d(d[1], 0), plot_d(d[2], 0), plot_d(d[3], 0),
        lst[[1]], lst[[2]], lst[[3]],
        ncol = 3,
        heights = c(0.10, 1),
        widths = c(1, 1, 1)
    )
}

pseudo_facet_3x3 <- function(lst, d = c(2, 4, 8), n = c(5, 10, 20)) {
    wrap_plots(
        plot_n(n[1]), plot_n(n[2]), plot_n(n[3]), plot_spacer(),
        lst[[1]],  lst[[2]],  lst[[3]],  plot_d(d[1]),
        lst[[4]],  lst[[5]],  lst[[6]],  plot_d(d[2]),
        lst[[7]],  lst[[8]],  lst[[9]],  plot_d(d[3]),
        nrow = 4,
        heights = c(0.15, 1, 1, 1),
        widths = c(1, 1, 1, 0.15)
    )
}
