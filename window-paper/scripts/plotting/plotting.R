suppressPackageStartupMessages({
    library(cetcolor)
    library(colorspace)
    library(dplyr)
    library(ggh4x)
    library(ggplot2)
    library(ggrepel)
    library(ggtext)
    library(glue)
    library(patchwork)
    library(readr)
    library(ragg)
    library(scales)
    library(tidyr)
})

base_theme <- theme_minimal(base_size = 12) %+replace%
    theme(
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        plot.title = element_markdown(
            size = 14,
            hjust = 0.0,
            margin = margin(0, 0, 0, 0)
        ),
        strip.text = element_markdown(
            size = 14,
            hjust = 0.0,
            margin = margin(0, 0, 0, 0)
        ),
        axis.title.x = element_markdown(
            size = 12,
            margin = margin(t = 4)
        ),
        axis.text.x = element_text(
            size = 10,
            colour = "grey30",
            margin = margin(t = 0)
        ),
        axis.title.y = element_markdown(
            size = 12,
            angle = 90,
            margin = margin(r = 2)
        ),
        axis.title.y.right = element_markdown(
            size = 12,
            angle = 270,
            margin = margin(l = 2)
        ),
        axis.text.y = element_text(
            size = 10,
            colour = "grey30",
            angle = 90,
            vjust = 0,
            hjust = 0.5,
            margin = margin(r = 0)
        ),
        axis.text.y.right = element_text(
            hjust = 0.5
        ),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(
            colour = "grey10",
            size = 0.1
        ),
        legend.title = element_markdown(size = 12),
        legend.text = element_text(
            size = 10,
            colour = "grey30"
        ),
        legend.margin = margin(0, 0, 0, 0, "pt"),
    )

theme_set(base_theme)

normalise <- function(x) x / sum(x)

condition_labeller <- function(x, sep = " ") {
    case_when(
        x == "realsfs" ~ "realSFS",
        x == "winsfs_b500_w100" ~ paste("winsfs", "(W=100)", sep = sep),
        x == "winsfs_b500_w250" ~ paste("winsfs", "(W=250)", sep = sep),
        x == "winsfs_b500_w500" ~ paste("winsfs", "(W=500)", sep = sep),
    )
}

label_sci <- function(x, digits = 3) {
    if_else(
        x == 0,
        "0",
        gsub("(?<=e)([-+])?[0]+", "\\1", scientific(x, digits), perl = TRUE)
    )
}

fmt_populations <- function(x, sep = " ") {
    fmt_pops <- c(
        "YRI" = paste0("Human 1D", sep, "(YRI)"),
        "CEU-YRI" = paste0("Human 2D", sep, "(CEU/YRI)"),
        "MasaiMara" = paste0("Impala 1D", sep, "(Maasai Mara)"),
        "Shangani-MasaiMara" = paste0("Impala 2D", sep, "(Shangani/Maasai Mara)")
    )
    ordered(fmt_pops[x], fmt_pops)
}

label_k <- function(x, accuracy = 1) {
    comma(x, accuracy = accuracy, scale = 1/1000, suffix = "K")
}
