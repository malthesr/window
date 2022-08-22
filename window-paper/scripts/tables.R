#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(kableExtra)
    library(knitr)
    library(readr)
})

options(scipen = 999)

rd <- function(x, digits) {
    x <- round(x, digits)
    sprintf(paste0("%.", digits, "f"), x)
}

latex_sci <- function(x, signif) {
    x <- signif(x, signif)
    x <- formatC(x, format = "e", digits = signif - 1)
    paste0("\\scinum{", x, "}")
}

latex_num <- function(x, digits) {
    x <- rd(x, digits = digits)
    paste0("\\num{", x, "}")
}

latex_depth <- function(x, digits) {
    paste0("\\SI{",  rd(x, digits) , "}{\\depth}")
}

latex_depth_range <- function(from, to, digits) {
    from <- rd(from, digits)
    to <- rd(to, digits)
    paste0("\\SIrange{", from, "}{", to, "}{\\depth}")
}

make_table <- function(df) {
    df |>
        transmute(
            "Population" = population,
            "Individuals" = individuals,
            "Sites" = latex_sci(sites, 3),
            "Median depth (range)" = paste0(
                latex_depth(depth_median, 1),
                " (",
                latex_depth_range(depth_min, depth_max, 1),
                ")"
            ),
            "Contigs $\\geq \\SI{100}{\\kilo\\bases}$" = contigs,
            "$n_{50}$" = latex_sci(n50, 2),
            "$\\fst$" = rd(fst, 2),
        ) |>
        kable(
            format = "latex",
            booktabs = TRUE,
            escape = FALSE,
            align = "lrrrrrr"
        ) |>
        pack_rows("Human", 1, 2) |>
        pack_rows("Impala", 3, 4) |>
        collapse_rows(5:7, latex_hline = "none")
}

df <- read_csv("data/data.csv", col_types = "cciidddiid")

df |>
    filter(split == "train") |>
    make_table() |>
    write_lines(file = "tables/data_train.tex")

df |>
    filter(split == "test") |>
    make_table() |>
    write_lines(file = "tables/data_test.tex")

