suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

parse_sfs_header <- function(header) {
    matches <- regmatches(header, gregexpr("[[:digit:]]+", header))
    as.numeric(unlist(matches))
}

parse_flat_sfs <- function(flat, shape) {
    split <- as.numeric(strsplit(flat, " ")[[1]])
    n <- length(split)
    dims <- length(shape)
    if (dims == 1) {
        split
    } else if (dims == 2) {
        stopifnot(prod(shape) == n)
        matrix(split, nrow = shape[1], ncol = shape[2], byrow = TRUE)
    } else {
        stop("unexpected shape in parse")
    }
}

read_sfs <- function(p) {
    lines <- read_lines(p)
    header <- lines[1]
    shape <- parse_sfs_header(header)
    flat <- lines[2]
    parse_flat_sfs(flat, shape)
}

normalise <- function(x) {
    x / sum(x)
}
