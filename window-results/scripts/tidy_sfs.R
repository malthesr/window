#!/usr/bin/env Rscript

source("scripts/sfs.R")

tidy_sfs <- function(sfs) {
    if (is.null(dim(sfs))) {
        tibble(i = seq_along(sfs) - 1L, x = sfs)
    } else if (length(dim(sfs)) == 2) {
        expand_grid(
            j = seq_len(ncol(sfs)) - 1L,
            i = seq_len(nrow(sfs)) - 1L
        ) |>
        mutate(x = c(sfs))
    } else {
        stop("unexpected shape in tidy")
    }
}

tajimas_theta <- function(sfs) {
    stopifnot(is.null(dim(sfs)))
    sfs <- normalise(sfs)
    n <- length(sfs) - 1
    seg_i <- 1:(n - 1)
    seg <- sfs[seg_i + 1]
    sum(seg_i * (n - seg_i) * seg) / choose(n, 2)
}

f2 <- function(sfs) {
    stopifnot(length(dim(sfs)) == 2)
    sfs <- normalise(sfs)
    n1 <- nrow(sfs) - 1
    n2 <- ncol(sfs) - 1
    p1 <- 0:n1 / n1
    p2 <- 0:n2 / n2
    f <- 0
    for (i in seq_along(p1)) {
        for (j in seq_along(p2)) {
            add <- (p1[i] - p2[j])^2 * sfs[i, j]
            f <- f + add
        }
    }

    f
}

hudsons_fst <- function(sfs) {
    stopifnot(length(dim(sfs)) == 2)
    mat <- sfs
    mat[c(1, length(mat))] <- 0
    mat <- normalise(mat)
    n1 <- nrow(mat) - 1
    n2 <- ncol(mat) - 1
    p1 <- 0:n1 / n1
    p2 <- 0:n2 / n2
    f_num <- function(p1, p2) {
        q1 <- 1 - p1
        q2 <- 1 - p2
        (p1 - p2)^2 - p1 * q1 / (n1 - 1) - p2 * q2 / (n2 - 1)
    }
    num <- outer(p1, p2, f_num)
    f_denom <- function(p1, p2) {
        q1 <- 1 - p1
        q2 <- 1 - p2
        p1 * q2 + p2 * q1
    }
    denom <- outer(p1, p2, f_denom)
    sum(mat * num) / sum(mat * denom)
}

sfs_stats <- function(sfs) {
    if (is.null(dim(sfs))) {
        list(theta = tajimas_theta(sfs))
    } else if (length(dim(sfs)) == 2) {
        list(fst = hudsons_fst(sfs), f2 = f2(sfs))
    }
}

epochs <- snakemake@params[["epochs"]]
sfs_paths <- snakemake@input[["sfs"]]
sfs <- lapply(sfs_paths, read_sfs)

sfs_df <- sfs |>
    lapply(tidy_sfs) |>
    setNames(epochs) |>
    bind_rows(.id = "epoch")

stat_df <- sfs |>
    lapply(sfs_stats) |>
    setNames(epochs) |>
    bind_rows(.id = "epoch")

sfs_csv <- snakemake@output[["sfs_csv"]]
write_csv(sfs_df, sfs_csv)

stat_csv <- snakemake@output[["stat_csv"]]
write_csv(stat_df, stat_csv)
