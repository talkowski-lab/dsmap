#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot DEL and DUP DSMap observed/expected values (O/Es) in 5kb bins
# vs. gnocchi (gnomAD noncoding z)


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly = TRUE)
require(matrixStats, quietly = TRUE)
require(optparse, quietly = TRUE)
require(stringr, quietly = TRUE)

# Set global options and constants
options(stringsAsFactors = FALSE, scipen = 1000)
dsmapR::load.constants(c("colors", "scales"))


##################
# Data functions #
##################
# Load data
load.scores <- function(del.tsv, dup.tsv, constraint.tsv, n.quantiles = 100) {
    # Load DEL and DUP data files
    del <- read.table(del.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(del)[1] <- "bin"
    dup <- read.table(dup.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(dup)[1] <- "bin"

    # Infer bin size
    del_coords <- as.data.frame(str_split_fixed(del[, 1], "_", 3))
    for (i in c(2, 3)) {
        del_coords[, i] <- as.numeric(del_coords[, i])
    }
    binsize_del <- infer.bin.size(del_coords)
    dup_coords <- as.data.frame(str_split_fixed(dup[, 1], "_", 3))
    for (i in c(2, 3)) {
        dup_coords[, i] <- as.numeric(dup_coords[, i])
    }
    binsize_dup <- infer.bin.size(dup_coords)
    if (binsize_del != binsize_dup) {
        stop("Bin sizes in DEL and DUP files are not identical")
    } else {
        binsize <- binsize_del
    }

    # Load noncoding constraint data file
    constraint.scores <- read.table(constraint.tsv, header = TRUE, sep = "\t", comment.char = "")

    # Assign each noncoding constraint 1kb bin to the corresponding DSMap 5kb bin
    constraint.scores$bin <- paste(
        constraint.scores$chrom,
        floor(constraint.scores$start / binsize) * binsize,
        ceiling(constraint.scores$start / binsize) * binsize,
        sep = "_"
    )
    constraint.scores <- constraint.scores[, c("element_id", "bin", "expected", "observed", "oe", "z")]

    # Compute different aggregations of noncoding constraint scores within each 5kb bin:
    # oe_overall = total obs / total exp over 5kb bin
    # oe_min = min O/E
    # oe_average = average O/E
    constraint.bin.scores.agg <- by(constraint.scores, constraint.scores$bin, function(x) {
        c(
            bin = unique(x$bin),
            noncoding_oe_overall = sum(x$observed) / sum(x$expected),
            noncoding_oe_min = min(x$oe),
            noncoding_oe_average = mean(x$oe)
        )
    })
    # Collect aggregation results in dataframe
    constraint.bin.scores <- as.data.frame(
        do.call("rbind", constraint.bin.scores.agg),
        stringsAsFactors = FALSE
    )
    for (i in 2:ncol(constraint.bin.scores)) {
        constraint.bin.scores[, i] <- as.numeric(constraint.bin.scores[, i])
    }

    # # Remove rows where mu == na.val or mu is infinite
    # # (na.val is introduced by athena as a placeholder for situations where
    # #  mutation rates are missing)
    # # TODO: Amend these in the model
    # mu <- mu[!(mu$mu == na.val) & !(is.infinite(mu$mu)), ]

    # Merge data
    bin.scores <- merge(del, dup,
        by = "bin",
        suffixes = c(".DEL", ".DUP"), all = FALSE, sort = FALSE
    )
    bin.scores <- merge(bin.scores, constraint.bin.scores, by = "bin", all = FALSE, sort = FALSE)
    print(
        paste(
            "Number of DSMap bins lost in merge with constraint scores:",
            nrow(del) - nrow(bin.scores)
        )
    )

    # Calculate quantiles of constraint score
    for (f in colnames(constraint.bin.scores)[-1]) {
        bin.scores[, paste(f, "quantile", sep = "_")] <- ceiling(n.quantiles *
            rank(bin.scores[, f], na.last = "keep", ties.method = "random") /
            nrow(bin.scores))
    }

    # Compute expected # of bin dels & dups given sample size
    n.samples <- 63046
    bin.scores$exp.DEL <- (10^bin.scores$mu.DEL) * 2 * n.samples
    bin.scores$exp.DUP <- (10^bin.scores$mu.DUP) * 2 * n.samples

    # Compute bin DEL and DUP O/E
    bin.scores$oe.DEL <- bin.scores$n_svs.DEL / bin.scores$exp.DEL
    bin.scores$oe.DUP <- bin.scores$n_svs.DUP / bin.scores$exp.DUP
    return(bin.scores)
}

# Compute bin O/E summary statistics with confidence intervals per constraint score quantile
compute.oe.stats.by.quantile <- function(bin.scores, n.quantiles = 100, bootstrap = FALSE, n.bootstraps = 1000, seed = 42) {
    # Retrieve all noncoding constraint quantile types
    constraint.aggs <- colnames(bin.scores)[grepl("quantile", colnames(bin.scores))]

    oe.stats <- do.call("rbind", lapply(constraint.aggs, function(f) {
        f.oe.stats <- cbind(1:n.quantiles, do.call("rbind", lapply(1:n.quantiles, function(i) {
            # Get bins in quantile
            idxs <- which(bin.scores[, f] == i)
            if (bootstrap) {
                # Bootstrap sampling of bins in quantile
                set.seed(seed)
                bootstrap.idxs <- replicate(n.bootstraps, sample(idxs, length(idxs), replace = TRUE))

                # Compute obs and exp values with bootstrap for DEL and DUP
                sumstats <- do.call("cbind", lapply(c("DEL", "DUP"), function(cnv) {
                    # Get obs and exp for bins in quantile and bootstrap samples
                    obs <- bin.scores[idxs, paste("n_svs", cnv, sep = ".")]
                    exp <- bin.scores[idxs, paste("exp", cnv, sep = ".")]
                    bootstrap.obs <- matrix(
                        bin.scores[bootstrap.idxs, paste("n_svs", cnv, sep = ".")],
                        ncol = n.bootstraps
                    )
                    bootstrap.exp <- matrix(
                        bin.scores[bootstrap.idxs, paste("exp", cnv, sep = ".")],
                        ncol = n.bootstraps
                    )

                    # Compute bootstrap confidence intervals on quantile summary statistics
                    bootstrap.medians <- colMedians(bootstrap.obs / bootstrap.exp)
                    median.ci <- quantile(bootstrap.medians, c(0.05, 0.95))
                    bootstrap.sums <- colSums(bootstrap.obs) / colSums(bootstrap.exp)
                    sums.ci <- quantile(bootstrap.sums, c(0.05, 0.95))

                    # Return quantile sum stat point estimates and CIs
                    data.frame(
                        median(obs / exp), median.ci[1], median.ci[2],
                        sum(obs) / sum(exp), sums.ci[1], sums.ci[2],
                        mad(obs / exp)
                    )
                }))
            } else {
                # Compute obs and exp values without bootstrap for DEL and DUP
                sumstats <- do.call("cbind", lapply(c("DEL", "DUP"), function(cnv) {
                    # Get obs and exp for bins in quantile and bootstrap samples
                    obs <- bin.scores[idxs, paste("n_svs", cnv, sep = ".")]
                    exp <- bin.scores[idxs, paste("exp", cnv, sep = ".")]
                    # Return quantile sum stat point estimates
                    data.frame(median(obs / exp), sum(obs) / sum(exp), mad(obs / exp))
                }))
            }
            return(sumstats)
        })))
        return(cbind(f, f.oe.stats))
    }))

    # Add column names
    if (bootstrap) {
        sumstats <- c(
            "median", "median.lower", "median.upper",
            "sum", "sum.lower", "sum.upper", "mad"
        )
    } else {
        sumstats <- c("median", "sum", "mad")
    }
    colnames(oe.stats) <- c("constraint.agg", "quantile", do.call("c", lapply(c("DEL", "DUP"), function(cnv) {
        paste(sumstats, cnv, sep = ".")
    })))
    return(oe.stats)
}

# Compute correlation coefficients and p-values for each relationship
compute.cor <- function(oe.stats) {
    constraint.aggs <- unique(oe.stats$constraint.agg)
    cors <- do.call("rbind", lapply(constraint.aggs, function(f) {
        f.cors <- do.call("rbind", lapply(colnames(oe.stats)[3:ncol(oe.stats)], function(s) {
            cor.vals <- cor.test(
                oe.stats[oe.stats$constraint.agg == f, "quantile"],
                oe.stats[oe.stats$constraint.agg == f, s],
                method = "spearman", exact = FALSE
            )
            return(data.frame(s, cor.vals$estimate, cor.vals$p.value))
        }))
        return(cbind(f, f.cors))
    }))
    colnames(cors) <- c("constraint.agg", "sumstat.type", "estimate", "p")
    return(cors)
}

# Compute robust Z-scores of O/Es for bins within each constraint score quantile
compute.z.by.quantile <- function(bin.scores, quantile.stats) {
    zs.by.cnv <- lapply(c("DEL", "DUP", "DUP.CG"), function(cnv) {
        # Annotate bins with quantile-level median and MAD O/E
        scores <- merge(
            bin.scores[, c(
                "bin", "quantile",
                paste("n_svs", cnv, sep = "."),
                paste("exp", cnv, sep = "."),
                paste("oe", cnv, sep = ".")
            )], quantile.stats[, c(
                "quantile",
                paste("median", cnv, sep = "."),
                paste("mad", cnv, sep = ".")
            )], "quantile"
        )
        colnames(scores)[(ncol(scores) - 1):ncol(scores)] <- c(
            paste("quantile.median", cnv, sep = "."),
            paste("quantile.mad", cnv, sep = ".")
        )

        # Calculate bin robust Z-scores
        scores[, paste("quantile.z", cnv, sep = ".")] <- (scores[, paste("oe", cnv, sep = ".")] -
            scores[, paste("quantile.median", cnv, sep = ".")]) /
            scores[, paste("quantile.mad", cnv, sep = ".")]
        return(scores)
    })
    return(Reduce(
        function(x, y) merge(x, y, by = c("bin", "quantile"), all = FALSE, sort = FALSE),
        zs.by.cnv
    ))
}

# train 0135 from boston on saturday july 13 at 1:35pm
# and arrive at 5:52pm


######################
# Plotting functions #
######################
# Plot O/Es of bins by their constraint score quantile
plot.constraint.quantile.oe <- function(dat, cors, constraint.label, bootstrap = FALSE) {
    par(mfrow = c(2, 2), mar = c(2, 3.5, 2, 1.5))

    pch_type <- list("sum" = 19, "median" = 21)
    sumstat_title <- list("sum" = "Quantile-wide Mean", "median" = "Median Bin")

    for (sumstat in c("sum", "median")) {
        for (cnv in c("DEL", "DUP")) {
            plot(
                dat[, paste(sumstat, cnv, sep = ".")],
                pch = pch_type[[sumstat]],
                col = get(paste(ifelse(cnv == "DUP.CG", "DUP", cnv), "colors", sep = "."))$main,
                xaxt = "n", xlab = "", ylab = "", las = 2,
                main = paste(cnv, " Obs/Exp ", "(", sumstat_title[[sumstat]], ")", sep = ""),
                panel.first = c(abline(h = 1, lty = 5))
            )
            if (bootstrap) {
                polygon(
                    x = c(1:nrow(dat), rev(1:nrow(dat))),
                    y = c(
                        dat[, paste(sumstat, "upper", cnv, sep = ".")],
                        rev(dat[, paste(sumstat, "lower", cnv, sep = ".")])
                    ),
                    col = adjustcolor(
                        get(paste(ifelse(cnv == "DUP.CG", "DUP", cnv), "colors", sep = "."))$main,
                        alpha.f = 0.3
                    ),
                    border = NA
                )
            }
            title(
                main = paste(
                    "Spearman rho =",
                    round(cors[cors$sumstat.type == paste(sumstat, cnv, sep = "."), "estimate"], 3),
                    "; p =",
                    formatC(cors[cors$sumstat.type == paste(sumstat, cnv, sep = "."), "p"],
                        format = "e", digits = 1
                    )
                ), outer = FALSE, line = -1, cex.main = 1, font.main = 1
            )
            mtext(1, text = paste("5kb bins grouped by", constraint.label, "quantile"))
            mtext(2, line = 2.5, text = "Observed / Expected")
        }
    }
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
    make_option("--n-quantiles",
        type = "integer", default = 100,
        help = "Number of quantiles to divide bins into [default %default]"
    )
)

# Get command-line arguments & options
arg_list <- c(
    "del.tsv", "dup.tsv", "constraint.tsv", "out.prefix"
    # "out_bins_del.tsv", "out_bins_dup.tsv",
)
args <- parse_args(
    OptionParser(
        usage = paste("%prog", paste0(arg_list, collapse = " ")),
        option_list = option_list
    ),
    positional_arguments = TRUE
)
opts <- args$options

# Checks for appropriate positional arguments
if (length(args$args) != length(arg_list)) {
    stop(
        paste(
            length(arg_list), "positional arguments required:",
            paste0(arg_list, collapse = ", ")
        )
    )
}

# Writes args & opts to vars
del.tsv <- args$args[1]
dup.tsv <- args$args[2]
constraint.tsv <- args$args[3]
out.prefix <- args$args[4]
# out.bins.del.tsv <- args$args[5]
# out.bins.dup.tsv <- args$args[6]
n.quantiles <- opts$`n-quantiles`

# Load bin obs, mus, exps, and constraint scores
bin.scores <- load.scores(del.tsv, dup.tsv, constraint.tsv, n.quantiles)

# Compute bin O/E summary statistics per constraint score quantile
oe.stats <- compute.oe.stats.by.quantile(bin.scores, n.quantiles, bootstrap = FALSE)

# Compute correlations for each O/E summary statistic type with constraint score percentile
cors <- compute.cor(oe.stats)

constraint.aggs <- unique(oe.stats$constraint.agg)
# Plot O/Es of bins by their constraint score quantile
for (f in constraint.aggs) {
    pdf(paste(out.prefix, f, "pdf", sep = "."), height = 5, width = 10)
    plot.constraint.quantile.oe(
        oe.stats[oe.stats$constraint.agg == f, ], cors[cors$constraint.agg == f, ],
        gsub("_quantile", "", f),
        bootstrap = FALSE
    )
    dev.off()
}

