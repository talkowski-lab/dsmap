#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot gene DEL and DUP observed/expected values (O/Es) by rCNV pHaplo or pTriplo bin
# and write out genes ranked by discordance of O/Es vs. rCNV scores


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly = TRUE)
require(matrixStats, quietly = TRUE)
require(optparse, quietly = TRUE)

# Set global options and constants
options(stringsAsFactors = FALSE, scipen = 1000)
dsmapR::load.constants(c("colors"))


##################
# Data functions #
##################
# Load data
load.gene.scores <- function(del.tsv, dup.tsv, rcnv.tsv) {
    del <- read.table(del.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(del)[1] <- "gene"
    dup <- read.table(dup.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(dup)[1] <- "gene"
    n.samples <- 62989
    rcnv.scores <- read.table(rcnv.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(rcnv.scores)[1] <- "gene"

    # Merge data
    gene.scores <- merge(del, dup, by = "gene", all = FALSE, suffixes = c(".DEL", ".DUP"), sort = FALSE)
    gene.scores <- merge(gene.scores, rcnv.scores, by = "gene", all = FALSE, sort = FALSE)

    # Percentiles of pHaplo & pTriplo
    gene.scores$DEL.pct <- ceiling(100 * rank(gene.scores$pHaplo) / nrow(gene.scores))
    gene.scores$DUP.pct <- ceiling(100 * rank(gene.scores$pTriplo) / nrow(gene.scores))

    # Compute expected # of dels & dups given sample size
    gene.scores$exp.DEL <- (10^gene.scores$mu.DEL) * 2 * n.samples
    gene.scores$exp.DUP <- (10^gene.scores$mu.DUP) * 2 * n.samples

    # Compute obs/exp
    gene.scores$oe.DEL <- gene.scores$n_svs.DEL / gene.scores$exp.DEL
    gene.scores$oe.DUP <- gene.scores$n_svs.DUP / gene.scores$exp.DUP
    return(gene.scores)
}

# Compute gene O/E summary statistics with confidence intervals per rCNV score bin
compute.oe.stats.by.rcnv.bin <- function(gene.scores, n.bootstraps = 1000, seed = 42) {
    oe.stats <- do.call("cbind", lapply(c("DEL", "DUP"), function(cnv) {
        do.call("rbind", lapply(1:100, function(i) {
            idxs <- which(gene.scores[, paste(cnv, "pct", sep = ".")] == i)
            obs <- gene.scores[idxs, paste("n_svs", cnv, sep = ".")]
            exp <- gene.scores[idxs, paste("exp", cnv, sep = ".")]

            # Bootstrap sampling of genes
            set.seed(seed)
            bootstrap.idxs <- replicate(n.bootstraps, sample(idxs, length(idxs), replace = TRUE))

            # Get obs and exp per gene in bootstrap samples
            bootstrap.obs <- matrix(
                gene.scores[bootstrap.idxs, paste("n_svs", cnv, sep = ".")],
                ncol = n.bootstraps
            )
            bootstrap.exp <- matrix(
                gene.scores[bootstrap.idxs, paste("exp", cnv, sep = ".")],
                ncol = n.bootstraps
            )

            # Compute confidence intervals on bin summary statistics
            bootstrap.medians <- colMedians(bootstrap.obs / bootstrap.exp)
            median.ci <- quantile(bootstrap.medians, c(0.05, 0.95))
            bootstrap.sums <- colSums(bootstrap.obs) / colSums(bootstrap.exp)
            sums.ci <- quantile(bootstrap.sums, c(0.05, 0.95))
            data.frame(
                median(obs / exp), median.ci[1], median.ci[2],
                sum(obs) / sum(exp), sums.ci[1], sums.ci[2]
            )
        }))
    }))
    colnames(oe.stats) <- c(
        "median.DEL", "median.DEL.lower", "median.DEL.upper",
        "sum.DEL", "sum.DEL.lower", "sum.DEL.upper",
        "median.DUP", "median.DUP.lower", "median.DUP.upper",
        "sum.DUP", "sum.DUP.lower", "sum.DUP.upper"
    )
    return(oe.stats)
}

# Compute correlation coefficients and p-values for each relationship
compute.cor <- function(oe.stats) {
    cors <- do.call("rbind", lapply(seq_along(oe.stats), function(i) {
        cor.vals <- cor.test(
            seq_along(oe.stats[, i]), unlist(oe.stats[, i]),
            method = "spearman"
        )
        data.frame(names(oe.stats)[i], cor.vals$estimate, cor.vals$p.value)
    }))
    colnames(cors) <- c("sum.type", "estimate", "p")
    return(cors)
}

# Compute robust Z-scores of O/Es for genes within each rCNV score bin
compute.z.by.rcnv.bin <- function(gene.scores) {
    zs.by.cnv <- lapply(c("DEL", "DUP"), function(cnv) {
        # Calculate median and MAD O/Es per bin
        bin.stats <- cbind(pct = 1:100, do.call("rbind", lapply(
            1:100,
            function(i) {
                oes <- gene.scores[
                    gene.scores[, paste(cnv, "pct", sep = ".")] == i,
                    paste("oe", cnv, sep = ".")
                ]
                data.frame(bin.median.oe = median(oes), bin.mad.oe = mad(oes))
            }
        )))
        colnames(bin.stats) <- paste(cnv, colnames(bin.stats), sep = ".")

        # Annotate genes with bin-level median and MAD O/E
        scores <- merge(
            gene.scores[, c(
                "gene",
                paste(cnv, "pct", sep = "."),
                paste("n_svs", cnv, sep = "."),
                paste("exp", cnv, sep = "."),
                paste("oe", cnv, sep = ".")
            )], bin.stats,
            by = paste(cnv, "pct", sep = ".")
        )

        # Calculate gene robust Z-scores
        scores <- cbind(
            scores,
            (scores[, paste("oe", cnv, sep = ".")] -
                scores[, paste(cnv, "bin.median.oe", sep = ".")]) /
                scores[, paste(cnv, "bin.mad.oe", sep = ".")]
        )
        colnames(scores)[length(colnames(scores))] <- paste(cnv, "bin.z", sep = ".")
        return(scores)
    })
    return(merge(zs.by.cnv[[1]], zs.by.cnv[[2]], by = "gene"))
}


######################
# Plotting functions #
######################
# Plot O/Es of genes by their rCNV score bin
plot.rcnv.bin.oe <- function(dat, cors) {
    par(mfrow = c(2, 2), mar = c(2, 3.5, 2, 1.5))
    plot(
        dat$sum.DEL,
        pch = 19, col = DEL.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DEL Obs/Exp (Bin-wide Mean)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    polygon(
        x = c(1:100, rev(1:100)),
        y = c(dat$sum.DEL.upper, rev(dat$sum.DEL.lower)),
        col = adjustcolor(DEL.colors$main, alpha.f = 0.3), border = NA
    )
    title(
        main = paste(
            "Spearman rho =", round(cors[cors$sum.type == "sum.DEL", "estimate"], 3),
            "; p =", formatC(cors[cors$sum.type == "median.DUP", "p"], format = "e", digits = 1)
        ), outer = FALSE, line = -1, cex.main = 1, font.main = 1
    )
    mtext(1, text = "Genes Binned by pHaplo Percentile")
    mtext(2, line = 2.5, text = "Observed / Expected")

    plot(
        dat$sum.DUP,
        pch = 19, col = DUP.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DUP Obs/Exp (Bin-wide Mean)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    polygon(
        x = c(1:100, rev(1:100)),
        y = c(dat$sum.DUP.upper, rev(dat$sum.DUP.lower)),
        col = adjustcolor(DUP.colors$main, alpha.f = 0.3), border = NA
    )
    title(
        main = paste(
            "Spearman rho =", round(cors[cors$sum.type == "sum.DUP", "estimate"], 3),
            "; p =", formatC(cors[cors$sum.type == "sum.DUP", "p"], format = "e", digits = 1)
        ), outer = FALSE, line = -1, cex.main = 1, font.main = 1
    )
    mtext(1, text = "Genes Binned by pTriplo Percentile")
    mtext(2, line = 2.5, text = "Observed / Expected")

    plot(
        dat$median.DEL,
        pch = 21, col = DEL.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DEL Obs/Exp (Median Gene)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    polygon(
        x = c(1:100, rev(1:100)),
        y = c(dat$median.DEL.upper, rev(dat$median.DEL.lower)),
        col = adjustcolor(DEL.colors$main, alpha.f = 0.3), border = NA
    )
    title(
        main = paste(
            "Spearman rho =", round(cors[cors$sum.type == "median.DEL", "estimate"], 3),
            "; p =", formatC(cors[cors$sum.type == "median.DEL", "p"], format = "e", digits = 1)
        ), outer = FALSE, line = -1, cex.main = 1, font.main = 1
    )
    mtext(1, text = "Genes Binned by pHaplo Percentile")
    mtext(2, line = 2.5, text = "Observed / Expected")

    plot(
        dat$median.DUP,
        pch = 21, col = DUP.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DUP Obs/Exp (Median Gene)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    polygon(
        x = c(1:100, rev(1:100)),
        y = c(dat$median.DUP.upper, rev(dat$median.DUP.lower)),
        col = adjustcolor(DUP.colors$main, alpha.f = 0.3), border = NA
    )
    title(
        main = paste(
            "Spearman rho =", round(cors[cors$sum.type == "median.DUP", "estimate"], 3),
            "; p =", formatC(cors[cors$sum.type == "median.DUP", "p"], format = "e", digits = 1)
        ), outer = FALSE, line = -1, cex.main = 1, font.main = 1
    )
    mtext(1, text = "Genes Binned by pTriplo Percentile")
    mtext(2, line = 2.5, text = "Observed / Expected")
}


###########
# RScript #
###########
# List of command-line options
option_list <- list()
# TODO: Parameterize CNV type, rCNV score type, bin quantile size?

# Get command-line arguments & options
args <- parse_args(
    OptionParser(
        usage = "%prog del.tsv dup.tsv rcnv.tsv out.pdf out_genes_del.tsv out_genes_dup.tsv",
        option_list = option_list
    ),
    positional_arguments = TRUE
)
opts <- args$options

# Checks for appropriate positional arguments
if (length(args$args) != 6) {
    stop(
        paste(
            "Six positional arguments required:",
            "del.tsv, dup.tsv, rcnv.tsv, out.pdf, out_genes_del.tsv, and out_genes_dup.tsv"
        )
    )
}

# Writes args & opts to vars
del.tsv <- args$args[1]
dup.tsv <- args$args[2]
rcnv.tsv <- args$args[3]
out.pdf <- args$args[4]
out_genes_del.tsv <- args$args[5]
out_genes_dup.tsv <- args$args[6]

# Load gene obs, mus, exps, and rCNV scores
gene.scores <- load.gene.scores(del.tsv, dup.tsv, rcnv.tsv)

# # Compute gene O/E summary statistics per rCNV score bin
# oe.stats <- compute.oe.stats.by.rcnv.bin(gene.scores)

# # Compute correlations for each O/E summary statistic type with rCNV score percentile
# cors <- compute.cor(oe.stats)

# # Plot O/Es of genes by their rCNV score bin
# pdf(out.pdf, height = 5, width = 10)
# plot.rcnv.bin.oe(oe.stats, cors)
# dev.off()

# Compute robust Z-scores for O/Es of genes within rCNV score bins
zs <- compute.z.by.rcnv.bin(gene.scores)

# Output genes ranked by robust Z-score
write.table(
    zs[
        order(abs(zs$DEL.bin.z), decreasing = TRUE),
        c("gene", colnames(zs)[grepl("DEL", colnames(zs))])
    ],
    out_genes_del.tsv,
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)
write.table(
    zs[
        order(abs(zs$DUP.bin.z), decreasing = TRUE),
        c("gene", colnames(zs)[grepl("DUP", colnames(zs))])
    ],
    out_genes_dup.tsv,
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)
