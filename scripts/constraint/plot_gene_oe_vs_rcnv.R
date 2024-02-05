#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot gene DEL and DUP observed/expected values (O/Es) by rCNV pHaplo or pTriplo bin


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly = TRUE)
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
    return(gene.scores)
}

# Compute gene O/E summary statistics per rCNV score bin
compute.oe.by.rcnv.bin <- function(gene.scores) {
    oe <- as.data.frame(do.call("cbind", lapply(c("DEL", "DUP"), function(cnv) {
        t(sapply(1:100, function(i) {
            idxs <- gene.scores[, paste(cnv, "pct", sep = ".")] == i
            obs <- gene.scores[idxs, paste("n_svs", cnv, sep = ".")]
            exp <- gene.scores[idxs, paste("exp", cnv, sep = ".")]
            data.frame(median(obs / exp), sum(obs) / sum(exp))
        }))
    })))
    colnames(oe) <- c("median.DEL", "sum.DEL", "median.DUP", "sum.DUP")
    return(oe)
}

# Compute correlation coefficients and p-values for each relationship
compute.cor <- function(oe) {
    cors <- as.data.frame(do.call("rbind", lapply(seq_along(oe), function(i) {
        cor.vals <- cor.test(seq_along(oe[, i]), unlist(oe[, i]), method = "spearman")
        data.frame(names(oe)[i], cor.vals$estimate, cor.vals$p.value)
    })))
    colnames(cors) <- c("sum.type", "estimate", "p")
    return(cors)
}


######################
# Plotting functions #
######################
# Plot O/Es of genes by their rCNV score bin
plot.rcnv.bin.oe <- function(dat, cors) {
    par(mfrow = c(2, 2), mar = c(2, 3.5, 2, 1.5))
    plot(as.numeric(dat$sum.DEL),
        pch = 19, col = DEL.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DEL Obs/Exp (Bin-wide Mean)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    title(
        main = paste(
            "Spearman rho =",
            round(cors[cors$sum.type == "sum.DEL", "estimate"], 3),
            "; p =",
            formatC(cors[cors$sum.type == "median.DUP", "p"], format = "e", digits = 1)
        ),
        outer = FALSE, line = -1, cex.main = 1, font.main = 1
    )
    mtext(1, text = "Genes Binned by pHaplo Percentile")
    mtext(2, line = 2.5, text = "Observed / Expected")
    plot(as.numeric(dat$sum.DUP),
        pch = 19, col = DUP.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DUP Obs/Exp (Bin-wide Mean)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    title(
        main = paste(
            "Spearman rho =",
            round(cors[cors$sum.type == "sum.DUP", "estimate"], 3),
            "; p =",
            formatC(cors[cors$sum.type == "sum.DUP", "p"], format = "e", digits = 1)
        ),
        outer = FALSE, line = -1, cex.main = 1, font.main = 1
    )
    mtext(1, text = "Genes Binned by pTriplo Percentile")
    mtext(2, line = 2.5, text = "Observed / Expected")
    plot(as.numeric(dat$median.DEL),
        pch = 21, col = DEL.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DEL Obs/Exp (Median Gene)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    title(
        main = paste(
            "Spearman rho =",
            round(cors[cors$sum.type == "median.DEL", "estimate"], 3),
            "; p =",
            formatC(cors[cors$sum.type == "median.DEL", "p"], format = "e", digits = 1)
        ),
        outer = FALSE, line = -1, cex.main = 1, font.main = 1
    )
    mtext(1, text = "Genes Binned by pHaplo Percentile")
    mtext(2, line = 2.5, text = "Observed / Expected")
    plot(as.numeric(dat$median.DUP),
        pch = 21, col = DUP.colors$main,
        xaxt = "n", xlab = "", ylab = "", las = 2, main = "DUP Obs/Exp (Median Gene)",
        panel.first = c(abline(h = 1, lty = 5))
    )
    title(
        main = paste(
            "Spearman rho =",
            round(cors[cors$sum.type == "median.DUP", "estimate"], 3),
            "; p =",
            formatC(cors[cors$sum.type == "median.DUP", "p"], format = "e", digits = 1)
        ),
        outer = FALSE, line = -1, cex.main = 1, font.main = 1
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
        usage = "%prog del.tsv dup.tsv rcnv.tsv out.pdf",
        option_list = option_list
    ),
    positional_arguments = TRUE
)
opts <- args$options

# Checks for appropriate positional arguments
if (length(args$args) != 4) {
    stop(
        paste(
            "Four positional arguments required: del.tsv, dup.tsv, rcnv.tsv, and out.pdf\n",
            sep = " "
        )
    )
}

# Writes args & opts to vars
del.tsv <- args$args[1]
dup.tsv <- args$args[2]
rcnv.tsv <- args$args[3]
out.pdf <- args$args[4]

# Load gene scores
gene.scores <- load.gene.scores(del.tsv, dup.tsv, rcnv.tsv)

# Compute gene O/E summary statistics per rCNV bin
oe <- compute.oe.by.rcnv.bin(gene.scores)

# Compute correlations for each O/E summary statistic type with rCNV score percentile
cors <- compute.cor(oe)

# Plot O/Es of genes by their rCNV score bin
pdf(out.pdf, height = 5, width = 10)
plot.rcnv.bin.oe(oe, cors)
dev.off()
