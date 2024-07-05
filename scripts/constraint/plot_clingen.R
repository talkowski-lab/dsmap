#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot gene DEL and CG DUP O/Es by ClinGen dosage sensitivity evidence tier


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly = TRUE)
require(optparse, quietly = TRUE)

# Set global options and constants
options(stringsAsFactors = FALSE, scipen = 1000)
dsmapR::load.constants(c("colors", "scales"))


##################
# Data functions #
##################
# Load data
load.gene.scores <- function(cnv.tsv, clingen.txt) {
    # Load gene CNV observed and expected counts
    cnv.scores <- read.table(cnv.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(cnv.scores)[1] <- "gene"

    # Load and preprocess ClinGen gene dosage sensitivity data
    clingen <- read.table(clingen.txt, header = TRUE, sep = ",", comment.char = "")
    clingen <- clingen[, c(2, 4, 5)]
    colnames(clingen) <- c("gene", "HI_evidence", "TS_evidence")
    # Remove HGNC ID info from gene name
    clingen$gene <- gsub("HGNC.*", "", clingen$gene)

    # Merge data
    gene.scores <- merge(cnv.scores, clingen, by = "gene", all.x = TRUE, sort = FALSE)
    # Relabel genes with no ClinGen data to "NA" for plotting
    for (dosage.sens in c("HI", "TS")) {
        gene.scores[, paste(dosage.sens, "evidence", sep = "_")] = ifelse(
            (gene.scores[, paste(dosage.sens, "evidence", sep = "_")] %in% "") |
                (is.na(gene.scores[, paste(dosage.sens, "evidence", sep = "_")])),
            "NA",
            gene.scores[, paste(dosage.sens, "evidence", sep = "_")]
        )
        gene.scores[, paste(dosage.sens, "evidence", sep = "_")] = factor(
            gene.scores[, paste(dosage.sens, "evidence", sep = "_")],
            levels = c(
                "SufficientEvidence", "EmergingEvidence", "LittleEvidence",
                "AutosomalRecessive", "NoEvidence", "SensitivityUnlikely", "NA"
            )
        )
    }

    # Compute expected # of gene dels & dups given sample size
    n.samples <- 63046
    gene.scores$exp <- (10^gene.scores$mu) * 2 * n.samples

    # Compute gene O/E
    gene.scores$oe <- gene.scores$n_svs / gene.scores$exp
    return(gene.scores)
}


######################
# Plotting functions #
######################
# Plot CNV O/E distribution in ClinGen dosage sensitive genes
# by evidence tier vs. unannotated genes
plot.cnv.oe <- function(dat, cnv, dosage.sens) {
    # Set plot parameters
    if (cnv %in% c("DEL", "DUP", "CNV")) {
        bar.color <- get(paste(cnv, "colors", sep = "."))$main
    } else {
        bar.color <- browns$main
    }

    # Set value to replace OE = 0 before taking log to avoid infinite log
    # oe.zero.val <- 10**(floor(log10(min(dat[dat$oe != 0, ]$oe))) - 1)
    oe.zero.val <- 10**-3

    # Set upper and lower x-limits for plots
    # x.upper.lim <- max(log10(logscale.major))
    # x.lower.lim <- floor(min(oes, na.rm = TRUE))
    x.lower.lim <- log10(oe.zero.val)
    x.upper.lim <- ceiling(max(log10(dat$oe)))
    x.ax.at <- x.lower.lim:x.upper.lim

    # Extract ClinGen evidence tiers for dosage sensitivity type
    tiers <- sort(as.numeric(unique(dat[, paste(dosage.sens, "evidence", sep = "_")])))

    # Generate histogram distributions per evidence tier
    h.breaks.by <- 0.2
    h.breaks <- seq(x.lower.lim, x.upper.lim, h.breaks.by)
    hs <- lapply(tiers, function(t) {
        # Extract tier gene O/Es
        oes <- dat[dat[, paste(dosage.sens, "evidence", sep = "_")] %in%
            levels(dat[, paste(dosage.sens, "evidence", sep = "_")])[t], ]$oe
        # Replace value for OE = 0 before taking log to avoid infinite log
        oes[oes == 0] <- oe.zero.val
        # Take log of O/Es
        oes <- log10(oes)
        # Compute histogram distribution
        return(hist(oes, plot = FALSE, breaks = h.breaks))
    })
    max.h.prop <- max(do.call("c", lapply(hs, function(h) h$counts / sum(h$counts))))

    # Set plot dimensions
    par(mfrow = c(length(tiers), 1), mar = c(2, 3.5, 2, 1.5))

    # Plot panels for each evidence tier
    for (i in seq_along(tiers)) {
        # Set x and y limits for panel
        prep.plot.area(c(x.lower.lim, x.upper.lim), c(0, max.h.prop), parmar = c(2.3, 3, 1.2, 1))

        # Add histogram
        rect(
            xleft = hs[[i]]$breaks[-length(hs[[i]]$breaks)], xright = hs[[i]]$breaks[-1],
            ybottom = 0, ytop = hs[[i]]$counts / sum(hs[[i]]$counts), col = bar.color, border = "white"
        )

        # Add X axis
        axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
        axis(1, at = x.ax.at, tck = -0.025, col = offblack, labels = NA)
        sapply(x.ax.at, function(x) {
            axis(1,
                at = x, tick = F, line = -0.65, labels = bquote(10^.(x)),
                cex.axis = 0.85
            )
        })

        # Add Y axis
        mtext(2, line = 1.9, text = "Proportion of genes")
        axis(2, tck = -0.025, col = offblack, labels = NA)
        axis(2, tick = F, line = -0.3, cex.axis = 0.85, las = 1)

        # Add title
        mtext(3,
            font = 2,
            text = paste(
                dosage.sens, "=",
                levels(dat[, paste(dosage.sens, "evidence", sep = "_")])[tiers[i]],
                paste("(n=", sum(hs[[i]]$counts), ")", sep = "")
            ), xpd = T
        )
    }
    # Label last X axis
    mtext(1, line = 1.25, text = paste(cnv, "O/E"))
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
    make_option(c("--cnv"),
        help = "Specify CNV type. Used for plotting colors only.",
        type = "character", default = NA
    )
)

# Get command-line arguments & options
arg_list <- c("cnv.tsv", "clingen.txt", "out.prefix")
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
cnv.tsv <- args$args[1]
clingen.txt <- args$args[2]
out.prefix <- args$args[3]
cnv <- opts$cnv

# Load gene obs, mus, exps, and constraint scores
gene.scores <- load.gene.scores(cnv.tsv, clingen.txt)

# Plot O/Es of genes by their ClinGen dosage sensitivity evidence tier
for (dosage.sens in c("HI", "TS")) {
    pdf(paste(out.prefix, cnv, "clingen", dosage.sens, "pdf", sep = "."), height = 13, width = 10)
    plot.cnv.oe(gene.scores, cnv, dosage.sens)
    dev.off()
}
