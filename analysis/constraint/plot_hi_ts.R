#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>

# Plot gene DEL and DUP DSMap observed/expected values (O/Es) against each other


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
load.gene.scores <- function(del.tsv, dup.cg.tsv, n.quantiles = 100) {
    # Load gene coding DEL and copy-gain DUP observed and expected counts
    del <- read.table(del.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(del)[1] <- "gene"
    dup.cg <- read.table(dup.cg.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(dup.cg)[1] <- "gene"

    # Add colname suffixes for each CNV type and merge
    colnames(del) <- ifelse(
        colnames(del) == "gene", colnames(del), paste(colnames(del), "DEL", sep = ".")
    )
    colnames(dup.cg) <- ifelse(
        colnames(dup.cg) == "gene", colnames(dup.cg), paste(colnames(dup.cg), "DUP.CG", sep = ".")
    )
    gene.scores <- merge(
        del, dup.cg,
        by = "gene",
        suffixes = c(".DEL", ".DUP.CG"), all = FALSE, sort = FALSE
    )

    # Compute expected # of gene dels & CG dups given sample size
    n.samples <- 63046
    gene.scores$exp.DEL <- (10^gene.scores$mu.DEL) * 2 * n.samples
    gene.scores$exp.DUP.CG <- (10^gene.scores$mu.DUP.CG) * 2 * n.samples

    # Compute gene O/E
    gene.scores$oe.DEL <- gene.scores$n_svs.DEL / gene.scores$exp.DEL
    gene.scores$oe.DUP.CG <- gene.scores$n_svs.DUP.CG / gene.scores$exp.DUP.CG
    gene.scores$oe.DEL.quantile <- ceiling(
        n.quantiles * rank(gene.scores$oe.DEL) / nrow(gene.scores)
    )
    gene.scores$oe.DUP.CG.quantile <- ceiling(
        n.quantiles * rank(gene.scores$oe.DUP.CG) / nrow(gene.scores)
    )
    return(gene.scores)
}

# Compute correlation coefficients and p-values for each relationship
compute.cor <- function(oes) {
    cor.estimate <- cor.test(oes$oe.DEL, oes$oe.DUP.CG, method = "spearman", exact = FALSE)
    return(list(estimate = cor.estimate$estimate, p = cor.estimate$p.value))
}


######################
# Plotting functions #
######################
# Plot DEL vs. CG DUP O/Es of genes
plot.hi.ts <- function(dat, cor) {
    par(mar = c(3.5, 3.5, 2.5, 1.5))
    plot(
        dat$oe.DEL.quantile, dat$oe.DUP.CG.quantile,
        pch = 19,
        col = adjustcolor(offblack, alpha = 0.25),
        xlab = "", ylab = "", las = 2
    )
    mtext(3, line = 1.5, text = "Gene HI vs. TS", cex = 1.2, font = 2)
    mtext(3, line = 0.25, text = paste(
        "Spearman rho =",
        round(cor[["estimate"]], 3),
        "; p =",
        formatC(cor[["p"]], format = "e", digits = 1)
    ))
    mtext(1, line = 2.5, text = paste("DEL O/E quantile"))
    mtext(2, line = 2.5, text = "CG DUP O/E quantile")
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
    make_option("--n-quantiles",
        type = "integer", default = 100,
        help = "Number of quantiles to divide genes into [default %default]"
    )
)

# Get command-line arguments & options
arg_list <- c("del.tsv", "cg_dup.tsv", "out.pdf")
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
dup.cg.tsv <- args$args[2]
out.pdf <- args$args[3]
n.quantiles <- opts$`n-quantiles`

# Load gene obs, mus, exps, and O/Es
gene.scores <- load.gene.scores(del.tsv, dup.cg.tsv, n.quantiles)
# Compute DEL O/E vs. CG DUP O/E correlation
cor.estimate <- compute.cor(gene.scores)
# Plot DEL O/E vs. CG DUP O/E
pdf(out.pdf, height = 5, width = 10)
plot.hi.ts(gene.scores, cor.estimate)
dev.off()

print(paste0(
    "# genes with 0 obs DELs: ", sum(gene.scores$n_svs.DEL == 0),
    " (", round(sum(gene.scores$n_svs.DEL == 0) / nrow(gene.scores) * 100, 1), "%)"
))
print(paste0(
    "# genes with 0 exp DELs: ", sum(gene.scores$exp.DEL == 0),
    " (", round(sum(gene.scores$exp.DEL == 0) / nrow(gene.scores) * 100, 1), "%)"
))
print(paste0(
    "# genes with 0 obs CG DUPs: ", sum(gene.scores$n_svs.DUP.CG == 0),
    " (", round(sum(gene.scores$n_svs.DUP.CG == 0) / nrow(gene.scores) * 100, 1), "%)"
))
print(paste0(
    "# genes with 0 exp CG DUPs: ", sum(gene.scores$exp.DUP.CG == 0),
    " (", round(sum(gene.scores$exp.DUP.CG == 0) / nrow(gene.scores) * 100, 1), "%)"
))
print(paste0(
    "# genes with 0 obs DELs and 0 obs CG DUPs: ",
    sum((gene.scores$n_svs.DEL == 0) & (gene.scores$n_svs.DUP.CG == 0)),
    " (",
    round(sum((gene.scores$n_svs.DEL == 0) & (gene.scores$n_svs.DUP.CG == 0)) /
        nrow(gene.scores) * 100, 1),
    "%)"
))
hist(gene.scores[gene.scores$n_svs.DEL == 0, "exp.DEL"])
hist(gene.scores[gene.scores$n_svs.DEL != 0, "exp.DEL"])
hist(gene.scores[gene.scores$n_svs.DUP.CG == 0, "exp.DUP.CG"])
hist(gene.scores[gene.scores$n_svs.DUP.CG != 0, "exp.DUP.CG"])
