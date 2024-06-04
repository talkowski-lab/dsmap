#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# 1. Plot gene DEL and DUP DSMap observed/expected values (O/Es) against a constraint metric
# 2. Write out genes ranked by discordance of DSMap O/Es vs. other constraint metrics


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
load.gene.scores <- function(del.tsv, dup.tsv, dup.cg.tsv, constraint.tsv, n.quantiles = 100) {
    # Load data files:
    # gene coding DEL, coding DUP, and copy-gain DUP observed and expected counts
    # gene constraint scores
    print("Loading DELs...")
    del <- read.table(del.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(del)[1] <- "gene"
    print("Loading DUPs...")
    dup <- read.table(dup.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(dup)[1] <- "gene"
    print("Loading CG DUPs...")
    dup.cg <- read.table(dup.cg.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(dup.cg)[1] <- "gene"
    print("Loading constraint scores...")
    constraint.scores <- read.table(constraint.tsv, header = TRUE, sep = "\t", comment.char = "")
    colnames(constraint.scores)[1] <- "gene"
    constraint.metric <- colnames(constraint.scores)[2]
    print(constraint.metric)

    print("Merging DEL, DUP, CG DUP data and constraint scores...")
    # Add colname suffixes for each CNV type and merge
    colnames(del) <- ifelse(
        colnames(del) == "gene", colnames(del), paste(colnames(del), "DEL", sep = ".")
    )
    colnames(dup) <- ifelse(
        colnames(dup) == "gene", colnames(dup), paste(colnames(dup), "DUP", sep = ".")
    )
    colnames(dup.cg) <- ifelse(
        colnames(dup.cg) == "gene", colnames(dup.cg), paste(colnames(dup.cg), "DUP.CG", sep = ".")
    )
    gene.scores <- Reduce(
        function(x, y) merge(x, y, by = "gene", all = FALSE, sort = FALSE),
        list(del, dup, dup.cg)
    )
    gene.scores <- merge(gene.scores, constraint.scores, by = "gene", all = FALSE, sort = FALSE)
    print(
        paste(
            "Number of DSMap genes lost in merge with constraint scores:",
            nrow(del) - nrow(gene.scores)
        )
    )
    print(
        paste(
            "Number of genes after merge with NA in constraint score:", sum(is.na(gene.scores[, constraint.metric]))
        )
    )
    gene.scores <- gene.scores[!is.na(gene.scores[, constraint.metric]), ]

    print("Calculating quantiles...")
    # Calculate quantiles of constraint score
    gene.scores$quantile <- ceiling(n.quantiles *
        rank(gene.scores[, constraint.metric], na.last = "keep", ties.method = "random") /
        nrow(gene.scores))

    print("Computing DEL, DUP, CG DUP exp and O/E...")
    # Compute expected # of gene dels & dups given sample size
    n.samples <- 63046
    gene.scores$exp.DEL <- (10^gene.scores$mu.DEL) * 2 * n.samples
    gene.scores$exp.DUP <- (10^gene.scores$mu.DUP) * 2 * n.samples
    gene.scores$exp.DUP.CG <- (10^gene.scores$mu.DUP.CG) * 2 * n.samples

    # Compute gene O/E
    gene.scores$oe.DEL <- gene.scores$n_svs.DEL / gene.scores$exp.DEL
    gene.scores$oe.DUP <- gene.scores$n_svs.DUP / gene.scores$exp.DUP
    gene.scores$oe.DUP.CG <- gene.scores$n_svs.DUP.CG / gene.scores$exp.DUP.CG
    return(list(gene.scores, constraint.metric))
}

# Compute gene O/E summary statistics with confidence intervals per constraint score quantile
compute.oe.stats.by.quantile <- function(gene.scores, n.quantiles = 100, n.bootstraps = 1000, seed = 42) {
    oe.stats <- cbind(1:n.quantiles, do.call("rbind", lapply(1:n.quantiles, function(i) {
        # Get genes in quantile
        idxs <- which(gene.scores$quantile == i)
        # Bootstrap sampling of genes in quantile
        set.seed(seed)
        bootstrap.idxs <- replicate(n.bootstraps, sample(idxs, length(idxs), replace = TRUE))

        # Compute obs and exp values with bootstrap for DEL and DUP
        do.call("cbind", lapply(c("DEL", "DUP", "DUP.CG"), function(cnv) {
            # Get obs and exp for genes in quantile and bootstrap samples
            obs <- gene.scores[idxs, paste("n_svs", cnv, sep = ".")]
            exp <- gene.scores[idxs, paste("exp", cnv, sep = ".")]
            bootstrap.obs <- matrix(
                gene.scores[bootstrap.idxs, paste("n_svs", cnv, sep = ".")],
                ncol = n.bootstraps
            )
            bootstrap.exp <- matrix(
                gene.scores[bootstrap.idxs, paste("exp", cnv, sep = ".")],
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
    })))
    colnames(oe.stats) <- c("quantile", do.call("c", lapply(c("DEL", "DUP", "DUP.CG"), function(cnv) {
        paste(c(
            "median", "median.lower", "median.upper",
            "sum", "sum.lower", "sum.upper", "mad"
        ), cnv, sep = ".")
    })))
    return(oe.stats)
}

# Compute correlation coefficients and p-values for each relationship
compute.cor <- function(oe.stats) {
    cors <- do.call("rbind", lapply(seq_along(oe.stats), function(i) {
        cor.vals <- cor.test(
            seq_along(oe.stats[, i]), unlist(oe.stats[, i]),
            method = "spearman", exact = FALSE
        )
        data.frame(names(oe.stats)[i], cor.vals$estimate, cor.vals$p.value)
    }))
    colnames(cors) <- c("sumstat.type", "estimate", "p")
    return(cors)
}

# Compute robust Z-scores of O/Es for genes within each constraint score quantile
compute.z.by.quantile <- function(gene.scores, quantile.stats) {
    zs.by.cnv <- lapply(c("DEL", "DUP", "DUP.CG"), function(cnv) {
        # Annotate genes with quantile-level median and MAD O/E
        scores <- merge(
            gene.scores[, c(
                "gene", "quantile",
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

        # Calculate gene robust Z-scores
        scores[, paste("quantile.z", cnv, sep = ".")] <- (scores[, paste("oe", cnv, sep = ".")] -
            scores[, paste("quantile.median", cnv, sep = ".")]) /
            scores[, paste("quantile.mad", cnv, sep = ".")]
        return(scores)
    })
    return(Reduce(
        function(x, y) merge(x, y, by = c("gene", "quantile"), all = FALSE, sort = FALSE),
        zs.by.cnv
    ))
}


######################
# Plotting functions #
######################
# Plot O/Es of genes by their constraint score quantile
plot.constraint.quantile.oe <- function(dat, cors, constraint.label) {
    par(mfrow = c(2, 3), mar = c(2, 3.5, 2, 1.5))

    pch_type <- list("sum" = 19, "median" = 21)
    sumstat_title <- list("sum" = "Quantile-wide Mean", "median" = "Median Gene")

    for (sumstat in c("sum", "median")) {
        for (cnv in c("DEL", "DUP", "DUP.CG")) {
            plot(
                dat[, paste(sumstat, cnv, sep = ".")],
                pch = pch_type[[sumstat]],
                col = get(paste(ifelse(cnv == "DUP.CG", "DUP", cnv), "colors", sep = "."))$main,
                xaxt = "n", xlab = "", ylab = "", las = 2,
                main = paste(cnv, " Obs/Exp ", "(", sumstat_title[[sumstat]], ")", sep = ""),
                panel.first = c(abline(h = 1, lty = 5))
            )
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
            mtext(1, text = paste("Genes Binned by", constraint.label, "quantile"))
            mtext(2, line = 2.5, text = "Observed / Expected")
        }
    }
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
    make_option("--constraint-label",
        type = "character",
        help = "Output label for constraint metric if different from input TSV column name"
    ),
    make_option("--n-quantiles",
        type = "integer", default = 100,
        help = "Number of quantiles to divide genes into [default %default]"
    )
)

# Get command-line arguments & options
arg_list <- c(
    "del.tsv", "dup.tsv", "cg_dup.tsv", "constraint.tsv", "out.pdf",
    "out_genes_del.tsv", "out_genes_dup.tsv", "out_genes_cg_dup.tsv"
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
dup.cg.tsv <- args$args[3]
constraint.tsv <- args$args[4]
out.pdf <- args$args[5]
out.genes.del.tsv <- args$args[6]
out.genes.dup.tsv <- args$args[7]
out.genes.dup.cg.tsv <- args$args[8]
n.quantiles <- opts$`n-quantiles`
constraint.label <- opts$`constraint-label`
print(constraint.label)

print("Loading scores...")
# Load gene obs, mus, exps, and constraint scores
gene.scores <- load.gene.scores(del.tsv, dup.tsv, dup.cg.tsv, constraint.tsv, n.quantiles)
constraint.metric <- gene.scores[[2]]
print(constraint.metric)
gene.scores <- gene.scores[[1]]
print(head(gene.scores))

print("Computing stats by quantile...")
# Compute gene O/E summary statistics per constraint score quantile
oe.stats <- compute.oe.stats.by.quantile(gene.scores, n.quantiles)
print(head(oe.stats))

print("Computing correlations...")
# Compute correlations for each O/E summary statistic type with constraint score percentile
cors <- compute.cor(oe.stats)

print("Plotting...")
# Plot O/Es of genes by their constraint score quantile
pdf(out.pdf, height = 5, width = 10)
plot.constraint.quantile.oe(
    oe.stats, cors,
    ifelse(is.null(constraint.label), constraint.metric, constraint.label)
)
dev.off()

print("Computing robust Z-scores...")
# Compute robust Z-scores for O/Es of genes within constraint score quantiles
zs <- compute.z.by.quantile(gene.scores, oe.stats)
print(head(zs))

print("Writing out tables...")
# Output genes ranked by robust Z-score for each CNV type
write.table(
    zs[
        order(abs(zs$quantile.z.DEL), decreasing = TRUE),
        c("gene", colnames(zs)[grepl("DEL$", colnames(zs))])
    ],
    out.genes.del.tsv,
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)
write.table(
    zs[
        order(abs(zs$quantile.z.DUP), decreasing = TRUE),
        c("gene", colnames(zs)[grepl("DUP$", colnames(zs))])
    ],
    out.genes.dup.tsv,
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)
write.table(
    zs[
        order(abs(zs$quantile.z.DUP.CG), decreasing = TRUE),
        c("gene", colnames(zs)[grepl("DUP.CG$", colnames(zs))])
    ],
    out.genes.dup.cg.tsv,
    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)
# print(paste0(
#     "# genes with 0 obs DELs: ", sum(zs$n_svs.DEL == 0),
#     " (", round(sum(zs$n_svs.DEL == 0) / nrow(zs) * 100, 1), "%)"
# ))
# print(paste0(
#     "# genes with 0 exp DELs: ", sum(zs$exp.DEL == 0),
#     " (", round(sum(zs$exp.DEL == 0) / nrow(zs) * 100, 1), "%)"
# ))
# print(paste0(
#     "# quantiles with DEL O/E median of 0: ", sum(oe.stats$median.DEL == 0),
#     " (", round(sum(oe.stats$median.DEL == 0) / nrow(oe.stats) * 100, 1), "%)"
# ))
# print(paste0(
#     "# quantiles with DEL O/E MAD of 0: ", sum(oe.stats$mad.DEL == 0),
#     " (", round(sum(oe.stats$mad.DEL == 0) / nrow(oe.stats) * 100, 1), "%)"
# ))
# print(paste0(
#     "# genes with abs(DEL O/E Z) > 3 in quantiles with median & MAD O/E > 0: ",
#     sum((abs(zs$quantile.z.DEL) > 3) & (zs$quantile.median.DEL > 0) & ((zs$quantile.mad.DEL > 0))),
#     " (", round(
#         sum((abs(zs$quantile.z.DEL) > 3) & (zs$quantile.median.DEL > 0) & ((zs$quantile.mad.DEL > 0))) /
#             sum((zs$quantile.median.DEL > 0) & ((zs$quantile.mad.DEL > 0))) * 100, 1
#     ), "%)"
# ))
# print("###########################")
# print(paste0(
#     "# genes with 0 obs DUPs: ", sum(zs$n_svs.DUP == 0),
#     " (", round(sum(zs$n_svs.DUP == 0) / nrow(zs) * 100, 1), "%)"
# ))
# print(paste0(
#     "# genes with 0 exp DUPs: ", sum(zs$exp.DUP == 0),
#     " (", round(sum(zs$exp.DUP == 0) / nrow(zs) * 100, 1), "%)"
# ))
# print(paste0(
#     "# quantiles with DUP O/E median of 0: ", sum(oe.stats$median.DUP == 0),
#     " (", round(sum(oe.stats$median.DUP == 0) / nrow(oe.stats) * 100, 1), "%)"
# ))
# print(paste0(
#     "# quantiles with DUP O/E MAD of 0: ", sum(oe.stats$mad.DUP == 0),
#     " (", round(sum(oe.stats$mad.DUP == 0) / nrow(oe.stats) * 100, 1), "%)"
# ))
# print(paste0(
#     "# genes with abs(DUP O/E Z) > 3 in quantiles with median & MAD O/E > 0: ",
#     sum((abs(zs$quantile.z.DUP) > 3) & (zs$quantile.median.DUP > 0) & ((zs$quantile.mad.DUP > 0))),
#     " (", round(
#         sum((abs(zs$quantile.z.DUP) > 3) & (zs$quantile.median.DUP > 0) & ((zs$quantile.mad.DUP > 0))) /
#             sum((zs$quantile.median.DUP > 0) & ((zs$quantile.mad.DUP > 0))) * 100, 1
#     ), "%)"
# ))
