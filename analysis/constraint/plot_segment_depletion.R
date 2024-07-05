#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distribution of DEL and DUP depletion on segments of interest


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
load.segment.scores <- function(del.tsv, dup.tsv, n.bins = 100) {
    # Load data files
    del <- read.table(del.tsv, header = TRUE, sep = "\t", comment.char = "")
    dup <- read.table(dup.tsv, header = TRUE, sep = "\t", comment.char = "")

    # # Remove rows where mu == na.val or mu is infinite
    # # (na.val is introduced by athena as a placeholder for situations where
    # #  mutation rates are missing)
    # # TODO: Amend these in the model
    # mu <- mu[!(mu$mu == na.val) & !(is.infinite(mu$mu)), ]

    colnames(del)[1] <- "segment"
    colnames(dup)[1] <- "segment"

    # Merge DEL and DUP data
    segment.scores <- merge(del, dup,
        by = colnames(del)[1],
        suffixes = c(".DEL", ".DUP"), all = FALSE, sort = FALSE
    )

    # Compute expected # of dels and dups given sample size
    n.samples <- 63046
    segment.scores$exp.DEL <- (10^segment.scores$mu.DEL) * 2 * n.samples
    segment.scores$exp.DUP <- (10^segment.scores$mu.DUP) * 2 * n.samples

    # Compute segment O/E
    segment.scores$oe.DEL <- segment.scores$n_svs.DEL / segment.scores$exp.DEL
    segment.scores$oe.DUP <- segment.scores$n_svs.DUP / segment.scores$exp.DUP
    segment.scores$oe.DEL.quantile <- ceiling(
        n.quantiles * rank(segment.scores$oe.DEL) / nrow(segment.scores)
    )
    segment.scores$oe.DUP.quantile <- ceiling(
        n.quantiles * rank(segment.scores$oe.DUP) / nrow(segment.scores)
    )
    return(segment.scores)
}

# Compute correlation coefficients and p-values for each relationship
compute.cor <- function(oes) {
    cor.estimate <- cor.test(oes$oe.DEL, oes$oe.DUP, method = "spearman", exact = FALSE)
    return(list(estimate = cor.estimate$estimate, p = cor.estimate$p.value))
}


######################
# Plotting functions #
######################
# Plot segment CNV O/E distribution
plot.cnv.oe <- function(dat, cnv, segment.name = "Segment") {
    # Set plot parameters
    if (cnv %in% c("DEL", "DUP", "CNV")) {
        bar.color <- get(paste(cnv, "colors", sep = "."))$main
    } else {
        bar.color <- browns$main
    }

    oes <- dat[, paste("oe", cnv, sep = ".")]
    # Replace value for OE = 0 before taking log to avoid infinite log
    oe.zero.val <- 10**-3
    oes[oes == 0] <- oe.zero.val
    # oes <- oes[oes != 0]
    # Take log of O/Es
    oes <- log10(oes)
    # Set limits for x-axis plotting
    # x.upper.lim <- max(log10(logscale.major))
    # x.lower.lim <- floor(min(oes, na.rm = TRUE))
    x.lower.lim <- -3
    x.upper.lim <- 5
    xlims <- c(x.lower.lim, x.upper.lim)

    # Prep plot area
    prep.plot.area(xlims, c(0, 1), parmar = c(2.5, 2.8, 1.2, 1))

    # Add histogram
    h <- hist(oes, plot = FALSE, breaks = 200)
    rect(
        xleft = h$breaks[-length(h$breaks)], xright = h$breaks[-1],
        ybottom = 0, ytop = h$counts / length(oes), col = bar.color, border = "white"
    )

    # Add X axis
    axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
    x.ax.at <- log10(logscale.major)
    axis(1, at = x.ax.at, tck = -0.025, col = offblack, labels = NA)
    sapply(x.ax.at, function(x) {
        axis(1,
            at = x, tick = F, line = -0.65, labels = bquote(10^.(x)),
            cex.axis = 0.85
        )
    })
    mtext(1, line = 1.25, text = paste(cnv, "O/E"))

    # Add Y axis
    mtext(2, line = 1.9, text = paste("Proportion of ", segment.name, "s", sep = ""))
    axis(2, tck = -0.025, col = offblack, labels = NA)
    axis(2, tick = F, line = -0.3, cex.axis = 0.85, las = 1)

    # Add title
    mtext(3, font = 2, text = paste(segment.name, cnv, "O/E"), xpd = T)
}

# Plot segment DEL vs. DUP O/E quantiles
plot.del.dup.quantile <- function(dat, cor, segment.name) {
    par(mar = c(3.5, 3.5, 2.5, 1.5))
    plot(
        dat$oe.DEL.quantile, dat$oe.DUP.quantile,
        pch = 19,
        col = adjustcolor(offblack, alpha = 0.25),
        xlab = "", ylab = "", las = 2
    )
    mtext(3, line = 1.5, text = paste(segment.name, "DEL vs. DUP O/E"), cex = 1.2, font = 2)
    mtext(3, line = 0.25, text = paste(
        "Spearman rho =",
        round(cor[["estimate"]], 3),
        "; p =",
        formatC(cor[["p"]], format = "e", digits = 1)
    ))
    mtext(1, line = 2.5, text = paste("DEL O/E quantile"))
    mtext(2, line = 2.5, text = "DUP O/E quantile")
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
    make_option("--n-quantiles",
        type = "integer", default = 100,
        help = "Number of quantiles to divide segments into [default %default]"
    ),
    make_option("--segment-name",
        type = "character", default = "Segment",
        help = "Segment type (for plot titles only) [default %default]"
    )
)

# Get command-line arguments & options
arg_list <- c("del.tsv", "dup.tsv", "out_prefix")
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
out.prefix <- args$args[3]
n.quantiles <- opts$`n-quantiles`
segment.name <- opts$`segment-name`

# Load segment obs, mus, exps, and O/Es
segment.scores <- load.segment.scores(del.tsv, dup.tsv, n.quantiles)
# Compute DEL O/E vs. DUP O/E correlation
cor.estimate <- compute.cor(segment.scores)

# Plot DEL O/E distribution
for (cnv in c("DEL", "DUP")) {
    pdf(paste(out.prefix, cnv, "oe.pdf", sep = "."), height = 5, width = 5)
    plot.cnv.oe(segment.scores, cnv, segment.name)
    dev.off()
}

# Plot DEL O/E vs. DUP O/E quantiles
pdf(paste(out.prefix, "DEL_DUP_quantiles.pdf", sep = "."), height = 5, width = 5)
plot.del.dup.quantile(segment.scores, cor.estimate, segment.name)
dev.off()

print(paste0(
    "# segments with 0 obs DELs: ", sum(segment.scores$n_svs.DEL == 0),
    " (", round(sum(segment.scores$n_svs.DEL == 0) / nrow(segment.scores) * 100, 1), "%)"
))
print(paste0(
    "# segments with 0 exp DELs: ", sum(segment.scores$exp.DEL == 0),
    " (", round(sum(segment.scores$exp.DEL == 0) / nrow(segment.scores) * 100, 1), "%)"
))
print(paste0(
    "# segments with 0 obs DUPs: ", sum(segment.scores$n_svs.DUP == 0),
    " (", round(sum(segment.scores$n_svs.DUP == 0) / nrow(segment.scores) * 100, 1), "%)"
))
print(paste0(
    "# segments with 0 exp DUPs: ", sum(segment.scores$exp.DUP == 0),
    " (", round(sum(segment.scores$exp.DUP == 0) / nrow(segment.scores) * 100, 1), "%)"
))
print(paste0(
    "# segments with 0 obs DELs and 0 obs DUPs: ",
    sum((segment.scores$n_svs.DEL == 0) & (segment.scores$n_svs.DUP == 0)),
    " (",
    round(sum((segment.scores$n_svs.DEL == 0) & (segment.scores$n_svs.DUP == 0)) /
        nrow(segment.scores) * 100, 1),
    "%)"
))
