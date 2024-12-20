#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot correlation of DEL and DUP mutation rate predictions in bin pairs


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
# Load mutation rate table
load.mu.tsv <- function(mu.in, na.val = -49.0) {
  # Load data
  mu <- read.table(mu.in, header = T, sep = "\t", comment.char = "", check.names = F)
  colnames(mu)[1] <- gsub("#", "", colnames(mu)[1])

  # Remove rows where mu == na.val or mu is infinite
  # (na.val is introduced by athena as a placeholder for situations where
  #  mutation rates are missing)
  # TODO: Amend these in the model
  mu <- mu[!(mu$mu == na.val) & !(is.infinite(mu$mu)), ]

  return(mu)
}

# Load bin pair table with DEL and DUP mutation rate quantiles
load.mu.scores <- function(del.tsv, dup.tsv, n.quantiles = 100) {
  # Load DEL and DUP mutation rate tables
  del.mu <- load.mu.tsv(del.tsv)
  dup.mu <- load.mu.tsv(dup.tsv)

  # Calculate quantiles
  del.mu$quantile <- ceiling(n.quantiles *
    rank(del.mu$mu, na.last = "keep", ties.method = "random") /
    nrow(del.mu))
  dup.mu$quantile <- ceiling(n.quantiles *
    rank(dup.mu$mu, na.last = "keep", ties.method = "random") /
    nrow(dup.mu))

  # Merge DEL and DUP data
  mu.scores <- merge(del.mu[, -4], dup.mu[, -4],
    by = colnames(del.mu)[1:3], suffixes = c(".DEL", ".DUP"),
    all = FALSE, sort = FALSE,
  )
  return(mu.scores)
}

# Compute bin-pair DUP mu quantile summary statistics per DEL mu quantile
compute.sum.stats.by.quantile <- function(mu.scores, n.quantiles = 100) {
  sum.stats <- cbind(1:n.quantiles, do.call("rbind", lapply(1:n.quantiles, function(i) {
    # Get bin-pairs in quantile
    idxs <- which(mu.scores$quantile.DEL == i)
    # Compute mean and CI on DUP mu quantiles
    quantile.mean <- mean(mu.scores[idxs, "quantile.DUP"])
    quantile.se <- sd(mu.scores[idxs, "quantile.DUP"]) / sqrt(length(idxs))
    data.frame(
      median(mu.scores[idxs, "quantile.DUP"]),
      quantile.mean, quantile.mean - (qnorm(0.975) * quantile.se),
      quantile.mean + (qnorm(0.975) * quantile.se)
    )
  })))
  colnames(sum.stats) <- c(
    "quantile.DEL", "quantile.DUP.median", "quantile.DUP.mean",
    "quantile.DUP.CI.lower", "quantile.DUP.CI.upper"
  )
  return(sum.stats)
}

# Compute correlation coefficients and p-values
compute.cor <- function(x, y) {
  cor.estimate <- cor.test(x, y, method = "pearson", exact = FALSE)
  return(list(estimate = cor.estimate$estimate, p = cor.estimate$p.value))
}


######################
# Plotting functions #
######################
# DEL vs. DUP mu quantile heatmap
plot.quantile.heatmap <- function(mu.scores, sum.stats, cor.in, n.quantiles = 100, 
                                  title.in = "Bin-Pair DEL vs. DUP Mu") {
  quantiles <- seq(1, n.quantiles)
  x.axis.title <- "DEL Mu Quantile"
  y.axis.title <- "DUP Mu Quantile"

  # Prep plot area
  # 1/2 buffer added to axis limits to accommodate each axis tick being
  # in the middle of a drawn row/column
  prep.plot.area(
    c(quantiles[1] - 1 / 2, quantiles[length(quantiles)] + 1 / 2),
    c(quantiles[1] - 1 / 2, quantiles[length(quantiles)] + 1 / 2),
    parmar = c(2.5, 2.5, 1.2, 1)
  )

  # For bin-pairs in each DEL mu quantile, compute density/count distribution
  # of DUP mu quantiles
  # Create matrix where rows are DEL mu quantiles and
  # columns are DUP mu quantiles (transposed when plotted with image)
  quantile.density <- table(mu.scores[, c("quantile.DEL", "quantile.DUP")])
  # Plot 0 counts as white (using mask) and nonzero counts as nonwhite
  pal <- hcl.colors(palette = "Purples 3", rev = TRUE, n = 100)
  image(
    quantiles, quantiles,
    quantile.density,
    add = TRUE,
    col = pal
  )
  quantile.density.masked <- quantile.density
  quantile.density.masked[quantile.density.masked != 0] <- NA
  image(
    quantiles, quantiles,
    quantile.density.masked,
    add = TRUE,
    col = "#FFFFFF"
  )
  # Plot medians of DUP mu quantiles for bin-pairs in each DEL quantile
  points(
    sum.stats$quantile.DEL,
    sum.stats$quantile.DUP.median,
    type = "l",
    col = adjustcolor(browns$dark2, alpha = 0.5),
    xaxt = "n", yaxt = "n", xlab = "", ylab = "",
  )
  # Plot y = x line
  abline(0, 1, col = adjustcolor(offblack, alpha = 0.5), lty = 2)

  # Add X axis
  x.ax.at <- axTicks(1)
  axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(1, at = x.ax.at, tck = -0.0125, col = offblack, labels = NA)
  axis(1, at = x.ax.at, tick = F, line = -0.8, cex.axis = 0.8)
  mtext(1, line = 1.25, text = x.axis.title)
  # Add Y axis
  y.ax.at <- axTicks(2)
  axis(2, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(2, at = x.ax.at, tck = -0.0125, col = offblack, labels = NA)
  axis(2, at = x.ax.at, tick = F, line = -0.65, cex.axis = 0.8, las = 2)
  mtext(2, line = 1.25, text = y.axis.title)

  # Add correlation stats
  title(
    main = paste(
      "Pearson rho =",
      round(cor.in$estimate, 3),
      "; p =",
      formatC(cor.in$p, format = "e", digits = 1)
    ), outer = FALSE, line = -1, cex.main = 1, font.main = 1
  )
  # Add title
  mtext(3, font = 2, text = title.in, xpd = T)
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
  make_option("--n-quantiles",
    type = "integer", default = 100,
    help = "Number of quantiles to divide genes into [default %default]"
  ),
  make_option(c("--title"),
    help = "Custom title",
    type = "character", default = NULL
  )
)

# Get command-line arguments & options
arg_list <- c("del.tsv", "dup.tsv", "out.prefix")
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
title.in <- opts$title

# Load DEL and DUP mutation rate quantiles
mu.scores <- load.mu.scores(del.tsv, dup.tsv, n.quantiles)

# Compute bin-pair DUP mu quantile summary statistics per DEL mu quantile
sum.stats <- compute.sum.stats.by.quantile(mu.scores, n.quantiles)

# Compute correlation between DEL and DUP mu quantiles over all bin pairs
del.dup.cor <- compute.cor(mu.scores$quantile.DEL, mu.scores$quantile.DUP)

# Plot DEL vs. DUP mu quantile distribution
pdf(paste(out.prefix, "DEL.DUP.mu.heatmap.pdf", sep = "."),
  height = 5, width = 5
)
plot.quantile.heatmap(
  mu.scores, sum.stats, del.dup.cor, n.quantiles,
  ifelse(is.null(title.in), "DEL vs. DUP Mu Over Bin-Pairs", title.in)
)
dev.off()
