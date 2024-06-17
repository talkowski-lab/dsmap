#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot histogram of mutation rate predictions over targets


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


######################
# Plotting functions #
######################
# Histogram of mutation rates
mu.hist <- function(mu, cnv = NULL, x.axis.title = "CNVs per allele per generation",
                    y.axis.title = "Loci", title = "Mutation rate") {
  # Set plot parameters
  if (cnv %in% c("DEL", "DUP", "CNV")) {
    bar.color <- get(paste(cnv, "colors", sep = "."))$main
  } else {
    bar.color <- browns$main
  }
  xlims <- c(floor(min(mu$mu)), ceiling(max(mu$mu)))
  h <- hist(mu$mu, plot = F, breaks = 100)

  # Prep plot area
  prep.plot.area(xlims, range(h$counts), parmar = c(2.5, 2.8, 1.2, 1))

  # Add bars
  rect(
    xleft = h$breaks[-length(h$breaks)], xright = h$breaks[-1],
    ybottom = 0, ytop = h$counts, col = bar.color, border = "white"
  )

  # Add X axis
  axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
  if (min(mu$mu) >= min(log10(logscale.major))) {
    x.ax.at <- log10(logscale.major)
    axis(1, at = log10(logscale.minor), tck = -0.0125, col = offblack, labels = NA)
  } else {
    x.ax.at <- floor(min(mu$mu)):log10(max(logscale.major))
    x.ax.at.demi <- as.numeric(sapply(10^(x.ax.at), function(e) {
      log10(c(1, 5) * e)
    }))
    axis(1, at = x.ax.at.demi, tck = -0.0125, col = offblack, labels = NA)
  }
  axis(1, at = x.ax.at, tck = -0.025, col = offblack, labels = NA)
  sapply(x.ax.at, function(x) {
    axis(1,
      at = x, tick = F, line = -0.65, labels = bquote(10^.(x)),
      cex.axis = 0.85
    )
  })
  mtext(1, line = 1.25, text = x.axis.title)

  # Add Y axis
  y.ax.at <- axTicks(2)
  if (max(y.ax.at) > 1000000) {
    denom <- 1000000
    units <- "millions"
  } else if (max(y.ax.at) > 1000) {
    denom <- 1000
    units <- "thousands"
  } else {
    denom <- 1
    units <- NULL
  }
  y.ax.labels <- prettyNum(y.ax.at / denom, big.mark = ",")
  axis(2, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(2, at = y.ax.at, tck = -0.025, col = offblack, labels = NA)
  axis(2, at = y.ax.at, tick = F, line = -0.65, labels = y.ax.labels, las = 2)
  if (!is.null(units)) {
    y.axis.title <- paste(y.axis.title, " (", units, ")", sep = "")
  }
  mtext(2, line = 1.9, text = y.axis.title)

  # Add title
  mtext(3, font = 2, text = title, xpd = T)
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
  make_option(c("--cnv"),
    help = "Specify CNV type. Used for plotting colors only.",
    type = "character", default = NA
  ),
  make_option(c("--title"),
    help = "Custom title [default '%default']",
    type = "character", default = "Mutation rate"
  ),
  make_option(c("--x-title"),
    help = "Custom X-axis title [default '%default']",
    type = "character", default = "CNVs per allele per generation"
  ),
  make_option(c("--y-title"),
    help = "Custom Y-axis title [default '%default']",
    type = "character", default = "Loci"
  )
)

# Get command-line arguments & options
args <- parse_args(
  OptionParser(
    usage = "%prog mu.tsv out_prefix",
    option_list = option_list
  ),
  positional_arguments = TRUE
)
opts <- args$options

# Checks for appropriate positional arguments
if (length(args$args) != 2) {
  stop(paste("Two positional arguments required: mu.tsv and output_prefix\n", sep = " "))
}

# Writes args & opts to vars
mu.in <- args$args[1]
out.prefix <- args$args[2]
cnv <- opts$cnv
title <- opts$title
x.title <- opts$`x-title`
y.title <- opts$`y-title`

# Load mutation rates
mu <- load.mu.tsv(mu.in)

# Compute bin size
binsize <- infer.bin.size(mu)

# Plot histogram of mutation rates
pdf(paste(out.prefix, "mutation_rate.hist.pdf", sep = "."),
  height = 3, width = 4.5
)
mu.hist(mu,
  cnv = cnv, x.axis.title = x.title,
  y.axis.title = y.title, title = title
)
dev.off()

