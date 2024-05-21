#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distribution of mutation rate estimates


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

  # Assign NA to rows where mu == drop.val
  # (These are introduced by athena as a placeholder for situations where
  #  mutation rates are missing)
  mu[which(mu$mu == na.val), "mu"] <- NA

  return(mu)
}


######################
# Plotting functions #
######################
# Histogram of mutation rates
mu.hist <- function(mu, cnv = NULL, x.axis.title = NULL, y.axis.title = "Loci",
                    title = "Mutation rate") {
  # Set plot parameters
  if (cnv %in% c("DEL", "DUP", "CNV")) {
    bar.color <- get(paste(cnv, "colors", sep = "."))$main
  } else {
    bar.color <- browns$main
  }
  if (is.null(x.axis.title)) {
    x.axis.title <- paste(cnv, "mutation rate")
  }
  xlims <- c(min(mu$mu), max(mu$mu))
  h <- hist(mu$mu, plot = F, breaks = 100)

  # Prep plot area
  prep.plot.area(xlims, range(h$counts), parmar = c(2.5, 2.8, 1.2, 1))

  # Add bars
  rect(
    xleft = h$breaks[-length(h$breaks)], xright = h$breaks[-1],
    ybottom = 0, ytop = h$counts, col = bar.color, border = "white"
  )

  # Add X axis
  if (min(mu$mu) > min(log10(logscale.major))) {
    x.ax.at <- log10(logscale.major)
    axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
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
  }
  y.ax.labels <- prettyNum(y.ax.at / denom, big.mark = ",")
  axis(2, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(2, at = y.ax.at, tck = -0.025, col = offblack, labels = NA)
  axis(2, at = y.ax.at, tick = F, line = -0.65, labels = y.ax.labels, las = 2)
  mtext(2, line = 1.9, text = paste(y.axis.title, " (", units, ")", sep = ""))

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
  make_option(c("--x-title"),
    help = "Custom X-axis title",
    type = "character", default = NULL
  ),
  make_option(c("--y-title"),
    help = "Custom Y-axis title",
    type = "character", default = NULL
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
x.title <- opts$`x-title`
y.title <- opts$`y-title`

# Load mutation rates
mu <- load.mu.tsv(mu.in)

# Plot histogram of mutation rates
pdf(paste(out.prefix, "mutation_rate_hist.pdf", sep = "."),
  height = 3, width = 4.5
)
mu.hist(mu, cnv, x.title, y.title)
dev.off()
