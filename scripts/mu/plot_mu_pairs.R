#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distribution of mutation rate predictions in bin pairs


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly = TRUE)
require(optparse, quietly = TRUE)
require(viridis, quietly = TRUE)

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
# Mutation rate by bin pair distance
mu.distance <- function(mu, binsize, cnv = NULL,
                        title = "Mutation rate density by pair distance") {
  # Set plot parameters
  if (cnv == "DEL") {
    pal <- "Reds"
  } else if (cnv == "DUP") {
    pal <- "Blues"
  } else {
    pal <- "Oranges"
  }
  x.axis.title <- "Pair distance"
  y.axis.title <- paste(ifelse(is.null(cnv), "CNV", cnv),
    "s per allele per generation",
    sep = ""
  )

  # Compute bin pair sizes
  mu$size <- mu$end - mu$start - binsize
  sizes <- seq(0, max(mu$size), binsize)

  # Prep plot area
  ylims <- c(floor(min(mu$mu)), ceiling(max(mu$mu)))
  prep.plot.area(
    c(min(mu$size) - binsize / 2, max(mu$size) + binsize / 2),
    ylims,
    parmar = c(2.5, 2.8, 1.2, 1)
  )

  # For each bin pair size, compute density distribution of mu in log space
  # Create matrix where rows are bin pair sizes and
  # columns are densities in mu bins (transposed when plotted with image)
  ybreaks.by <- 0.5
  ybreaks <- seq(ylims[1], ylims[2], ybreaks.by)
  mu.size.density <- do.call("rbind", lapply(sizes, function(s) {
    hist(mu[mu$size == s, "mu"], ybreaks, plot = F)$density
  }))
  # Plot at midpoints of histogram bars
  image(
    sizes, ybreaks[-length(ybreaks)] + ybreaks.by / 2,
    mu.size.density,
    add = TRUE,
    col = hcl.colors(palette = pal, rev = TRUE, n = 50)
  )

  # Add X axis
  x.ax.at <- axTicks(1)
  if (max(x.ax.at) > 1000000) {
    denom <- 1000000
    units <- "Mb"
  } else if (max(x.ax.at) > 1000) {
    denom <- 1000
    units <- "kb"
  } else {
    denom <- 1
    units <- NULL
  }
  x.ax.labels <- prettyNum(x.ax.at / denom, big.mark = ",")
  axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(1, at = x.ax.at, tck = -0.025, col = offblack, labels = NA)
  axis(1, at = x.ax.at, tick = F, line = -0.65, labels = x.ax.labels)
  if (!is.null(units)) {
    x.axis.title <- paste(x.axis.title, " (", units, ")", sep = "")
  }
  mtext(1, line = 1.25, text = x.axis.title)

  # Add Y axis
  axis(2, at = c(-10e10, 10e10), col = offblack, tck = 0)
  if (min(mu$mu) >= min(log10(logscale.major))) {
    y.ax.at <- log10(logscale.major)
    axis(2, at = log10(logscale.minor), tck = -0.0125, col = offblack, labels = NA)
  } else {
    y.ax.at <- floor(min(mu$mu)):log10(max(logscale.major))
  }
  axis(2, at = y.ax.at, tck = -0.025, col = offblack, labels = NA)
  sapply(y.ax.at, function(y) {
    axis(2,
      at = y, tick = F, line = -0.65, labels = bquote(10^.(y)),
      cex.axis = 0.5, las = 2
    )
  })
  mtext(2, line = 1.25, text = y.axis.title)

  # Add title
  mtext(3, font = 2, text = title, xpd = T)
}

# Diagonal heatmap of mutation rate in each bin pair along chromosome
mu.heatmap <- function(mu, binsize, title = "Mutation rate") {
  # Temporary stand-in fix for dsmapR rotate.points
  # TODO: Fix in dsmapR
  rotate.points <- function(coords, angle, x.origin = 0, y.origin = 0) {
    # Process each point
    res <- as.data.frame(t(apply(coords, 1, function(xy) {
      # Assign coordinates to variables
      x <- as.numeric(xy[1]); y <- as.numeric(xy[2])

      # Convert new angle degrees to radians
      new.angle <- angle * (pi / 180)

      # Infer current angle from (x - x.o, y - y.o)
      old.angle <- atan2(y - y.origin, x - x.origin)

      # Theta = old.angle + new.angle
      theta <- old.angle + new.angle

      # Compute distance from point to origin
      r <- sqrt(((x - x.origin)^2) + ((y - y.origin)^2))

      # New x = (r * cos(theta)) + x.origin
      # New y = (r * sin(theta)) + y.origin
      c(
        r * cos(theta) + x.origin,
        r * sin(theta) + y.origin
      )
    })))

    colnames(res) <- c("x", "y")

    return(res)
  }

  # Set plot parameters
  x.axis.title <- "Coordinates"
  legend.title <- "log10(mu)"

  # Set plot limits (in rotated coordinates)
  # Note: X axis starts at first bin pair starting coordinate, Y axis starts at 0
  x.min <- min(mu$start)
  raw.x.max <- max(mu$end)
  x.max <- rotate.points(matrix(c(raw.x.max, raw.x.max - x.min), nrow = 1), angle = -45, x.origin = x.min)$x
  diag.coords <- mu[min(which((mu$end - mu$start) == max(mu$end - mu$start))), c("start", "end")]
  y.max <- rotate.points(matrix(c(diag.coords[1], diag.coords[2] - x.min), nrow = 1), angle = -45, x.origin = x.min)$y

  # Choose color palette
  n.colors <- 101
  heat.pal <- viridis(n.colors)
  # Rescale mu values to match palette
  val.range <- range(mu[!is.na(mu$mu) & !is.infinite(mu$mu), ]$mu)
  mu.vals <- ceiling((mu$mu - val.range[1]) * (100 / diff(val.range))) + 1

  # Create one panel for legend and one for heatmap
  layout(t(1:2), widths = c(0.02, 1))

  # Prep plot area for legend
  prep.plot.area(xlims = c(0, 1), ylims = val.range, parmar = c(1.3, 1.8, 1.3, 0.5))

  # Add legend
  mu.legend <- seq(val.range[1], val.range[2], length.out = n.colors + 1)
  mu.legend.at <- pretty(val.range, 8)
  rect(0, mu.legend[-(n.colors + 1)], 1, mu.legend[-1], col = heat.pal, border = NA)
  axis(2, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(2, at = mu.legend.at, tck = -0.1, col = offblack, labels = NA)
  axis(2, at = mu.legend.at, tick = F, line = -0.65, labels = mu.legend.at, cex.axis = 0.9, las = 2)
  mtext(3, line = 0.25, text = legend.title, cex = 0.9)

  # Prep plot area for heatmap
  prep.plot.area(xlims = c(x.min, x.max), ylims = c(0, y.max), parmar = c(2.5, 1, 1.5, 1))

  # Add rectangles for each bin pair
  rect.x <- unlist(lapply(mu$start, function(s) c(s, s + binsize, s + binsize, s, NA)))
  rect.y <- unlist(lapply(mu$end - x.min, function(e) c(e - binsize, e - binsize, e, e, NA)))
  rotated.coords <- rotate.points(data.frame("x" = rect.x, "y" = rect.y), angle = -45, x.origin = x.min)
  polygon(
    x = rotated.coords$x, y = rotated.coords$y,
    border = heat.pal[mu.vals],
    col = heat.pal[mu.vals],
    lwd = 0.25
  )

  # Add X axis
  x.ax.at <- axTicks(1)
  if (max(x.ax.at) > 1000000) {
    denom <- 1000000
    units <- "Mb"
  } else if (max(x.ax.at) > 1000) {
    denom <- 1000
    units <- "kb"
  } else {
    denom <- 1
    units <- NULL
  }
  x.ax.labels <- prettyNum(x.ax.at / denom, big.mark = ",")
  axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(1, at = x.ax.at, tck = -0.025, col = offblack, labels = NA)
  axis(1, at = x.ax.at, tick = F, line = -0.65, labels = x.ax.labels)
  if (!is.null(units)) {
    x.axis.title <- paste(x.axis.title, " (", units, ")", sep = "")
  }
  mtext(1, line = 1.25, text = x.axis.title)

  # Add title
  mtext(3, line = 0.25, font = 2, text = title, xpd = T)
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
  make_option(c("--distance"),
    help = paste(
      "Plot mutation rate distribution grouped by bin pair distance.",
      "Otherwise, plot heatmap along chromosome."
    ),
    action = "store_true", default = FALSE
  ),
  make_option(c("--cnv"),
    help = "Specify CNV type. Used for plotting colors and titles only.",
    type = "character", default = NULL
  ),
  make_option(c("--title"),
    help = "Custom title [default 'Mutation rate' or 'Mutation rate density by pair distance']",
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
distance <- opts$distance
cnv <- opts$cnv
title <- opts$title

# Load mutation rates
mu <- load.mu.tsv(mu.in)

# Compute bin size
binsize <- infer.bin.size(mu)

if (distance) {
  # Plot mutation rate by bin pair distance
  pdf(paste(out.prefix, "mutation_rate.size.pdf", sep = "."),
    height = 3, width = 4.5
  )
  mu.distance(mu,
    binsize = binsize, cnv = cnv,
    title = ifelse(is.null(title), "Mutation rate density by pair distance", title)
  )
  dev.off()
} else {
  # Plot mutation rate heatmap along chromosome
  pdf(paste(out.prefix, "mutation_rate.heatmap.pdf", sep = "."),
    height = 3, width = 40
  )
  mu.heatmap(mu,
    binsize = binsize,
    title = ifelse(is.null(title), "Mutation rate", title)
  )
  dev.off()
}
