#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>

# Plot correlation of DEL and DUP model weights


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
# Load DEL and DUP model weights on PCs into table
load.weights <- function(del.tsv, dup.tsv) {
  del.weights <- read.table(del.tsv)[, 1]
  dup.weights <- read.table(dup.tsv)[, 1]
  if (length(del.weights) != length(dup.weights)) {
    stop(
      paste(
        "Length of DEL model weights (", length(del.weights), ") and DUP model weights (",
        length(dup.weights), ") do not match",
        sep = ""
      )
    )
  }
  return(data.frame(
    pc = seq_along(del.weights), weight.DEL = del.weights, weight.DUP = dup.weights
  ))
}

# Identify PCs with outlier differences between DEL and DUP model weights
annot.outliers <- function(weights.in) {
  # Define outliers as having distance to y=x >2 SDs away from the mean across PCs
  weights.in$dist <- abs(weights.in$weight.DEL - weights.in$weight.DUP) / sqrt(2)
  weights.in$outlier <- (
    weights.in$dist >= (2 * sd(weights.in$dist) + mean(weights.in$dist))
  )
  return(weights.in)
}

# Compute correlation coefficients and p-values
compute.cor <- function(x, y) {
  cor.estimate <- cor.test(x, y, method = "pearson", exact = FALSE)
  return(list(estimate = cor.estimate$estimate, p = cor.estimate$p.value))
}


######################
# Plotting functions #
######################
# Plot correlation between DEL and DUP model weights on each PC
plot.weights.cor <- function(weights.in, cor.in, do.label.outliers,
                             title.in = "PC Weights in DEL vs. DUP Model") {
  x.axis.title <- "DEL"
  y.axis.title <- "DUP"

  # Prep plot area
  par(mar = c(2.5, 2.5, 1.2, 1))
  ax.lims <- c(
    floor(min(weights.in$weight.DEL, weights.in$weight.DUP) * 10) / 10,
    ceiling(max(weights.in$weight.DEL, weights.in$weight.DUP) * 10) / 10
  )

  # Plot points
  plot(
    weights.in$weight.DEL, weights.in$weight.DUP,
    pch = 19, xaxt = "n", xlab = "", ylab = "", yaxt = "n", las = 2,
    xlim = ax.lims, ylim = ax.lims
  )
  # Label outlier points if optioned
  if (do.label.outliers) {
    outlier.weights <- weights.in[weights.in$outlier, ]
    if (nrow(outlier.weights) > 0) {
      text(
        outlier.weights$weight.DEL, outlier.weights$weight.DUP,
        labels = rownames(outlier.weights),
        pos = ifelse(outlier.weights$weight.DEL > outlier.weights$weight.DUP, 4, 2),
        cex = 0.7, col = "red", offset = 0.3
      )
    }
  }

  # Plot y = x line
  abline(0, 1, col = adjustcolor(offblack, alpha = 0.5), lty = 2)

  # Add X axis
  x.ax.at <- axTicks(1)
  axis(1, at = x.ax.at, tck = -0.0125, col = offblack, labels = NA)
  axis(1, at = x.ax.at, tick = F, line = -0.8, cex.axis = 0.8)
  mtext(1, line = 1.25, text = x.axis.title)
  # Add Y axis
  y.ax.at <- axTicks(2)
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
  make_option(c("--label-outliers"),
    help = "Label outlier points",
    action = "store_true", default = FALSE
  ),
  make_option(c("--title"),
    help = "Custom title [default '%default']",
    type = "character", default = "PC Weights in DEL vs. DUP Model"
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
do.label.outliers <- opts$`label-outliers`
title.in <- opts$title

# Load DEL and DUP model scalar weights on PCs
model.weights <- load.weights(del.tsv, dup.tsv)

# Annotate outliers if any
model.weights <- annot.outliers(model.weights)

# Compute correlation between DEL and DUP model weights on PCs
del.dup.cor <- compute.cor(model.weights$weight.DEL, model.weights$weight.DUP)

# Plot DEL vs. DUP mu quantile distribution
pdf(paste(out.prefix, "DEL.DUP.model_weights.pdf", sep = "."),
  height = 5, width = 5
)
plot.weights.cor(
  model.weights, del.dup.cor, do.label.outliers,
  ifelse(is.null(title.in), "PC Weights in DEL vs. DUP Model", title.in)
)
dev.off()

# Write out PCs sorted by difference in DEL and DUP weight
write.table(
  model.weights[order(model.weights$dist, decreasing = TRUE), ],
  paste(out.prefix, "DEL.DUP.model_weights.tsv", sep = "."),
  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE
)
