#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2024-Present Lily Wang and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>

# Collect diagnostics for probabilities of observed CNVs in genome bins

# TODO: Handle quantitative counts that are not binary/probabilities


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly = TRUE)
require(optparse, quietly = TRUE)
require(stringr, quietly = TRUE)

# Set global options and constants
options(stringsAsFactors = FALSE, scipen = 1000)
dsmapR::load.constants("colors")


##################
# Data functions #
##################
# Summarize distributions of bins with CNV probabilities
# for an input BED loaded with dsmapR::load.bins()
summarize.bins <- function(bins) {
  contigs <- unique(bins$coords[, 1])
  has_sv <- bins$feats[, 1] >= 0.5

  # Compute dataframe of counts per contig based on SV overlap
  df <- cbind(contigs, as.data.frame(
    do.call("rbind", lapply(contigs, function(contig) {
      c(
        length(which(bins$coords[, 1] == contig & !(has_sv))),
        length(which(bins$coords[, 1] == contig & has_sv))
      )
    }))
  ))
  colnames(df) <- c("contig", "no_sv", "has_sv")
  df$pct_has_sv <- df$has_sv / (df$has_sv + df$no_sv)
}


######################
# Plotting functions #
######################
# Barplots of bin positive vs. negative counts or positive percentage
# Optionally colored by CNV type
plot.counts <- function(df, pct = FALSE, title = NA, x.axis.title = NA, cnv = NA,
                        label.all.x.ticks = FALSE, x.label.cex = 1, x.label.las = 1) {
  # Set plotting values
  all.x.labels <- df[, 1]
  n.bars <- nrow(df)
  legend.labs <- c()
  if (pct) {
    if (cnv %in% c("DEL", "DUP", "CNV")) {
      bar.colors <- c(get(paste(cnv, "colors", sep = "."))$main)
    } else {
      bar.colors <- c(browns$main)
    }
    ylims <- c(0, max(df$pct_has_sv))
  } else {
    if (cnv %in% c("DEL", "DUP", "CNV")) {
      bar.colors <- c(
        get(paste(cnv, "colors", sep = "."))$light1,
        get(paste(cnv, "colors", sep = "."))$dark1
      )
      legend.labs <- paste(c("No", "Has"), cnv)
    } else {
      bar.colors <- c(browns$light1, browns$dark1)
      legend.labs <- c("No SV", "Has SV")
    }
    ylims <- c(0, max(df$has_sv + df$no_sv))
  }

  # Prep plotting area
  prep.plot.area(
    xlims = c(0, n.bars), ylims = ylims,
    parmar = c(2.1, 3.2, 1.2, 0.3)
  )

  # Add bars
  if (pct) {
    rect(
      xleft = (1:n.bars) - 1, xright = 1:n.bars, ybottom = 0,
      ytop = df$pct_has_sv, border = "white", col = bar.colors[1]
    )
  } else {
    rect(
      xleft = (1:n.bars) - 1, xright = 1:n.bars, ybottom = 0, ytop = df$has_sv,
      border = NA, col = bar.colors[2]
    )
    rect(
      xleft = (1:n.bars) - 1, xright = 1:n.bars, ybottom = df$has_sv,
      ytop = df$has_sv + df$no_sv, border = NA, col = bar.colors[1]
    )
    rect(
      xleft = (1:n.bars) - 1, xright = 1:n.bars,
      ybottom = 0, ytop = df$has_sv + df$no_sv, border = "white", col = NA
    )
  }

  # Add X axis
  if (label.all.x.ticks) {
    x.ticks <- 0:length(all.x.labels)
  } else {
    x.ticks <- axTicks(1)
  }
  x.ax.at <- x.ticks[-length(x.ticks)] + 0.5
  x.ax.labels <- all.x.labels[x.ax.at + 0.5]
  axis(1, at = c(-10e10, 10e10), col = offblack, tck = 0)
  axis(1, at = x.ax.at, tck = -0.025, col = offblack, labels = NA)
  sapply(seq_along(x.ax.at), function(x) {
    axis(1,
      at = x.ax.at[x], tick = F, line = -0.65, labels = x.ax.labels[x],
      cex.axis = x.label.cex, las = x.label.las
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
  if (pct) {
    y.text <- paste("Proportion bins with", cnv)
  } else {
    y.text <- "Bins"
    if (!is.null(units)) {
      y.text <- paste(y.text, " (", units, ")", sep = "")
    }
  }
  mtext(2, line = 2.3, text = y.text)

  # Add title
  mtext(3, font = 2, text = title, xpd = T)

  # Add legend
  if (!pct) {
    legend("topright",
      legend = legend.labs, fill = bar.colors,
      cex = 0.85, border = offblack, bg = offwhite, xpd = T
    )
  }
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
arg_list <- c("bins.bed", "out.prefix")
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
bins.in <- args$args[1]
out.prefix <- args$args[2]
cnv <- opts$cnv

# Load bins
bins <- load.bins(bins.in)

# Summarize counts
dat <- summarize.bins(bins)

# Merge counts per contig & write to output file
write.table(dat[["contig"]], paste(out.prefix, "counts_per_contig.tsv", sep = "."),
  row.names = F, col.names = T, sep = "\t", quote = F
)

# Plot bin counts and percentages with CNVs per contig
for (count_type in c("counts", "pcts")) {
  pdf(
    paste(
      out.prefix, paste(count_type, "per_contig", sep = "_"),
      "pdf",
      sep = "."
    ),
    height = 2.5, width = 4.25
  )
  plot.counts(
    dat,
    pct = (count_type == "pcts"),
    title = "Bins",
    x.axis.title = "Chromosome", cnv = cnv,
    label.all.x.ticks = T, x.label.cex = 0.85, x.label.las = 2
  )
  dev.off()
}
