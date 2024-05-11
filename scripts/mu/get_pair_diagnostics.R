#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Collect diagnostics for bin-pairs prior to training mutation rate model


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly = TRUE)
require(optparse, quietly = TRUE)

# Set global options and constants
options(stringsAsFactors = FALSE, scipen = 1000)
dsmapR::load.constants("colors")


##################
# Data functions #
##################
# Summarize distributions of pairs for an input BED loaded with dsmapR::load.pairs()
summarize.pairs <- function(pairs) {
  contigs <- unique(pairs$coords[, 1])
  binsize <- infer.bin.size(pairs$coords)
  sizes <- pairs$coords[, 3] - pairs$coords[, 2] - binsize
  has_sv <- pairs$feats[, 1] >= 0.5

  # Compute dataframe of counts per contig based on SV overlap
  df.by.contig <- do.call("rbind", lapply(contigs, function(contig) {
    c(
      length(which(pairs$coords[, 1] == contig & !(has_sv))),
      length(which(pairs$coords[, 1] == contig & has_sv))
    )
  }))
  rownames(df.by.contig) <- contigs
  colnames(df.by.contig) <- c("no_sv", "has_sv")

  # Compute dataframe of counts by bin size
  size.range <- seq(0, max(sizes), by = binsize)
  df.by.size <- do.call("rbind", lapply(size.range, function(size) {
    c(
      length(which(sizes == size & !(has_sv))),
      length(which(sizes == size & has_sv))
    )
  }))
  rownames(df.by.size) <- size.range / 1000
  colnames(df.by.size) <- c("no_sv", "has_sv")

  return(list(
    "df.by.contig" = df.by.contig,
    "df.by.size" = df.by.size
  ))
}


######################
# Plotting functions #
######################
# Barplots of bin-pair positive vs. negative counts or positive percentage
# Optionally colored by CNV type
plot.counts <- function(df, pct = FALSE, title = NA, x.axis.title = NA, cnv = NA,
                        label.all.x.ticks = FALSE, x.label.cex = 1, x.label.las = 1) {
  # Set plotting values
  all.x.labels <- rownames(df)
  n.bars <- nrow(df)
  legend.labs <- c()
  if (pct) {
    df$pct_has_sv <- df$has_sv / (df$has_sv + df$no_sv)
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
    parmar = c(2.1, 2.8, 1.2, 0.3)
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
  sapply(1:length(x.ax.at), function(x) {
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
  y.text <- "Bin-pairs"
  if (!is.null(units)) {
    y.text <- paste(y.text, " (", units, ")", sep = "")
  }
  mtext(2, line = 1.9, text = ytext)

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
args <- parse_args(
  OptionParser(
    usage = "%prog pairs.bed training.bed out_prefix",
    option_list = option_list
  ),
  positional_arguments = TRUE
)
opts <- args$options

# Checks for appropriate positional arguments
if (length(args$args) != 3) {
  stop(paste("Three positional arguments required: pairs.bed, training.bed, and output_prefix\n", sep = " "))
}

# Writes args & opts to vars
pairs.in <- args$args[1]
training.in <- args$args[2]
out.prefix <- args$args[3]
cnv <- opts$cnv

# Load pairs
pairs <- load.bins(pairs.in)
training <- load.bins(training.in)

# Summarize counts
pairs.dat <- summarize.pairs(pairs)
training.dat <- summarize.pairs(training)

# Merge counts per contig & write to output file
df.by.contig <- merge(pairs.dat$df.by.contig, training.dat$df.by.contig,
  by = 0, suffixes = c("_all", "_training")
)
colnames(df.by.contig)[1] <- "contig"
write.table(df.by.contig, paste(out.prefix, "pair_counts_per_contig.tsv", sep = "."),
  row.names = F, col.names = T, sep = "\t", quote = F
)

# Merge counts by size & write to output file
df.by.size <- merge(pairs.dat$df.by.size, training.dat$df.by.size,
  by = 0, suffixes = c("_all", "_training"), sort = F
)
colnames(df.by.size)[1] <- "pair_distance_kb"
df.by.size <- df.by.size[order(as.numeric(df.by.size$pair_distance_kb)), ]
write.table(df.by.size, paste(out.prefix, "pair_counts_vs_distance.tsv", sep = "."),
  row.names = F, col.names = T, sep = "\t", quote = F
)

# Plot counts per contig
pdf(paste(out.prefix, "pair_counts_per_contig.all.pdf", sep = "."),
  height = 2.5, width = 4.25
)
plot.counts(pairs.dat$df.by.contig,
  pct = FALSE, title = "All bin-pairs",
  x.axis.title = "Chromosome", cnv = cnv,
  label.all.x.ticks = T, x.label.cex = 0.85, x.label.las = 2
)
dev.off()
pdf(paste(out.prefix, "pair_pcts_per_contig.all.pdf", sep = "."),
  height = 2.5, width = 4.25
)
plot.counts(pairs.dat$df.by.contig,
  pct = TRUE, title = "All bin-pairs",
  x.axis.title = "Chromosome", cnv = cnv,
  label.all.x.ticks = T, x.label.cex = 0.85, x.label.las = 2
)
dev.off()
pdf(paste(out.prefix, "pair_counts_per_contig.training.pdf", sep = "."),
  height = 2.5, width = 4.25
)
plot.counts(training.dat$df.by.contig,
  pct = FALSE, title = "Training bin-pairs",
  x.axis.title = "Chromosome", cnv = cnv,
  label.all.x.ticks = T, x.label.cex = 0.85, x.label.las = 2
)
dev.off()
pdf(paste(out.prefix, "pair_pcts_per_contig.training.pdf", sep = "."),
  height = 2.5, width = 4.25
)
plot.counts(training.dat$df.by.contig,
  pct = TRUE, title = "Training bin-pairs",
  x.axis.title = "Chromosome", cnv = cnv,
  label.all.x.ticks = T, x.label.cex = 0.85, x.label.las = 2
)
dev.off()

# Plot counts vs pair distance contig
pdf(paste(out.prefix, "pair_counts_vs.distance.all.pdf", sep = "."),
  height = 2.5, width = 4.25
)
plot.counts(pairs.dat$df.by.size,
  title = "All bin-pairs",
  x.axis.title = "Pair distance (kb)", cnv = cnv
)
dev.off()
pdf(paste(out.prefix, "pair_counts_vs.distance.training.pdf", sep = "."),
  height = 2.5, width = 4.25
)
plot.counts(training.dat$df.by.size,
  title = "Training bin-pairs",
  x.axis.title = "Pair distance (kb)", cnv = cnv
)
dev.off()
