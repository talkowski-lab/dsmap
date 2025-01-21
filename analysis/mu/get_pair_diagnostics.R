#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Lily Wang <lily_wang@hms.harvard.edu>

# Collect diagnostics for probabilities of observed CNVs in bin-pairs


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
# Summarize distributions of pairs with CNV probabilities
# for an input BED loaded with dsmapR::load.bins()
summarize.pairs.binary <- function(pairs) {
  contigs <- unique(pairs$coords[, 1])
  binsize <- infer.bin.size(pairs$coords)
  sizes <- pairs$coords[, 3] - pairs$coords[, 2] - binsize
  has_sv <- pairs$feats[, 1] >= 0.5

  # Compute number of pairs per contig with and without CNVs
  df.by.contig <- cbind(contigs, as.data.frame(
    do.call("rbind", lapply(contigs, function(contig) {
      c(
        length(which(pairs$coords[, 1] == contig & !(has_sv))),
        length(which(pairs$coords[, 1] == contig & has_sv))
      )
    }))
  ))
  colnames(df.by.contig) <- c("contig", "no_sv", "has_sv")

  # Compute number of pairs per bin size with and without CNVs
  size.range <- seq(0, max(sizes), by = binsize)
  df.by.size <- cbind(size.range / 1000, as.data.frame(
    do.call("rbind", lapply(size.range, function(size) {
      c(
        length(which(sizes == size & !(has_sv))),
        length(which(sizes == size & has_sv))
      )
    }))
  ))
  colnames(df.by.size) <- c("pair_distance_kb", "no_sv", "has_sv")

  # Compute number of pairs per size type (same-bin, adjacent-bin, other)
  # with and without CNVs
  adjacencies <- ifelse(sizes == 0, "same", ifelse(sizes == binsize, "adj", "none"))
  adjacency.types <- unique(adjacencies)
  df.by.adjacency <- cbind(adjacency.types, as.data.frame(
    do.call("rbind", lapply(adjacency.types, function(adjacency) {
      c(
        length(which(adjacencies == adjacency & !(has_sv))),
        length(which(adjacencies == adjacency & has_sv))
      )
    }))
  ))
  colnames(df.by.adjacency) <- c("adjacency", "no_sv", "has_sv")
  df.by.adjacency$adjacency <- factor(
    df.by.adjacency$adjacency,
    levels = c("same", "adj", "none")
  )

  # Compute total number of pairs with and without CNVs
  df.overall <- data.frame(no_sv = sum(df.by.contig$no_sv), has_sv = sum(df.by.contig$has_sv))

  # Compute proportions of pairs with and without CNVs
  dfs <- lapply(
    list(
      "contig" = df.by.contig,
      "size" = df.by.size,
      "adjacency" = df.by.adjacency,
      "overall" = df.overall
    ), function(df) {
      df$prop_has_sv <- df$has_sv / (df$has_sv + df$no_sv)
      return(df)
    }
  )

  return(dfs)
}

# Summarize distributions of pairs with CNV integer counts
# for an input BED loaded with dsmapR::load.bins()
summarize.pairs.integer <- function(pairs) {
  contigs <- unique(pairs$coords[, 1])
  binsize <- infer.bin.size(pairs$coords)
  sizes <- pairs$coords[, 3] - pairs$coords[, 2] - binsize
  n.svs <- pairs$feats[, 1]

  # Compute number of pairs per contig with each number of CNVs
  # Do not fill in zeros
  df.by.contig <- as.data.frame(do.call("rbind", lapply(contigs, function(contig) {
    cbind(contig, as.data.frame(table(n.svs[pairs$coords[, 1] == contig])))
  })))
  colnames(df.by.contig) <- c("contig", "n_svs", "n_pairs")

  # Compute number of pairs by bin size with each number of CNVs
  size.range <- seq(0, max(sizes), by = binsize)
  df.by.size <- as.data.frame(do.call("rbind", lapply(size.range, function(size) {
    cbind(size / 1000, as.data.frame(table(n.svs[sizes == size])))
  })))
  colnames(df.by.size) <- c("pair_distance_kb", "n_svs", "n_pairs")

  # Compute number of pairs per size type (same-bin, adjacent-bin, other)
  # with each number of CNVs
  adjacencies <- ifelse(sizes == 0, "same", ifelse(sizes == binsize, "adj", "none"))
  adjacency.types <- unique(adjacencies)
  df.by.adjacency <- as.data.frame(do.call(
    "rbind",
    lapply(adjacency.types, function(adjacency) {
      cbind(adjacency, as.data.frame(table(n.svs[adjacencies == adjacency])))
    })
  ))
  colnames(df.by.adjacency) <- c("adjacency", "n_svs", "n_pairs")
  df.by.adjacency$adjacency <- factor(
    df.by.adjacency$adjacency,
    levels = c("same", "adj", "none")
  )

  df.overall <- as.data.frame(table(n.svs))
  colnames(df.overall) <- c("n_svs", "n_pairs")
  df.overall$prop_pairs <- df.overall$n_pairs / sum(df.overall$n_pairs)

  # Compute proportions of pairs with each number of CNVs
  dfs <- c(lapply(list(
    "contig" = df.by.contig,
    "size" = df.by.size,
    "adjacency" = df.by.adjacency
  ), function(df) {
    df$prop_pairs <- do.call("c", by(df, df[, 1], function(x) {
      x$n_pairs / sum(x$n_pairs)
    }))
    return(df)
  }), list("overall" = df.overall))

  return(dfs)
}


######################
# Plotting functions #
######################
# Barplots of bin-pair positive vs. negative counts or positive proportion
# Optionally colored by CNV type
plot.counts.binary <- function(df, prop = FALSE, title.in = NA, x.axis.title = NA, cnv = NA,
                               label.all.x.ticks = FALSE, x.label.cex = 1, x.label.las = 1) {
  # Set plotting values
  all.x.labels <- df[, 1]
  n.bars <- nrow(df)
  legend.labs <- c()
  if (prop) {
    if (cnv %in% c("DEL", "DUP", "CNV")) {
      bar.colors <- c(get(paste(cnv, "colors", sep = "."))$main)
    } else {
      bar.colors <- c(browns$main)
    }
    ylims <- c(0, max(df$prop_has_sv))
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
    parmar = c(2.75, 3.2, 1.2, 0.3)
  )

  # Add bars
  if (prop) {
    rect(
      xleft = (1:n.bars) - 1, xright = 1:n.bars, ybottom = 0,
      ytop = df$prop_has_sv, border = "white", col = bar.colors[1]
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
  mtext(1, line = 1.8, text = x.axis.title)

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
  if (prop) {
    y.text <- paste("Prop. pairs with", cnv)
  } else {
    y.text <- "Pairs"
    if (!is.null(units)) {
      y.text <- paste(y.text, " (", units, ")", sep = "")
    }
  }
  mtext(2, line = 2.3, text = y.text)

  # Add title
  mtext(3, font = 2, text = title.in, xpd = T)

  # Add legend
  if (!prop) {
    legend("topright",
      legend = legend.labs, fill = bar.colors,
      cex = 0.85, border = offblack, bg = offwhite, xpd = T
    )
  }
}

# Barplots of integer CNV counts by bin-pair grouping
# Assumes that groupings are defined as factor levels of first column
# and correspond to group.labels
# Optionally colored by CNV type
plot.counts.integer <- function(df, title.in = NA, group.title = NA, group.labels = NA,
                                x.axis.title = NA, cnv = NA, label.all.x.ticks = FALSE,
                                x.label.cex = 1, x.label.las = 1) {
  # Set plotting values
  groups <- levels(df[, 1])
  if (cnv %in% c("DEL", "DUP", "CNV")) {
    bar.colors <- c(get(paste(cnv, "colors", sep = "."))$main)
  } else {
    bar.colors <- c(browns$main)
  }
  all.x.labels <- sort(unique(df$n_svs))
  x.vals <- length(unique(df$n_svs))

  # Prep plotting area
  par(mfrow = c(length(groups), 1), mar = c(2.5, 3.3, 2.6, 1.5), bty = "n")

  # For each grouping, plot number of pairs with each number of CNVs
  # with proportion on y-value, number above bar
  # NOTE: X-labels are potentially non-consecutive integers
  for (i in seq_along(groups)) {
    # Subset to data in this grouping
    group <- groups[i]
    group.df <- df[df[, 1] == group, ]
    group.x.labels <- group.df$n_svs
    group.x.vals <- which(all.x.labels %in% group.x.labels) - 1

    # Create bar plot
    plot(NA,
      xlim = c(0, max(x.vals)), ylim = c(-0.05, 1), type = "n", xaxs = "i",
      xlab = "", xaxt = "n", yaxs = "i", ylab = "", yaxt = "n"
    )
    rect(
      xleft = group.x.vals, xright = group.x.vals + 1, ybottom = 0,
      ytop = group.df$prop_pairs, border = "white", col = bar.colors[1]
    )
    text(group.x.vals + 0.5, group.df$prop_pairs + 0.05, labels = group.df$n_pairs, xpd = T)

    # Add X axis
    x.ticks <- 0:length(all.x.labels)
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
    if (i == length(groups)) {
      mtext(1, line = 1.5, text = x.axis.title)
    }

    # Add Y axis
    y.ax.at <- axTicks(2)
    axis(2, at = c(-10e10, 10e10), col = offblack, tck = 0)
    axis(2, at = y.ax.at, tck = -0.025, col = offblack, labels = NA)
    axis(2, at = y.ax.at, tick = F, line = -0.65, labels = y.ax.at, las = 2)
    y.text <- "Prop. pairs"
    mtext(2, line = 2, text = y.text)

    # Add title
    if (i == 1) {
      mtext(3, font = 2, line = 1.3, text = title.in, xpd = T)
    }
    mtext(3, font = 2, text = paste(group.title, ": ", group.labels[i], sep = ""), xpd = T, cex = 0.85)
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
  ),
  make_option(c("--integer"),
    help = "Count data are integers, rather than probabilities or binary 0/1.",
    action = "store_true", default = FALSE
  )
)

# Get command-line arguments & options
arg_list <- c("pairs.bed", "out.prefix")
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
pairs.in <- args$args[1]
out.prefix <- args$args[2]
cnv <- opts$cnv
counts.are.integers <- opts$integer

# Load pairs
pairs <- load.bins(pairs.in)

if (!counts.are.integers) {
  # Summarize counts
  dat <- summarize.pairs.binary(pairs)

  # Plot bin-pair counts and proportions with CNVs per contig
  for (count_type in c("counts", "props")) {
    pdf(
      paste(
        out.prefix, paste(count_type, "per_contig", sep = "_"),
        "pdf",
        sep = "."
      ),
      height = 2.5, width = 4.25
    )
    plot.counts.binary(
      dat[["contig"]],
      prop = (count_type == "props"),
      title.in = "Bin-pairs",
      x.axis.title = "Chromosome", cnv = cnv,
      label.all.x.ticks = T, x.label.cex = 0.85, x.label.las = 2
    )
    dev.off()
  }

  # Plot bin-pair counts and proportions with CNVs vs. pair distance
  for (count_type in c("counts", "props")) {
    pdf(
      paste(
        out.prefix, paste(count_type, "vs_distance", sep = "_"),
        "pdf",
        sep = "."
      ),
      height = 2.5, width = 4.25
    )
    plot.counts.binary(
      dat[["size"]],
      prop = (count_type == "props"),
      title.in = "Bin-pairs",
      x.axis.title = "Pair distance (kb)", cnv = cnv
    )
    dev.off()
  }

  # Plot bin-pair counts and proportions with CNVs vs. bin adjacency in pair
  for (count_type in c("counts", "props")) {
    pdf(
      paste(
        out.prefix, paste(count_type, "by_adjacency", sep = "_"),
        "pdf",
        sep = "."
      ),
      height = 2.5, width = 2
    )
    # Adjust x-labels
    adj.dat <- dat[["adjacency"]]
    levels(adj.dat$adjacency) <- c("Same", "Adj.", "None")
    plot.counts.binary(
      adj.dat,
      prop = (count_type == "props"),
      title.in = "Bin-pairs",
      x.axis.title = "Pair adjacency", cnv = cnv,
      label.all.x.ticks = T, x.label.cex = 0.85
    )
    dev.off()
  }
} else {
  # Summarize counts
  dat <- summarize.pairs.integer(pairs)

  # Plot counts by size type
  pdf(
    paste(
      out.prefix, "counts_by_adjacency",
      "pdf",
      sep = "."
    ),
    height = 5, width = 7
  )
  plot.counts.integer(
    dat[["adjacency"]],
    title.in = paste(cnv, "counts by bin-pair adjacency"),
    group.title = "Adjacency",
    group.labels = c("Same", "Adjacent", "None"),
    x.axis.title = paste(cnv, "count in bin-pair"), cnv = cnv,
    label.all.x.ticks = T, x.label.cex = 0.85
  )
  dev.off()

  # TODO: Plot counts by pair repeat content
  # NOTE: No plot for bin-pair counts with CNVs per chromosome -
  # check text file to find these counts
}

group.filenames <- list(
  "overall" = "overall_counts.tsv",
  "contig" = "counts_per_contig.tsv",
  "size" = "counts_vs_distance.tsv",
  "adjacency" = "counts_by_adjacency.tsv"
)

# Write counts by each grouping to output TSVs
for (group in names(dat)) {
  write.table(dat[[group]], paste(out.prefix, group.filenames[[group]], sep = "."),
    row.names = F, col.names = T, sep = "\t", quote = F
  )
}
