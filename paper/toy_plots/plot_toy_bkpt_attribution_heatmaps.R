#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Make toy mutation rate heatmaps for breakpoint precision & attribution schematics


# Load libraries
require(dsmapR, quietly=TRUE)
require(optparse, quietly=TRUE)
require(viridis, quietly=TRUE)

# Set global options and constants
options(stringsAsFactors=FALSE, scipen=1000)
colors <- dsmapR::get.constants("colors")


# Distribute breakpoint probability across 1D bins
attribute.bkpt <- function(bin.starts, bin.width, bkpt.pos, bkpt.stdev){
  sapply(bin.starts, function(x){
    lower <- pnorm(x, mean=bkpt.pos, sd=bkpt.stdev)
    upper <- pnorm(x + bin.width, mean=bkpt.pos, sd=bkpt.stdev)
    upper - lower
  })
}

# Simulate a matrix of breakpoint probabilities
sim.bkpt.probs <- function(matrix.start, matrix.end, bin.width, sv.df){
  # sv.df must be three-column data.frame with start, end, and bkpt stdev for each SV

  # Build bins
  bin.starts <- seq(matrix.start, matrix.end - bin.width, by=bin.width)
  n.bins <- length(bin.starts)

  # Make single matrix with P(no SV) for each cell
  p.mat <- 1 - Reduce("*", lapply(1:nrow(sv.df), function(i){
    p.start <- attribute.bkpt(bin.starts, bin.width,
                              bkpt.pos=sv.df[i, 1], bkpt.stdev=sv.df[i, 3])
    p.end <- attribute.bkpt(bin.starts, bin.width,
                            bkpt.pos=sv.df[i, 2], bkpt.stdev=sv.df[i, 3])
    1 - crossprod(matrix(p.start, nrow=1), t(matrix(p.end, ncol=1)))
  }))
  rownames(p.mat) <- colnames(p.mat) <- bin.starts

  return(p.mat)
}

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to vars
out.prefix <- args$args[1]

# Heatmap for single SV example
single.bkpt.prob <- sim.bkpt.probs(matrix.start=0, matrix.end=10000, bin.width=1000,
                                   sv.df=data.frame("start"=3500, "end"=6500, "stdev"=1000))
pdf(paste(out.prefix, "single_bkpt_example", "diag_heat", "pdf", sep="."),
    width=2.8, height=1.4)
plot.diag.heat(single.bkpt.prob, orient="down", parmar=rep(0, 4))
add.scale.bar(bin.size=1000, units="kb")
dev.off()

# Heatmap for multi-SV example:
set.seed(2021)
multi.sv.df <- data.frame("start"=c(1300, 2700, 3900, 6200, 12250, 13700, 15300),
                          "end"=c(8750, 10800, 7250, 17500, 14600, 20700, 22800),
                          "stdev"=rnorm(mean=1000, sd=90, n=7))
multi.bkpt.prob <- sim.bkpt.probs(matrix.start=0, matrix.end=25000, bin.width=1000,
                                  sv.df=multi.sv.df)
pdf(paste(out.prefix, "multi_bkpt_example", "diag_heat", "pdf", sep="."),
    width=4, height=2)
plot.diag.heat(multi.bkpt.prob, orient="down", parmar=rep(0, 4))
add.scale.bar(bin.size=1000, units="kb")
dev.off()

# Guide for multi-SV example:
pdf(paste(out.prefix, "multi_bkpt_example", "guide", "pdf", sep="."),
    width=4, height=1)
prep.plot.area(xlims=c(0, 25000), ylims=c(0, nrow(multi.sv.df)), parmar=rep(0, 4))
segments(x0=multi.sv.df$start, x1=multi.sv.df$end,
         y0=(1:nrow(multi.sv.df))-0.5, y1=(1:nrow(multi.sv.df))-0.5)
rect(xleft=multi.sv.df$start-(2*multi.sv.df$stdev),
     xright=multi.sv.df$start+(2*multi.sv.df$stdev),
     ybottom=(1:nrow(multi.sv.df))-0.65,
     ytop=(1:nrow(multi.sv.df))-0.35,
     col=adjustcolor("black", alpha=0.5))
rect(xleft=multi.sv.df$end-(2*multi.sv.df$stdev),
     xright=multi.sv.df$end+(2*multi.sv.df$stdev),
     ybottom=(1:nrow(multi.sv.df))-0.65,
     ytop=(1:nrow(multi.sv.df))-0.35,
     col=adjustcolor("black", alpha=0.5))
dev.off()
