#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot mutation rate model calibration from a .tsv of predicted rates and observed SVs


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly=TRUE)
require(optparse, quietly=TRUE)

# Set global options and constants
options(stringsAsFactors=FALSE, scipen=1000)
dsmapR::load.constants(c("colors", "scales"))


##################
# Data functions #
##################
# Load data and compute calibration
load.calibration <- function(calibration.tsv, min.points.per.bin=100, max.bins=100){
  # Read data & log-transform
  dat <- read.table(calibration.tsv, sep="\t", comment.char="", header=T)
  colnames(dat) <- c("predicted", "actual")
  dat$predicted <- log10(dat$predicted + (1/nrow(dat)))

  # Determine number of bins and bin breakpoints
  n.points <- nrow(dat)
  points.per.bin <- max(c(n.points / max.bins, min.points.per.bin))
  pct.step <- points.per.bin / n.points
  n.breaks <- (1 / pct.step) + 1
  breaks <- quantile(dat$predicted, probs=seq(0, 1, length.out=n.breaks))

  # Format calibration dataframe
  calibration.df <- as.data.frame(t(sapply(1:(n.breaks-1), function(i){
    keep.idx <- which(dat$predicted >= breaks[i] & dat$predicted < breaks[i+1])
    c(mean(dat$predicted[keep.idx], na.rm=T),
      log10(mean(c(dat$actual[keep.idx], 1), na.rm=T)))
  })))
  colnames(calibration.df) <- c("predicted", "actual")
  return(calibration.df)
}


######################
# Plotting functions #
######################
# Plot model calibration
plot.calibration <- function(dat, cnv=NULL){
  # Get parameters
  min.val <- floor(min(dat, na.rm=T))
  if(!is.null(cnv)){
    colors <- get(paste(cnv, "colors", sep="."))
  }else{
    colors <- brown
  }

  # Set up plot area
  prep.plot.area(xlims=c(min.val, 0), ylims=c(min.val, 0), parmar=c(3, 3, 1, 1),
                 xaxs="r", yaxs="r")
  abline(0, 1, col=colors$dark2, lty=5, lwd=2)

  # Add axes
  sapply(1:2, function(side){
    axis(side, at=log10(logscale.major), labels=NA, col=offblack, tck=-0.02)
    sapply(logscale.major, function(k){
      label <- if(nchar(k) > 4){bquote(10^.(log10(k)))}else{k}
      axis(side, at=log10(k), tick=F, line=-0.65, labels=label,
           cex.axis=5.5/6, las=if(side==2){2}else{1})
    })
    mtext(side, line=2,
          text=paste(if(side==1){"Predicted"}else{"Observed"}, "Probability"))
  })
  mtext(3, text="Model Calibration", font=2)

  # Add points
  polygon(x=c(dat$predicted, rev(dat$predicted)),
          y=c(dat$actual, rev(dat$predicted)),
          col=adjustcolor(colors$light2, alpha=0.5), border=NA)
  points(dat, col=colors$main, type="l")
  points(dat, pch=19, cex=0.7, col=colors$main)

  # Add legend
  legend("bottomright", legend=c("Actual", "Ideal"),
         col=c(colors$main, colors$dark2), lty=c(NA, 5), lwd=c(0, 2), pch=c(19, NA),
         pt.cex=1.3, border=offblack, bg=offwhite, cex=0.85, xpd=T)
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
  make_option(c("--cnv"), help="Specify CNV type. Used for plotting colors only.",
              type="character", default=NA),
  make_option(c("--min-points-per-bin"), default=100, type="numeric",
              help="Minimum number of points per calibration bin."),
  make_option(c("--max-bins"), default=100, type="numeric",
              help="Maximum number of calibration bins.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog calibration.tsv out.pdf",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: calibration.tsv and output.pdf\n", sep=" "))
}

# Writes args & opts to vars
calibration.tsv <- args$args[1]
out.pdf <- args$args[2]
cnv <- opts$cnv
min.points.per.bin <- opts$`min-points-per-bin`
max.bins <- opts$`max-bins`

# Load calibration data
dat <- load.calibration(calibration.tsv, min.points.per.bin, max.bins)

# Plot calibration
pdf(out.pdf, height=4, width=4)
plot.calibration(dat, cnv)
dev.off()
