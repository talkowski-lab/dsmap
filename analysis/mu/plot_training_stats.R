#!/usr/bin/env Rscript

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot training/testing performance for an athena mutation rate model


#########
# Setup #
#########
# Load necessary libraries
require(dsmapR, quietly=TRUE)
require(beeswarm, quietly=TRUE)
require(optparse, quietly=TRUE)

# Set global options and constants
options(stringsAsFactors=FALSE, scipen=1000)
dsmapR::load.constants(c("colors", "names"))


##################
# Data functions #
##################
# Load a table of training statistics from athena mu-train
load.stats <- function(stats.in){
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")
  colnames(stats)[1] <- gsub("^X\\.", "", colnames(stats)[1])
  cols.to.keep <- c("accuracy", "sensitivity", "specificity", "ppv", "npv", "f1")
  list("train.cv" = stats[stats$stage == "train" & stats$test_chrom != "None", cols.to.keep],
       "train.all" = stats[stats$stage == "train" & stats$test_chrom == "None", cols.to.keep],
       "test.cv" = stats[stats$stage == "test", cols.to.keep])
}


######################
# Plotting functions #
######################
# Plot a beeswarm of a single performance metric
swarm.single <- function(vals, x.at, color=offblack, pch=1, cex=1, width=0.8,
                         mean.bar=TRUE, orient.mean="right", mean.color=offblack){
  # Set parameters
  if(orient.mean == "right"){
    mean.pos <- 4
    mean.adj <- -1
  }else if(orient.mean == "left"){
    mean.pos <- 2
    mean.adj <- 1
  }else{
    stop(paste("swarm.single() does not recognize '", orient.mean,
               "' as a valid value for orient.mean", sep=""))
  }

  # Add mean bar, if optioned
  if(mean.bar){
    half.bar.width <- width/3
    vals.mean <- mean(vals, na.rm=T)
    rect(xleft=x.at-half.bar.width, xright=x.at+half.bar.width,
         ybottom=0, ytop=vals.mean, col=mean.color, xpd=T, border=NA)
    # segments(x0=x.at-half.bar.width, x1=x.at+half.bar.width,
    #          y0=vals.mean, y1=vals.mean,
    #          lwd=3, lend="round", col=mean.color, xpd=T)
  }

  # Swarm of values
  par(xpd=TRUE)
  beeswarm(vals, at=x.at, add=T, xlab="", ylab="", xaxt="n", yaxt="n",
           pch=pch, corral="wrap", corralWidth=width, col=color, cex=cex)
  par(xpd=FALSE)

  # Print mean
  if(mean.bar){
    text(x=x.at+(mean.adj*width), y=vals.mean+0.03,
         pos=mean.pos, cex=0.85, xpd=T,
         labels=format(round(vals.mean, 2), nsmall=2))
  }
}

# Wrapper function to plot all training statistics
plot.stats <- function(stats, cnv=NULL){
  # Set parameters
  n.metrics <- ncol(stats$train.cv)
  if(is.null(cnv)){
    colorset <- browns
  }else{
    if(cnv == "DEL"){
      colorset <- DEL.colors
    }else{
      colorset <- DUP.colors
    }
  }

  # Prep plot area
  prep.plot.area(xlims=c(0, n.metrics+2), ylims=c(0, 1), parmar=c(1.25, 1.5, 2.5, 0.1))

  # Add axes & title
  abline(v=0:n.metrics, col=offwhite)
  axis(1, at=c(-10e10, n.metrics), tck=0, col=offwhite, labels=NA)
  axis(3, at=c(-10e10, n.metrics), tck=0, col=offwhite, labels=NA)
  axis(2, at=seq(0, 1, 0.2), col=offblack, labels=NA, tck=-0.025)
  axis(2, at=seq(0,1, 0.2), tick=F, las=2, cex.axis=0.85, line=-0.65)
  axis(3, at=n.metrics/2, line=0.5, tick=F, font=2, labels="Model Performance", xpd=T)

  # Add metrics
  sapply(1:n.metrics, function(i){
    # Add group X-axis
    axis(1, at=i-c(0.85, 0.15), tck=0, labels=NA, col=offblack)
    axis(1, at=i-0.5, tick=F, line=-0.65,
         labels=performance.metric.abbrevs[colnames(stats$train.cv)[i]])

    # Get values
    train.vals <- stats$train.cv[, i]
    train.all.val <- stats$train.all[1, i]
    test.vals <- stats$test.cv[, i]

    # Plot values
    swarm.single(train.vals, x.at=i-(2/3), width=(1/4), cex=0.75,
                 color=colorset$light1, pch=1, mean.color=colorset$dark1,
                 orient.mean="left")
    points(x=i-(2/3), y=train.all.val, pch=18, col=offblack)
    swarm.single(test.vals, x.at=i-(1/3), width=(1/4), cex=0.75,
                 color=colorset$light2, pch=19, mean.color=colorset$dark1,
                 orient.mean="right")
  })

  # Add legend
  legend("right", legend=c("Train (CV)", "Train (All)", "Test (CV)", "Mean (CV)"),
         col=c(colorset$light1, offblack, colorset$light2, colorset$dark1),
         pch=c(1, 18, 19, 15), pt.cex=1.3, border=offblack, bg=offwhite,
         cex=0.85, xpd=T)
}


###########
# RScript #
###########
# List of command-line options
option_list <- list(
  make_option(c("--cnv"), help="Specify CNV type. Used for plotting colors only.",
              type="character", default=NA)
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog stats.tsv out.pdf",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: stats.tsv and output.pdf\n", sep=" "))
}

# Writes args & opts to vars
stats.in <- args$args[1]
out.pdf <- args$args[2]
cnv <- opts$cnv

# Load stats
stats <- load.stats(stats.in)

# Plot stats
pdf(out.pdf, height=3, width=5)
plot.stats(stats, cnv)
dev.off()
