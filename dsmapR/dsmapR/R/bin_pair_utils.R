#!/usr/bin/env R

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Functions for handling I/O and simple manipulations of athena bins and bin-pairs


#' Load athena bins or bin-pairs
#'
#' Load a BED file of athena 1D bins or 2D bin-pairs
#' @param bed_in Path to input BED file
#' @param n_keep_features Number of feature columns to keep. Defaults to keeping
#' all feature columns present in `bed_in`.
#' @param feats_are_numeric Boolean indicating if all features should be coerced
#' to numeric type. `TRUE` by default.
#' @return Returns a list of two dataframes:
#' 1. `$coords`: a three-column dataframe of bin coordinates
#' 2. `$feats`: a dataframe of bin features
#' @export
load.bins <- function(bed_in, n_keep_features=Inf, feats_are_numeric=TRUE){
  # Read full bed_in
  x <- read.table(bed_in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(x)[1] <- gsub("#", "", colnames(x)[1], fixed=T)

  # Split out coordinates
  coords <- x[, 1:3]
  coords[, 2:3] <- apply(coords[, 2:3], 2, as.numeric)

  # Split out features (if any)
  n_feat_cols <- ncol(x) - 3
  if(n_feat_cols == 0){
    feats <- NULL
  }else{
    feat_cidxs_to_keep <- seq(4, 3 + min(c(n_feat_cols, n_keep_features)))
    if(feats_are_numeric == TRUE){
      feats <- as.data.frame(apply(data.frame(x[, feat_cidxs_to_keep]), 2, as.numeric))
    }else{
      feats <- data.frame(x[, feat_cidxs_to_keep])
    }
    colnames(feats) <- colnames(x)[feat_cidxs_to_keep]
  }

  return(list("coords" = coords, "feats" = feats))
}


#' Infer bin size
#'
#' Infer bin size from BED file loaded as dataframe
#' @param coords Dataframe styled as BED (first three columns are coordinates)
#' @param n_samp number of rows to sample from dataframe. Default = 1,000.
#' @examples
#' # Load coordinates from a BED file & infer bin size
#' bins <- dsmapR::load.bins("path/to/my/bins.bed")
#' binsize <- infer.bin.size(bins$coords)
#' @return Integer estimate of inferred bin size
#' @export
infer.bin.size <- function(coords, n_samp=1000){
  contig <- coords[1, 1]
  starts <- head(sort(unique(coords[which(coords[, 1] == contig), 2])), n_samp)
  min(starts[2:length(starts)] - starts[(2:length(starts))-1])
}

