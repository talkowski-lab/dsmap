#!/usr/bin/env R

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Constants for DSMap project


#' Load DSMap constants
#'
#' Load a subset of constants for DSMap analyses
#'
#' @param subset Vector of constant groups to load See `Details` for options.
#'
#' @details Recognized values for `subset` include:
#' * `colors` : all color palettes used
#' * `cnv_colors` : color palettes for CNVs (subset of `colors`)
#' * `other_colors` : all other color constants not included in CNV palettes (subset of `colors`)
#' * `scales` : all scales and scale labels
#' * `names` : names of various variables
#' * `all` : load all constants
#'
#' @examples
#' # Load list of color palettes
#' get.constants("colors");
#'
#' # Load scales and CNV colors
#' get.constants(c("scales", "cnv_colors"))
#'
#' @export load.constants
#' @export
load.constants <- function(subset) {
  # Define colors
  cnv.colors <- list(
    "DEL.colors" = list(
      "dark2" = "#4F1C14",
      "dark1" = "#9F2B1C",
      "main" = "#D43925",
      "light1" = "#DD6151",
      "light2" = "#E5887C"
    ),
    "DUP.colors" = list(
      "dark2" = "#003F6A",
      "dark1" = "#1A5985",
      "main" = "#2376B2",
      "light1" = "#4F91C1",
      "light2" = "#7BADD1"
    ),
    "CNV.colors" = list(
      "dark2" = "#3F2759",
      "dark1" = "#5F3B85",
      "main" = "#7E4EB2",
      "light1" = "#9B6081",
      "light2" = "#B488A1"
    )
  )
  other.colors <- list(
    "greens" = list(
      "dark2" = "#0C4923",
      "dark1" = "#116D34",
      "main" = "#179145",
      "light1" = "#45A76A",
      "light2" = "#74BD8F"
    ),
    "oranges" = list(
      "dark2" = "#7E3809",
      "dark1" = "#BD530E",
      "main" = "#FC6F12",
      "light1" = "#FD8C41",
      "light2" = "#FDA971"
    ),
    "browns" = list(
      "dark2" = "#6A4838",
      "dark1" = "#9F6C54",
      "main" = "#D49070",
      "light1" = "#DDA68D",
      "light2" = "#E5BCA9"
    ),
    "offblack" = "#3D2921",
    "offwhite" = "#F9F0EB"
  )

  # Define scales
  logscale.major <- 10^(-10:10)
  scales <- list(
    "logscale.major" = logscale.major,
    "logscale.major.bp" = 10^(0:9),
    "logscale.major.bp.labels" = c(
      sapply(
        c("bp", "kb", "Mb"),
        function(suf) {
          paste(c(1, 10, 100), suf, sep = "")
        }
      ),
      "1 Gb"
    ),
    "logscale.demi" = as.numeric(sapply(logscale.major, function(e) {
      c(1, 5) * e
    })),
    "logscale.demi.bp" = as.numeric(sapply(10^(0:9), function(e) {
      c(1, 5) * e
    })),
    "logscale.demi.bp.labels" = c(
      paste(c(1, 5, 10, 50, 100, 500), "bp", sep = ""),
      paste(c(1, 5, 10, 50, 100, 500), "kb", sep = ""),
      paste(c(1, 5, 10, 50, 100, 500), "Mb", sep = ""),
      paste(c(1, 5), "Gb", sep = "")
    ),
    "logscale.minor" = as.numeric(sapply(logscale.major, function(e) {
      (1:9) * e
    }))
  )

  # Define names
  dsmap.names <- list(
    "performance.metric.abbrevs" = c(
      "accuracy" = "Acc.",
      "sensitivity" = "Sens.",
      "specificity" = "Spec.",
      "ppv" = "PPV",
      "npv" = "NPV",
      "f1" = "F1"
    )
  )

  # Assign constants to global environment
  if (length(intersect(subset, c("cnv_colors", "colors", "all"))) > 0) {
    for (variable in names(cnv.colors)) {
      assign(variable, cnv.colors[[variable]], envir = .GlobalEnv)
    }
  }
  if (length(intersect(subset, c("other_colors", "colors", "all"))) > 0) {
    for (variable in names(other.colors)) {
      assign(variable, other.colors[[variable]], envir = .GlobalEnv)
    }
  }
  if (length(intersect(subset, c("scales", "all"))) > 0) {
    for (variable in names(scales)) {
      assign(variable, scales[[variable]], envir = .GlobalEnv)
    }
  }
  if (length(intersect(subset, c("names", "all"))) > 0) {
    for (variable in names(dsmap.names)) {
      assign(variable, dsmap.names[[variable]], envir = .GlobalEnv)
    }
  }
}
