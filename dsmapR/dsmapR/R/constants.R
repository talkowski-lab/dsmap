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
#' @param susbet Name of constant subset to return. See `Details` for options.
#' @details Recognized values for `subset` include:
#' * `colors` : all color palettes used
#' * `cnv_colors` : color palettes for CNVs (subset of `colors`)
#' * `other_colors` : all other color constants not included in CNV palettes (subset of `colors`)
#' @examples
#' # Load list of color palettes
#' dsmap.colors <- get_constants("colors");
#' @return Returns a vector of constants corresponding to the value of `subset`
#' @export
get_constants <- function(subset){
  # Define colors
  cnv.colors <- list("DEL" = list("dark2" = "#4F1C14",
                                  "dark1" = "#9F2B1C",
                                  "main" = "#D43925",
                                  "light1" = "#DD6151",
                                  "light2" = "#E5887C"),
                     "DUP" = list("dark2" = "#003F6A",
                                  "dark1" = "#1A5985",
                                  "main" = "#2376B2",
                                  "light1" = "#4F91C1",
                                  "light2" = "#7BADD1"),
                     "CNV" = list("dark2" = "#3F2759",
                                  "dark1" = "#5F3B85",
                                  "main" = "#7E4EB2",
                                  "light1" = "#9B6081",
                                  "light2" = "#B488A1"))
  other.colors <- list("green" = list("dark2" = "#0C4923",
                                      "dark1" = "#116D34",
                                      "main" = "#179145",
                                      "light1" = "#45A76A",
                                      "light2" = "#74BD8F"),
                       "orange" = list("dark2" = "#7E3809",
                                       "dark1" = "#BD530E",
                                       "main" = "#FC6F12",
                                       "light1" = "#FD8C41",
                                       "light2" = "#FDA971"),
                       "brown" = list("dark2" = "#6A4838",
                                      "dark1" = "#9F6C54",
                                      "main" = "#D49070",
                                      "light1" = "#DDA68D",
                                      "light2" = "#E5BCA9"),
                       "offblack" = "#3D2921",
                       "offwhite" = "#F9F0EB")

  # Return values based on input
  if(subset == "colors"){
    c(cnv.colors, other.colors)
  }else if(subset == "cnv_colors"){
    cnv.colors
  }else if(subset == "other_colors"){
    other.colors
  }
}

