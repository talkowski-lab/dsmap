#!/usr/bin/env R

#######################
#    DSMap Project    #
#######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# DSMap plotting basics


#' Prepare plot area
#'
#' Prepare a standardized & formatted plot area
#' @param xlims Range of values for X axis
#' @param ylims Range of values for Y axis
#' @param parmar Margin values passed to par()
#' @examples
#' prep_plot_area(xlims=c(0, 5), ylims=(-10, 10), parmar=rep(3, 4));
#' @export
prep_plot_area <- function(xlims, ylims, parmar){
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims, type="n",
       xaxs="i", xlab="", xaxt="n",
       yaxs="i", ylab="", yaxt="n")
}
