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
#' prep.plot.area(xlims=c(0, 5), ylims=(-10, 10), parmar=rep(3, 4));
#' @export prep.plot.area
#' @export
prep.plot.area <- function(xlims, ylims, parmar){
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims, type="n",
       xaxs="i", xlab="", xaxt="n",
       yaxs="i", ylab="", yaxt="n")
}


#' Rotate points
#'
#' Rotate one or more points in 2D Cartesian space
#' @param coords Two-column matrix or data.frame of original x & y coordinates
#' @param angle Angle to rotate, in degrees
#' @param x.origin X coordinate to rotate about (default: 0)
#' @param y.origin Y coordinate to rotate about (default: 0)
#' @examples
#' rotate.points(coords=data.frame("x"=c(2, 5), "y"=c(0, 3)), angle=90)
#' @return Returns a two-column data.frame of rotated x & y coordinates
#' @export
rotate.points <- function(coords, angle, x.origin=0, y.origin=0){
  # Process each point
  res <- as.data.frame(t(apply(coords, 1, function(xy){
    # Assign coordinates to variables
    x <- as.numeric(xy[1]); y <- as.numeric(xy[2])

    # Convert new angle degrees to radians
    new.angle <- angle * (pi / 180)

    # Infer current angle from (x - x.o, y - y.o)
    old.angle <- atan2(x - x.origin, y - y.origin)

    # Theta = old.angle + new.angle
    theta <- old.angle + new.angle

    # Compute distance from point to origin
    r <- sqrt(((x - x.origin) ^ 2) + ((y - y.origin) ^ 2))

    # New x = (r * cos(theta)) + x.origin
    # New y = (r * sin(theta)) + y.origin
    c(r * cos(theta) + x.origin,
      r * sin(theta) + y.origin)
  })))

  colnames(res) <- c("x", "y")

  return(res)
}

#' Plot diagonal heatmap
#'
#' Plot the upper-right diagonal of a numeric square matrix as a rotated heatmap
#' @param mat Numeric square matrix (or data.frame) of values to be plotted
#' @param orient Specify orientation of plot as "up" diagonal pointing up or
#' "down" for diagonal pointing down (default: "up")
#' @param gridlines Boolean indicator whether to add gridlines (default: true)
#' @param gridline.lwd Width of gridlines (default: 0.1)
#' @param parmar Margin values passed to par()
#' @seealso [dsmapR::prep.plot.area()], [dsmapR::rotate.points()]
#' @export plot.diag.heat
#' @export
plot.diag.heat <- function(mat, orient="up", gridlines=TRUE, gridline.lwd=0.25,
                           parmar=c(0.5, 0.5, 0.5, 0.5)){
  # Get plot parameters
  n.rows <- nrow(mat)
  x.max <- rotate.points(matrix(c(n.rows, n.rows), nrow=1), angle=-45)$x
  y.max <- rotate.points(matrix(c(0, n.rows), nrow=1), angle=45)$y
  if(orient == "up"){
    ylims <- c(0, y.max)
  }else{
    ylims <- c(y.max, 0)
  }
  heat.pal <- viridis(101)

  # Rescale values to match palette
  val.range <- range(mat[which(!is.na(mat) & !is.infinite(mat))])
  val.scaled <- ceiling((mat - val.range[1]) * (100 / diff(val.range))) + 1

  # Prep plot area
  prep.plot.area(xlims=c(0, x.max), ylims=ylims, parmar=parmar)

  # Add rectangles for each value
  sapply(1:n.rows, function(x){
    sapply(1:n.rows, function(y){
      rect.coords <- data.frame("x"=c(x-1, x, x, x-1),
                                "y"=c(y-1, y-1, y, y))
      rotated.coords <- rotate.points(rect.coords, angle=-45)
      polygon(x=rotated.coords$x, y=rotated.coords$y,
              border=heat.pal[val.scaled[y, x]],
              col=heat.pal[val.scaled[y, x]],
              lwd=gridline.lwd)
    })
  })

  # Add thin gridlines, if optioned
  if(gridlines){
    sapply(0:n.rows, function(i){
      gridline.coords <- rotate.points(data.frame("x"=c(i, i, 0, n.rows),
                                                  "y"=c(0, n.rows, i, i)),
                                       angle=-45)
      segments(x0=gridline.coords$x[c(1, 3)],
               x1=gridline.coords$x[c(2, 4)],
               y0=gridline.coords$y[c(1, 3)],
               y1=gridline.coords$y[c(2, 4)],
               col=colors$offwhite, lwd=gridline.lwd)
    })
  }
}

#' Add scale bar to diagonal heatmap
#'
#' Add a 1D scale bar to a diagonal heatmap produced by [plot.diag.heat()]
#' @param bin.size Width of each bin
#' @param units Suffix for scale bar (default: "kb")
#' @param x.at Relative position on x-axis for center of scale bar (default: 0.8)
#' @param y.at Relative position on y-axis for center of scale bar (default: 0.8)
#' @param width Relative width of scale bar (default: 0.2)
#' @param lwd Line width for scale bar (default: 3)
#' @param label.vadj Relative vertical adjustment factor for label spacing (default: 0.025)
#' @details `width` will be rounded to the nearest whole value of bin.size
#' @examples
#' plot.diag.heat(matrix(rnorm(100), nrow=10))
#' add.scale.bar(bin.size=1000, units="kb", label.vadj=-0.01)
#' @seealso [dsmapR::plot.diag.heat()]
#' @export add.scale.bar
#' @export
add.scale.bar <- function(bin.size, units, x.at=0.8, y.at=0.8, width=0.2, lwd=3,
                          label.vadj=0.025){
  rounded.width <- round(width * diff(par("usr")[1:2]))
  x.mid <- par("usr")[1] + (x.at * diff(par("usr")[1:2]))
  x.lr <- x.mid + (c(-0.5, 0.5) * rounded.width)
  y <- min(par("usr")[3:4]) + (y.at * abs(diff(par("usr")[3:4])))
  segments(x0=x.lr[1], x1=x.lr[2], y0=y, y1=y, lwd=lwd, lend="butt")
  label <- paste(prettyNum(rounded.width, big.mark=","), units)
  text(x=x.mid, y=y+(label.vadj * abs(diff(par("usr")[3:4]))), pos=3, labels=label)
}

