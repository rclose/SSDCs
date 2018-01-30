#### LOAD PACKAGES AND FUNCTIONS ####
library(paleobioDB)
library(abind)
library(fossil)
library(ape)
library(geosphere)
library(vioplot)
library(lattice)
library(RColorBrewer)
library(moments)
library(viridis)
library(stats4)
library(countrycode)
library(geoscale)
library(RCurl)
library(animation)
library(plotrix)
library(plyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(stringr)
library(gridExtra)
library(alphahull)
library(iNEXT)
library(Hmisc)
library(dggridR)
library(png)
library(RCurl)
library(tidyverse)
library(ggrepel)
select <- dplyr::select
slice <- dplyr::slice

file.sources <- list.files("./functions", pattern = "*.R", full.names = TRUE)
sapply(file.sources, source)

#### function for finding ylim when missing data is present ####
ylim.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm = T), 1)

#### function for dropping outdated occurrence identifications for each time slice ####
dropObsoleteOccs <- function(X) {
	getObsoleteOccs <- function(X, reID) {
		Y <- X[X$occurrence_no == reID, ]
		if (length(unique(Y$ref_pubyr)) > 1) {
			to.drop <- Y[-which.max(Y$ref_pubyr), ]$rowID
		} else {
			to.drop <- Y[-which.max(Y$reid_no), ]$rowID
		}
		return(to.drop)
	}

	reids <- unique(X[duplicated(X$occurrence_no), ]$occurrence_no)
	if (length(reids) > 0) {
		to.drop <- mclapply(reids, function(x) {getObsoleteOccs(X, x)})
		to.drop <- unique(as.vector(unlist(to.drop)))
		X <- X[(!X$rowID %in% to.drop), ]
	}
	return(X)
}

#### functions to lighten and darken colours ####
darken <- function(color, factor = 1.4){
	col <- col2rgb(color)
	col <- col/factor
	col <- rgb(t(col), maxColorValue = 255)
	col
}

lighten <- function(color, factor = 1.4){
	col <- col2rgb(color)
	col <- col*factor
	col <- rgb(t(col), maxColorValue = 255)
	col
}


#### function to make stacked plot ####
#plot.stacked makes a stacked plot where each y series is plotted on top
#of the each other using filled polygons
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first")
#  "as.is" - plot in order of y column
#  "max" - plot in order of when each y series reaches maximum value
#  "first" - plot in order of when each y series first value > 0
#'col' - fill colors for polygons corresponding to y columns (will recycle)
#'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#'lwd' - border line width for polygons corresponding to y columns (will recycle)
#'...' - other plot arguments
plot.stacked <- function(
	x, y,
	order.method = "as.is",
	ylab="", xlab="",
	border = NULL, lwd=1,
	col=rainbow(length(y[1,])),
	ylim=NULL,
	...
){
	if(sum(y < 0) > 0) error("y cannot contain negative numbers")
	if(is.null(border)) border <- par("fg")
	border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
	col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
	lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))
	if(order.method == "max") {
		ord <- order(apply(y, 2, which.max))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}
	if(order.method == "first") {
		ord <- order(apply(y, 2, function(x) min(which(x>0))))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}
	top.old <- x*0
	polys <- vector(mode="list", ncol(y))
	for(i in seq(polys)){
		top.new <- top.old + y[,i]
		polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
		top.old <- top.new
	}
	if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
	plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
	for(i in seq(polys)){
		polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
	}
}

pause = function()
{
	if (interactive()) {
		invisible(readline(prompt = "Press <Enter> to continue..."))
	}
	else {
		cat("Press <Enter> to continue...")
		invisible(readLines(file("stdin"), 1))
	}
}

myLetters <- function(length.out) { # make long sequence of letters for subpanel labels
	a <- rep(letters, length.out = length.out)
	grp <- cumsum(a == "a")
	vapply(seq_along(a),
	       function(x) paste(rep(a[x], grp[x]), collapse = ""),
	       character(1L))
}

