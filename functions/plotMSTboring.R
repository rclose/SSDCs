#' Plot Boring Blank MST Visualisation
#'
#' Plot MST points and connections on blank background
#' @param palaeocoords Palaeocoordinates from PBDB occurrence file
#' @param mst.tree MST object
#' @param dots Colour of nodes in MST; defaults to red
#' @param lines Colour of lines in MST; defaults to grey
#' @keywords key 
#' @export
#' @examples
#' plotMSTboring(palaeocoords, mst.tree, dots = "red", lines = "grey")

plotMSTboring <- function(palaeocoords, mst.tree, dots = "red", lines = "grey") {

	palaeocoords <- palaeocoords[, c(2,1)] #this just ensures that lat and long are plotted on the right axes
	pairwise.palaeocoords <- mst.tree[,2:3]
	pairwise.palaeocoords[mst.tree$from, ]

	plot(palaeocoords, pch = 19, cex = 0.5, col = dots)
	text(palaeocoords, label = rownames(palaeocoords), cex = 0.2, pos = 2)

	for (i in 1:nrow(mst.tree)) {

		x0 <- palaeocoords[pairwise.palaeocoords[i,1], 1]
		y0 <- palaeocoords[pairwise.palaeocoords[i,1], 2]

		x1 <- palaeocoords[pairwise.palaeocoords[i,2], 1]
		y1 <- palaeocoords[pairwise.palaeocoords[i,2], 2]

		segments(x0, y0, x1, y1, col = lines)

	}

}


