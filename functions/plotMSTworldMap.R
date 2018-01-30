#' Plot MST on Palaeomap
#'
#' Plot MST And Connections With Map (not working at all)
#' @param arg Description
#' @keywords key 
#' @export
#' @examples
#' example()

plotMSTworldMap <- function(palaeocoords, mst.tree) {

	library(maps)

	pairwise.palaeocoords <- mst.tree[,2:3]

	map('world', col = "grey90", fill = TRUE)
	points(palaeocoords, pch = 19, cex = 0.5, col = "red")

	for (i in 1:nrow(mst.tree)) {

		x0 <- palaeocoords[pairwise.palaeocoords[i,1], 1]
		y0 <- palaeocoords[pairwise.palaeocoords[i,1], 2]

		x1 <- palaeocoords[pairwise.palaeocoords[i,2], 1]
		y1 <- palaeocoords[pairwise.palaeocoords[i,2], 2]

		segments(x0, y0, x1, y1)

	}

}

