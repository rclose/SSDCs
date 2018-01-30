#' Calculate Minimum Spanning Tree from Palaeocoordinates
#'
#' Function to run all steps and return a Minimum Spanning Tree and summed distances
#' @param palaeocoordinates Matrix of Palaeocoordinates from PBDB occurrence file
#' @param gcd.matrix Pairwise matrix of Great Circle Distances; calculated if not supplied
#' @keywords MST palaeocoordinates
#' @export
#' @examples
#' example()

MST <- function(palaeocoords, gcd.matrix = NULL) {

	if (c("paleolat","paleolng") %in% colnames(palaeocoords) && !(c("paleolatdec","paleolngdec") %in% colnames(palaeocoords))) {
		palaeocoords$paleolatdec <- palaeocoords$paleolat
		palaeocoords$paleolngdec <- palaeocoords$paleolng
		palaeocoords <- palaeocoords[,c("paleolatdec","paleolngdec")]
	} else {
		palaeocoords <- palaeocoords[,c("paleolatdec","paleolngdec")]
	}

	# # Get median palaeolatitudes
	# centroid <- medianPalaeolat(palaeocoords)

	# Convert Nx2 matrix of palaeocoordinates to NxN pairwise matrix of Great Circle Distances
	if (is.null(gcd.matrix)) {
		gcd.matrix <- pcoord2gcd(palaeocoords)
	}

	# Convert GCD matrix into Nx3 ordered interlocality distance matrix and then into to Mx3 minimum spanning tree
	mst.tree <- dist2MST(gcd.matrix)

	return(mst.tree)

}
