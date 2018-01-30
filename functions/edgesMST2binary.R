#' Function To Convert Distance/From/To MST Format To Binary Matrix Format Used By APE
#'
#' Function that Takes a geoMST Format MST and Returns a Binary MST
#' @param mst.tree A minimum-spanning tree obect in edge-matrix format
#' @keywords mst gcd gcd.matrix
#' @export
#' @examples
#' edgesMST2binary(mst.tree)

edgesMST2binary <- function(mst.tree) {

	localities <- unique(c(mst.tree[ ,2], mst.tree[ ,3]))

	mst.binary <- matrix(0, nrow = length(localities), ncol = length(localities))

	colnames(mst.binary) <- rownames(mst.binary) <- localities

	for (i in 1:nrow(mst.tree)) {

		mst.binary[mst.tree[i,2], mst.tree[i,3]] <- 1

	}

	class(mst.binary) <- "mst"

	return(mst.binary)

}


