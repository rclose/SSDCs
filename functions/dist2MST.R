#' Convert Matrix of Geographic Great Circle Distances to Minimum Spanning Tree
#'
#' Function that Takes a Pairwise Matrix of Great-Circle Distances and Returns a Minimum Spanning Tree.
#' @param gcd.matrix A geographic great-circle distance matrix containing pairwise distances between points.
#' @keywords mst gcd gcd.matrix
#' @export
#' @examples
#' binaryMST2edges(mst.tree = mst.tree, gcd.matrix = gcd.matrix)

dist2MST <- function(gcd.matrix) {

	if (length(gcd.matrix) > 1) { #@TODO better conditional handling required

		# cat(paste("Calculating MST for ", nrow(gcd.matrix), " collections... "), sep = "")

		gcd.matrix <- as.dist(t(gcd.matrix))
		mst.bin <- ape::mst(gcd.matrix)
		mst.bin[upper.tri(mst.bin, diag = TRUE)] <- 0 #zap the upper triangle since it duplicates the lower
		mst.tree <- binaryMST2edges(mst.bin, gcd.matrix)
		mst.tree <- mst.tree[order(mst.tree$distance), ]

		# cat("Done.\n")

	} else {mst.tree <- NA}

	return(mst.tree)

}


