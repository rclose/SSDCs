#' Convert Binary MST to Edge Matrix Format
#'
#' This function converts MSTs in the binary format as used by package{ape} and package{fossil} to the edge matrix format used by geoMST containing the columns 'distance', 'to' and 'from'.
#' @param mst.tree A minimum spanning tree object in binary format.
#' @param gcd.matrix A geographic great-circle distance matrix containing pairwise distances between points.
#' @keywords mst binary binaryMST2edges
#' @export
#' @examples
#' binaryMST2edges(mst.tree = mst.tree, gcd.matrix = gcd.matrix)

binaryMST2edges <- function(mst.tree, gcd.matrix) {

	mst.mat <- matrix(nrow = 0, ncol = 2); colnames(mst.mat) <- c("from","to")

	for (i in 1:nrow(mst.tree)) { # @TODO: only copy the lower or upper triangle into new object!

		for (j in 1:ncol(mst.tree)) {

			if (mst.tree[i,j] == 1) {

				mst.mat <- rbind(mst.mat, cbind(colnames(mst.tree)[i], rownames(mst.tree)[j]))

			}
		}

	}


	gcd.mat <- as.matrix(as.dist(t(as.table(gcd.matrix))))
	distance <- vector(mode = "numeric")

	for (i in 1:nrow(mst.mat)) {

		distance[i] <- gcd.mat[mst.mat[i,1], mst.mat[i,2]]

	}

	mst.mat <- as.data.frame(cbind(distance, mst.mat), stringsAsFactors = FALSE)
	mst.mat$distance <- as.numeric(mst.mat$distance)
	mst.mat <- mst.mat[order(mst.mat$distance), ]
	rownames(mst.mat) <- NULL

	return(mst.mat)

}


