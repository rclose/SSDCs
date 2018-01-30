#' Split Minimum Spanning Tree At Link
#'
#' Function that takes a minimum-spanning tree and a pair of collection numbers representing a link and splits the object into two, returning both subtrees.
#' @param mst.tree MST object
#' @param connection Vector containing two collection names
#' @keywords key 
#' @export
#' @examples
#' example()

splitMST <- function(mst.tree, connection = NULL, remove.longest.connection = FALSE) {

	if (length(mst.tree) > 1) {

		mst.tree1 <- mst.tree

		if (remove.longest.connection == TRUE) {

			longest.row <- which.max(mst.tree$distance)
			connection <- mst.tree[longest.row, 2:3]

		}

		split.mst <- mst.tree[!(mst.tree[ , 3] %in% connection & mst.tree[ , 2] %in% connection), ]

		mst.subtree.1 <- seedMST(mst.tree = split.mst, collection = connection[1])
		mst.subtree.2 <- seedMST(mst.tree = split.mst, collection = connection[2])

	}

	out <- list(mst.subtree.1, mst.subtree.2); names(out) <- c(paste(connection[1], "subtree", sep = "_"), paste(connection[2], "subtree", sep = "_"))

	return(out)

}


