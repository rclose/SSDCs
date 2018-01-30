#' Find All Nodes Within Mst Connected To Starting Collection
#'
#' Function that takes a minimum spanning tree and a collection name and returns the subtree of connected nodes
#' @param mst.tree MST object
#' @param collection Name of collection from which to start growing MST
#' @keywords key 
#' @export
#' @examples
#' example()

seedMST <- function(mst.tree, collection) {

	if (length(mst.tree) > 1) {

		sub.mst <- matrix(nrow = 0, ncol = 3); colnames(sub.mst) <- colnames(mst.tree)

		difference <- 1

		while (difference > 0) {

			total <- length(collection)
			sub.mst <- mst.tree[unique(c(which(mst.tree[ ,2] %in% collection), which(mst.tree[ ,3] %in% collection))), ]
			collection <- unique(c(collection, sub.mst[ ,2], sub.mst[ ,3]))
			difference <- length(collection) - total

		}

		return(sub.mst)

	} else {print("MST object is probably blank")}


}

