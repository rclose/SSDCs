#' Find Subtrees in Fragmented MST Object
#'
#' Takes an MST object in geoMST format and returns a list object containing all disconnected subtrees.
#' @param mst.tree MST object in geoMST format
#' @keywords MST subtrees 
#' @export
#' @examples
#' findAllSubtrees(mst.tree)

#### findAllSubtrees: FUNCTION THAT FINDS ALL DISCONNECTED SUBTREES WITHIN AN MST OBJECT ####
findAllSubtrees <- function(mst.tree) {
	if (ncol(mst.tree) == 3 && nrow(mst.tree) > 1) {
		i <- 0
		mst.tree.bak <- mst.tree
		mst.subtrees <- list()
		while (nrow(mst.tree) > 1) {
			i <- i + 1
			mst.subtrees[[i]] <- seedMST(mst.tree = mst.tree, collection = mst.tree[1,2])
			remaining.nodes <- unique(c(mst.tree[,2], mst.tree[,3]))
			subtree.nodes <- unique(c(mst.subtrees[[i]][,2], mst.subtrees[[i]][,3]))
			keep <- setdiff(remaining.nodes, subtree.nodes)
			mst.tree <- mst.tree[mst.tree[ , 3] %in% keep & mst.tree[ , 2] %in% keep, ]
		}
		mst.subtrees <- mst.subtrees[sapply(mst.subtrees, nrow) > 1]
		return(mst.subtrees)
	} else {cat("MST malformed")}
}


