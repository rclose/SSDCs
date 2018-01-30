#' Convert Interlocality Distance Table to Minimum Spanning Tree
#'
#' Function to convert an Nx3 ordered interlocality distance table (matrix) into an Mx3 minimum spanning tree
#' @param dist.table Nx3 ordered interlocality distance table
#' @keywords interlocality distance table MST 
#' @export
#' @examples
#' intDist2mst(dist.table)

intDist2mst <- function(dist.table) {

	mst <- NA
	mst.tree <- matrix(nrow = 0, ncol = 3) #To contain a summary of what the network actually looks like (distances, froms, tos)

	if (!is.null(nrow(dist.table))) {

		mst <- 0 #To contain the running total
		i <- 1 #counter

		while (!is.na(i)) {

			mst <- mst + dist.table[i,1] #add the current smallest pairwise distance to the running total

			mst.tree <- rbind(mst.tree, dist.table[i, ]) #rbind the current row of dist.table to the bottom of mst.tree

			dist.table[dist.table == dist.table[i,3]] <- dist.table[i,2] #get the 'from' locality for the current pairwise distance and use it to overwrite any of the 'to' localities in dist.table, as that connection has been made

			i <- which(dist.table[ ,2] != dist.table[ ,3])[1] #get row number of the next shortest connection that hasn't yet been made

		}

	}

	mst.out <- list(mst, mst.tree); names(mst.out) <- c("mst.dist","mst.tree")

	return(mst.out)

}


