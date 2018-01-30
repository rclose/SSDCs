#function to sever links between countries or regions in MSTs
severRegionalLinks <- function(mst.tree, dat, from, to) {

	countries <- as.data.frame(unique(dat[, c("collection_no","country")])); rownames(countries) <- countries$collection_no; countries$collection_no <- NULL
	severed_tree <- mst.tree[!((countries[mst.tree$from,] %in% from & countries[mst.tree$to,] %in% to) | (countries[mst.tree$from,] %in% to & countries[mst.tree$to,] %in% from)), ]
	sub_trees <- findAllSubtrees(mst.tree = severed_tree)
	names(sub_trees) <- paste(round(sapply(sub_trees, function(x) sum(x$distance))), "km", sep = "")
	return(sub_trees)

}
