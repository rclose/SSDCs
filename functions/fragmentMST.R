fragmentMST <- function(mst.tree, max.subtree.size = 12000, min.subtree.size = 0, method = c("iteratively.remove.longest.branch","split.iteratively","remove.outlier.branches"), outlier.sd = 2) {

	method <- match.arg(method, c("iteratively.remove.longest.branch","split.iteratively","remove.outlier.branches"))

	if (sum(mst.tree$distance) < max.subtree.size) {

		cat("MST already below max.subtree.size\n")

		mst.tree <- setNames(list(mst.tree), paste(round(sum(mst.tree$distance)), "km", sep = ""))

		return(mst.tree)

	} else {

		subtrees <- list()

		cat("Splitting MST of size ", sum(mst.tree$distance), " to obtain a maximum subtree size of ", max.subtree.size, " km and minumum of ", min.subtree.size, " km\n", sep = "")

		# OPTION 1: break up global MST by figuring out how many of the longest branches to remove to reduce subtrees below max.subtree.size
		if (method == "iteratively.remove.longest.branch") {

			prune <- 0
			subtree.lengths <- sum(mst.tree$distance)

			while (max(subtree.lengths) > max.subtree.size) {

				prune <- prune + 1
				subtrees <- findAllSubtrees(mst.tree[1:(nrow(mst.tree) - prune),])
				subtree.lengths <- sapply(subtrees, function(x) sum(x$distance))
				cat(paste(length(subtrees), "subtrees, largest is", round(max(subtree.lengths)), "km\n", sep = " "))

			}

		}


		#OPTION 2: break up global MST by iteratively splitting into subtrees until all are under max.subtree.size
		if (method == "split.iteratively") {

			subtrees <- splitMST(mst.tree, remove.longest.connection = TRUE) #split the global MST into two subtrees
			subtree.lengths <- sapply(subtrees, function(x) sum(x$distance))
			cat(paste(length(subtrees), "subtrees, largest is", round(max(subtree.lengths)), "km\n", sep = " "))

			while (max(subtree.lengths) > max.subtree.size ) {

				subtrees <- unlist(lapply(subtrees, function(x) splitMST2(x, split.if.above = max.subtree.size)), recursive = FALSE) #flatten nested list to make each element an MST object
				keep <- as.vector(sapply(subtrees, nrow)) #get number of segments in subtrees
				subtrees <- subtrees[keep > 2] #discard any subtrees smaller than 3 segments
				subtree.lengths <- sapply(subtrees, function(x) sum(x$distance)) #get lengths of subtrees
				names(subtrees) <- NULL #wipe the names because otherwise they grow to silly lengths
				cat(paste(length(subtrees), "subtrees, largest is", round(max(subtree.lengths)), "km\n\n", sep = " "))

			}

		}


		# #OPTION 3: break up global MST by removing branches greater than two standard deviations from the mean
		if (method == "remove.outlier.branches") {

			x <- mst.tree$distance
			outliers <- x - median(x) > (outlier.sd*sd(x)) #identify large outlier segments via standard deviations from mean
			X <- mst.tree[!outliers, ] #drop large outlier segments
			subtrees <- findAllSubtrees(X)

		}

		subtree.lengths <- sapply(subtrees, function(x) sum(x$distance)) #calculate subtree lengths
		names(subtrees) <- as.character(paste(round(subtree.lengths), "km", sep = "")) #name subtrees by length
		subtrees <- subtrees[subtree.lengths > min.subtree.size] #drop MSTs below the target spread size
		keep <- sapply(subtrees, function(x) nrow(x)) #get numbers of subtree connections
		subtrees <- subtrees[keep > 1] #drop subtrees with fewer than 2 connections

		return(subtrees)

	}

}
