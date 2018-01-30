#' Convert GCD Matrix to Ordered Interlocality Distance Table
#'
#' Function to convert NxN pairwise Great Circle Distance matrix to Nx3 ordered interlocality distance table (matrix)
#' @param arg Description
#' @keywords key 
#' @export
#' @examples
#' example()

gcd2intDist <- function(gcd.matrix) {

	if (nrow(gcd.matrix) > 2) {

		dist.table <- as.data.frame(as.table(gcd.matrix, stringsAsFactors = FALSE))
		dist.table <- setNames(dist.table[complete.cases(dist.table) & dist.table[ ,"Freq"] != 0, c(3,1:2)], c("distance","from","to"))
		dist.table <- dist.table[order(dist.table$distance), ]

	} else {dist.table <- NA}

	return(dist.table)

}

