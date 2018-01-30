#' Convert Palaeocoordinates to Pairwise Matrix of Great Circle Distances
#'
#' Function to convert Nx2 matrix of palaeocoordinates to NxN pairwise matrix of Great Circle Distances
#' @param palaeocoords Palaeocoordinates in Nx2 matrix
#' @keywords key
#' @export
#' @examples
#' example()

pcoord2gcd <- function(palaeocoords) {

	if (nrow(palaeocoords) > 1) {

		cat(paste("Calculating GCD matrix for ", nrow(palaeocoords), " palaeocoordinates... ", sep = ""))

		gcd.matrix <- matrix(nrow = nrow(palaeocoords), ncol = nrow(palaeocoords))

		rownames(gcd.matrix) <- colnames(gcd.matrix) <- rownames(palaeocoords)

		for (j in 1:(nrow(palaeocoords) - 1)) {

			for (k in (j + 1):nrow(palaeocoords)) {

				gcd.matrix[j,k] <- gcd(
					x = palaeocoords[[c(1,j)]],
					y = palaeocoords[[c(1,k)]],
					z = abs(palaeocoords[[c(2,j)]] - palaeocoords[[c(2,k)]])
				)

			}
		}
		cat("Done.\n")
	} else {gcd.matrix <- NA}

	return(gcd.matrix)

}

