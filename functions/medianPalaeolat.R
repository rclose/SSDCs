#' Median Palaeolatitudes
#'
#' Get median palaeolatitudes
#' @param palaeocoords Matrix of palaeocoordinates from PBDB occurrence file
#' @keywords palaeocoordiantes median palaeolatitude 
#' @export
#' @examples
#' medianPalaeolat(palaeocoords)

medianPalaeolat <- function(palaeocoords) {

	if (nrow(palaeocoords) > 1) {

		centroid <- median(palaeocoords[ ,"paleolatdec"])

	} else {centroid <- NA}

	return(centroid)

}

