#' Get Palaeocoordinates from PBDB Occurrence Data
#'
#' Extract Nx2 matrix of palaeocoordinates from PBDB data for all intervals
#' @param arg data PBDB occurrence data
#' @param arg intervals Matrix of time bins with LAD and FAD dates
#' @param arg unique.pcoords If TRUE, only return unique palaeocoordinates; otherwise return all collections even if some have identical palaeocoords.
#' @keywords key 
#' @export
#' @examples
#' example()

getPalaeoCoords <- function(data, intervals, unique.pcoords = FALSE, occs.entirely.within.bin = TRUE) {

	palaeocoords <- vector("list", length = nrow(intervals)) #list to receive Nx2 matrices of palaeocoordinates for each interval

	data <- as.data.frame(data)

	for (i in 1:nrow(intervals)) {
		
		if (occs.entirely.within.bin == TRUE) {
		int_data <- data[data$ma_min >= intervals[i,"LAD"] & data$ma_max <= intervals[i,"FAD"], ]
		} else {int_data <- data[data$ma_min < intervals[i,"FAD"] & data$ma_max > intervals[i,"LAD"], ]}
				
		unique_colls <- unique(int_data[c("collection_no","paleolatdec","paleolngdec")])
		rownames(unique_colls) <- unique_colls$collection_no; unique_colls$collection_no <- NULL
		if (unique.pcoords == TRUE) {unique_colls <- unique(unique_colls)}
		palaeocoords[[i]] <- unique_colls

	}

	names(palaeocoords) <- intervals$bin #name each list element with interval abbreviation

	return(palaeocoords)

}

