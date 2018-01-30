PGS <- function(pcoords, round.mst.pcoords = 1, grid.cell.size = 2) {

	#DEPENDENCIES
	#This function relies on a number of other packages and scripts.
	#These include the alphahull, fossil and ape packages, Roger Close's geoMST functions and the max_great_circle.R function from Matt Clapham.
	#The subfunctions in the latter were shifted outside his function to make them available for use here.

	#remove rows that don't contain palaeocoordinates
	pcoords <- unique(as.data.frame(na.omit(pcoords)))

	if (nrow(pcoords) <= 1) {
		PGSout <- cbind(
			# NA,
			NA,
			NA,
			NA,
			NA,
			NA,
			NA)
		colnames(PGSout) <- c(
			# "convex_hull",
			"mst_dist",
			"grid_cell_occupancy",
			"max_gcd",
			"median_pairwise_GCD",
			"mean_pairwise_GCD",
			"STD")
		return(PGSout)
	}

	# #calculate minimum convex polygon using the fossil package
	# if (nrow(pcoords) >= 3) {
	# 	earth_poly_out <- try(fossil::earth.poly(lats = cbind(pcoords$paleolngdec,pcoords$paleolatdec)))
	# 	if (class(earth_poly_out) == "try-error") {
	# 		earth_poly_out <- list(area = NA, vertices = NA)
	# 	}
	# } else if (nrow(pcoords == 2)) {earth_poly_out <- list(area = NA, vertices = NA)}

	#calculate MST
	#round the palaeocoordinates to bin by lat/lng cell, if required

	if ("paleolat" %in% colnames(pcoords) & !("paleolatdec" %in% colnames(pcoords))) {
	pcoords$paleolngdec <- pcoords$paleolng
	pcoords$paleolatdec <- pcoords$paleolat
	}

	mst.pcoords <- pcoords
	if (!is.null(round.mst.pcoords)) {
		mst.pcoords$paleolngdec <- plyr::round_any(pcoords$paleolngdec, round.mst.pcoords)
		mst.pcoords$paleolatdec <- plyr::round_any(pcoords$paleolatdec, round.mst.pcoords)
		mst.pcoords <- unique(mst.pcoords)
	}

	if (nrow(mst.pcoords) >= 2) {
		mst_out <- MST(mst.pcoords)
		mst_dist <- sum(mst_out$distance)
	} else {
		mst_dist <- NA
	}

	#count occupied grid-cells (number of unique rounded palaeocoordinates)
	grid.cell.pcoords <- pcoords
	grid.cell.pcoords$paleolngdec <- plyr::round_any(pcoords$paleolngdec, grid.cell.size)
	grid.cell.pcoords$paleolatdec <- plyr::round_any(pcoords$paleolatdec, grid.cell.size)
	grid_occupancy_out <- nrow(unique(mst.pcoords))

	#calculate median pairwise great-circle distances
	pairwise_GCDs <- pcoord2gcd(grid.cell.pcoords)
	median_pairwise_GCD <- median(pairwise_GCDs, na.rm = T)
	mean_pairwise_GCD <- mean(pairwise_GCDs, na.rm = T)

	# #calculate maximum great-circle distance using Matt Clapham's functions
	# max_gcd_out <- max.dist(latdata = pcoords$paleolat, longdata = pcoords$paleolng)

	#calculate maximum great-circle distance using my functions
	max_gcd_out <- max(pairwise_GCDs, na.rm = T)

	#calculate standard distance deviation (using formula from Burt et al. 2009, 'Elementary Statistics for Geographers')
	centroid <- apply(X = grid.cell.pcoords, MARGIN = 2, FUN = mean)
	STD <- sqrt((sum(apply(X = grid.cell.pcoords, MARGIN = 1, FUN = function(x) (x[1] - centroid[1])^2)) / nrow(grid.cell.pcoords)) + (sum(apply(X = grid.cell.pcoords, MARGIN = 1, FUN = function(x) (x[2] - centroid[2])^2)) / nrow(grid.cell.pcoords)))

	#return results
	# return(list(earth_poly = earth_poly_out$area, mst = mst_dist$mst.distance, grid_occupancy = grid_occupancy_out$occupied.grids, max_gcd = max_gcd_out, alpha_hull_area = alpha_hull_area))
	# return(list(convex_hull = earth_poly_out$area, mst = mst_dist$mst.distance, no_occupied_grid_cells = grid_occupancy_out$occupied.grids, max_gcd = max_gcd_out))
	PGSout <- cbind(
		# earth_poly_out$area,
		mst_dist,
		grid_occupancy_out,
		max_gcd_out,
		median_pairwise_GCD,
		mean_pairwise_GCD,
		STD)
	colnames(PGSout) <- c(
		# "convex_hull",
		"mst_dist",
		"grid_cell_occupancy",
		"max_gcd",
		"median_pairwise_GCD",
		"mean_pairwise_GCD",
		"STD")
	return(PGSout)

}
