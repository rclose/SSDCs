#' Subsample Minimum Spanning Trees of Specified Length
#'
#' Function that takes a pairwise matrix of great-circle distances and calculates a minimum-spanning tree of specified length
#' @param palaeocoordinates Matrix of Palaeocoordinates from PBDB occurrence file
#' @param gcd.matrix Pairwise matrix of Great Circle Distances; calculated if not supplied
#' @param desired.spread Length of desired subsampled MST
#' @param starting.point Use a random starting point, or choose the collection closest to the geographic centroid?
#' @param stochastic.sampling Grow the MST stochastically or deterministically? Defaults to FALSE.
#' @param check.max.spread Verify that maximum possible geographic spread exceeds the desired length. Defaults to TRUE.
#' @param error.margin Tolerance for accepting subsampled MSTs relative to desired spread; defaults to +/- 20%.
#' @keywords MST palaeocoordinates
#' @export
#' @examples
#' subsampleMST(palaeocoords, gcd.matrix, desired.spread = 5000, starting.point = "random", stochastic.sampling = FALSE, check.max.spread = FALSE, error.margin = 0.3)

subsampleMST <- function(palaeocoords, gcd.matrix = NULL, desired.spread = 5000, starting.point = c("random", "centroid"), stochastic.sampling = FALSE, check.max.spread = TRUE, error.margin = 0.2) {

	starting.point <- match.arg(starting.point, c("random", "centroid")) #match arguments

	okgo <- TRUE #do we run the analyses or not?

	sub.mst <- data.frame(distance = numeric(), from = numeric(), to = numeric())
	sub.pcoords <- data.frame(paleolatdec = numeric(), paleolngdec = numeric())
	totaldist <- numeric()
	median.palaeolat <- numeric()
	range.palaeolat <- numeric()
	iqr.palaeolat <- numeric()


	if (check.max.spread == TRUE) {

		cat("Checking that desired spread can be obtained from input data...\n")
		maximum.spread <- sum(dist2MST(pcoord2gcd(palaeocoords))$distance)
		if (desired.spread > maximum.spread) {
			cat("Longest possible MST (", maximum.spread, ") is smaller than desired subsampled distance!\n")
			okgo <- FALSE
		}

	}

	if (okgo == TRUE) {

		if (is.null(gcd.matrix)) {
			library(geosphere)
			gcd.matrix <- pcoord2gcd(palaeocoords) #calculate geographic distance matrix if not supplied
		}

		gcd.matrix.bak <- gcd.matrix #create backup of gcd matrix
		gcd.matrix <- as.matrix(stats::as.dist(t(gcd.matrix), diag = TRUE, upper = TRUE)) #convert gcd matrix to distance-matrix data type
		diag(gcd.matrix) <- NA #fill diagonal with NAs

		sub.mst <- matrix(nrow = 0, ncol = 3); colnames(sub.mst) <- c("distance", "from", "to") #create blank MST object

		if (starting.point == "random") { #choose a random collection as a starting point

			cat("Using random starting point.\n")

			origin <- sample(rownames(gcd.matrix), 1)

		}

		if (starting.point == "centroid") { #find the collection closest to centroid of palaeocoordinates

			cat("Using nearest collection to centroid as starting point\n")
			cat("Calculating geographic centroid...\n")
			sample.centroid <- t(as.matrix(centroid(palaeocoords[, c(2,1)])[, c(2,1)])) #calculate the geographic centroid
			colnames(sample.centroid) <- colnames(palaeocoords); row.names(sample.centroid) <- "centroid" #name them to conform to palaeocoords
			pcoords1 <- rbind(sample.centroid, palaeocoords) #bind the centroid to the palaeocoords
			centroid.distances <- pcoord2gcd(pcoords1) #calculate distance matrix for palaeocoords including centroid
			origin <- names(which.min(centroid.distances["centroid", ])) #find collection closest to centroid
			cat("Done.\n")

		}

		posscon <- gcd.matrix[rownames(gcd.matrix) %in% origin, !colnames(gcd.matrix) %in% origin] #get named vector of distances to unadded collections
		posscon <- posscon[!is.na(posscon)] #strip away NAs
		closest <- names(which.min(posscon)) #find the closest collection(s) to the origin
		closest.pcoords <- palaeocoords[closest,] #get pcoords of closest collection
		closest <- rownames(palaeocoords[palaeocoords$paleolatdec == closest.pcoords[, "paleolatdec"] & palaeocoords$paleolngdec == closest.pcoords[, "paleolngdec"], ]) #add all collections with identical pcoords
		colls <- unique(c(origin, closest)) #collate the first two collection names
		last.totaldist <- totaldist <- gcd.matrix[rownames(gcd.matrix) %in% colls[1], colnames(gcd.matrix) %in% colls[2]] #calculate distance

		i <- 1

		if (totaldist < desired.spread) {

			while (nrow(palaeocoords) > 2 && totaldist < desired.spread) {

				cat(paste("Calculating candidate MST (step ", i, "). ", sep = ""))

				posscon <- gcd.matrix[rownames(gcd.matrix) %in% colls, !colnames(gcd.matrix) %in% colls] #get possible next-closest collections

				if (is.null(nrow(posscon))) {

					cat("No more possible connections left; calculating final MST!\n")

					tmp.dist <- gcd.matrix[rownames(gcd.matrix) %in% colls, colnames(gcd.matrix) %in% colls]

					tmp.dist[lower.tri(tmp.dist)] <- NA

					sub.mst <- dist2MST(tmp.dist)

					totaldist <- sum(sub.mst$distance)

					break

				}

				if (!is.null(nrow(posscon)) && nrow(posscon) > 1) {posscon <- apply(posscon, 2, min)} #find minimum connection distance for each possible new connection
				posscon <- posscon[unique(names(posscon))] #strip out dupes due to symmetrisation of matrix
				posscon <- posscon[!is.na(posscon)] #remove NAs
				posscon <- posscon[order(posscon)]

				if (length(posscon) == 0 && totaldist < desired.spread) {stop("We ran out of places to add before the desired spread was reached! :'(\n")}

				if (length(posscon) == 1) {closest <- names(posscon)[1]}

				if (stochastic.sampling == FALSE && length(posscon) > 1) {

					closest <- names(which.min(posscon)) #get the name of the next closest connection
					closest.pcoords <- palaeocoords[closest,]
					closest <- rownames(palaeocoords[palaeocoords$paleolatdec == closest.pcoords[, "paleolatdec"] & palaeocoords$paleolngdec == closest.pcoords[, "paleolngdec"], ]) #add all collections with identical pcoords

				}

				if (stochastic.sampling == TRUE && length(posscon) > 1) { #this option chooses the next closest collection(s) to add with a probability proportional to their proximity

					sample.probs <- (1/posscon)/sum(1/posscon)
					closest <- names(sample(posscon, size = 1, prob = sample.probs))

				}

				last.totaldist <- totaldist
				totaldist <- sum(posscon[closest], totaldist, na.rm = TRUE)

				last.colls <- colls
				colls <- unique(c(colls, closest))

				if (totaldist >= desired.spread) {

					cat("Actually calculating MST this time... ")

					unique.colls <- rownames(unique(palaeocoords[colls, ]))

					tmp.dist <- gcd.matrix[rownames(gcd.matrix) %in% unique.colls, colnames(gcd.matrix) %in% unique.colls]

					tmp.dist[lower.tri(tmp.dist)] <- NA

					sub.mst <- dist2MST(tmp.dist)

					totaldist <- sum(sub.mst$distance)

					cat("Done.\n")

					if (totaldist >= desired.spread && length(unique.colls < length(colls))) {

						cat("Wahay, calculating final MST with all collections... ")

						tmp.dist <- gcd.matrix[rownames(gcd.matrix) %in% colls, colnames(gcd.matrix) %in% colls]

						tmp.dist[lower.tri(tmp.dist)] <- NA

						sub.mst <- dist2MST(tmp.dist)

						totaldist <- sum(sub.mst$distance)

						cat("Done.\n")

					}

				}

				i <- i + 1

				cat(paste("Total length = ", round(totaldist), " km \n", sep = ""))

			}

		}


		if (totaldist > (desired.spread + (desired.spread * error.margin)) && abs(totaldist - desired.spread) > abs(last.totaldist - desired.spread)) {

			sub.mst.bak <- sub.mst; totaldist.bak <- totaldist; colls.bak <- colls

			cat(paste("Final distance of ", round(totaldist), " km is at least 20% longer than desired distance of ", desired.spread, " km... calculating MST for previous step.\n", sep = ""))

			colls <- last.colls

			tmp.dist <- gcd.matrix[rownames(gcd.matrix) %in% colls, colnames(gcd.matrix) %in% colls]
			tmp.dist[lower.tri(tmp.dist)] <- NA

			last.sub.mst <- dist2MST(tmp.dist)

			last.totaldist <- sum(last.sub.mst$distance)

			if (abs(last.totaldist - desired.spread) < abs(totaldist - desired.spread)) {totaldist <- last.totaldist; sub.mst <- last.sub.mst}

		}

		sub.pcoords <- palaeocoords[rownames(palaeocoords) %in% colls, ]

		median.palaeolat <- medianPalaeolat(sub.pcoords)

		range.palaeolat <- abs(range(sub.pcoords$paleolatdec)[1] - range(sub.pcoords$paleolatdec)[2])

		iqr.palaeolat <- IQR(sub.pcoords$paleolatdec)

	}

	out <- setNames(
		list(sub.mst, sub.pcoords, totaldist, median.palaeolat, range.palaeolat, iqr.palaeolat),
		c("sub.mst", "sub.pcoords", "totaldist", "median.palaeolat", "range.palaeolat", "iqr.palaeolat")
	)

	return(out)

}


