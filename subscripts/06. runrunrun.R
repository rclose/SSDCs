#### RUN ANALYSES ####

#make sublists
if (run_SQS_Perl_1 == TRUE) {null.sqs.1[[rgn]] <- real.sqs.1[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))}
if (run_SQS_Perl_2 == TRUE) {null.sqs.2[[rgn]] <- real.sqs.2[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))}
if (run_iNEXT == TRUE) {null.inext[[rgn]] <- real.inext[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))}
if (run_lambda_5 == TRUE) {null.lambda5[[rgn]] <- real.lambda5[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))}
if (run_chao == TRUE) {null.chao[[rgn]] <- real.chao[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))}
if (run_CR == TRUE) {null.CR[[rgn]] <- real.CR[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))}
if (run_TRiPS == TRUE) {null.trips[[rgn]] <- real.trips[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))}
histoccRowIDs[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times))

goodsu[[rgn]] <- matrix(nrow = nrow(intervals), ncol = length(times), dimnames = list(intervals$bin, names(times)))

multiton_count[[rgn]] <- singleton_count[[rgn]] <- newoccs[[rgn]] <- newtax[[rgn]] <- numoccs[[rgn]] <- list()
for (int in 1:nrow(intervals)) { #for counting noveltons/oldtons
	singleton_count[[rgn]][[int]] <- multiton_count[[rgn]][[int]] <- newoccs[[rgn]][[int]] <- newtax[[rgn]][[int]] <- numoccs[[rgn]][[int]] <- setNames(rep(NA, length = length(times)), names(times))
}; names(multiton_count[[rgn]]) <- names(singleton_count[[rgn]]) <- names(newtax[[rgn]]) <- names(newoccs[[rgn]]) <- names(numoccs[[rgn]]) <- intervals$bin

if (run_PGS_calcs == TRUE) {
	pgs[[rgn]] <- setNames(vector(mode = "list", length = length(times)), names(times)) #palaeogeographic spread
	for (yr in 1:length(times)) {
		pgs[[rgn]][[yr]] <- setNames(vector(mode = "list", length = nrow(intervals)), intervals$bin)
	}
}


#### chronologically time-slice the data ####
for (yr in 1:length(times)) {

	#subset data at selected publication time-points
	histoccs <- subset(occs[[rgn]], ref_pubyr <= times[yr])

	#drop obsolete occurrence identifications
	cat(paste("Dropping obsolete occurrence identifications for ", clade, " ", names(times)[yr], "\n", sep = ""))
	if (any(duplicated(histoccs$occurrence_no))) {
		histoccs <- dropObsoleteOccs(histoccs)
	}

	#print an error if there are any obsolete occurrences remaining
	if (any(duplicated(histoccs$occurrence_no)) == TRUE) {stop("Some duplicate occurrence records remain")}

	#record the row IDs for later subsetting
	histoccRowIDs[[rgn]][[yr]] <- histoccs$rowID

}

#### tally enton counts, numbers of occurrences, etc
for (yr in 1:length(times)) {

	histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])

	#### tally occurrences for each publication year per interval ####
	if (rank == "species") {
		tmptable <- sort(table(histoccs$occurrence.binomial))
	} else if (rank == "genus") {
		tmptable <- sort(table(histoccs$occurrence.genus_name))
	}

	for (int in 1:nrow(intervals)) {
		#subset occurrences for interval
		intoccs <- histoccs[histoccs$ma_min >= intervals[int,"LAD"] & histoccs$ma_max <= intervals[int,"FAD"], ]

		#tally total occurrences known for each publication year per interval
		numoccs[[rgn]][[int]][yr] <- setNames(nrow(intoccs), names(times)[yr])

		latest_occs <- subset(intoccs, ref_pubyr > times[yr + 1])

		#tally new/old taxa described each year
		if (rank == "species") {
			newtax[[rgn]][[int]][yr] <- setNames(nrow(latest_occs[latest_occs$occurrence.binomial %in% names(tmptable[tmptable == 1]), ]), names(times)[yr])
			newoccs[[rgn]][[int]][yr] <- setNames(nrow(latest_occs[latest_occs$occurrence.binomial %in% names(tmptable[tmptable > 1]), ]), names(times)[yr])
			singleton_count[[rgn]][[int]][yr] <- sum(table(intoccs$occurrence.binomial) == 1)
			multiton_count[[rgn]][[int]][yr] <- sum(table(intoccs$occurrence.binomial) > 1)
		} else if (rank == "genus") {
			newtax[[rgn]][[int]][yr] <- setNames(nrow(latest_occs[latest_occs$occurrence.genus_name %in% names(tmptable[tmptable == 1]), ]), names(times)[yr])
			newoccs[[rgn]][[int]][yr] <- setNames(nrow(latest_occs[latest_occs$occurrence.genus_name %in% names(tmptable[tmptable > 1]), ]), names(times)[yr])
			singleton_count[[rgn]][[int]][yr] <- sum(table(intoccs$occurrence.genus_name) == 1)
			multiton_count[[rgn]][[int]][yr] <- sum(table(intoccs$occurrence.genus_name) > 1)
		}

	}; names(numoccs[[rgn]]) <- names(newtax[[rgn]]) <- names(newoccs[[rgn]]) <- intervals$bin

}

#### calculate palaeogeographic spread for interval ####
if (run_PGS_calcs == T) {

	for (yr in 1:length(times)) {

		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])

		#### extract NX2 matrices of palaeocoordinates from PBDB data for all regions/intervals ####
		pcoords <- getPalaeoCoords(data = histoccs[!is.na(histoccs$paleolatdec), ], intervals = intervals, unique.pcoords = TRUE, occs.entirely.within.bin = TRUE)

		#### calculate palaeogeographic spreads ####
		for (int in 1:nrow(intervals)) {

			cat(paste("Calculating palaeogeographic spreads for region ", names(regions)[rgn], ", interval ", intervals[int,"bin"],"\n", sep = ""))
			pgs[[rgn]][[yr]][[int]] <- PGS(pcoords = pcoords[[int]], round.mst.pcoords = 2, grid.cell.size = 2)

		}; pgs[[rgn]][[yr]] <- do.call('rbind', pgs[[rgn]][[yr]]); rownames(pgs[[rgn]][[yr]]) <- intervals$bin

	}
}

#### run SQS Perl script version 1 on each historical time-point ####
if (run_SQS_Perl_1 == TRUE) {
	for (yr in 1:length(times)) {
		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
		real.sqs.1[[rgn]][[yr]] <- list()
		for (q in 1:length(quorum.levels)) {
			cat("Running SQS with quota of 3 collections per reference at quorum", quorum.levels[q], "for ", names(regions)[rgn], " occurrences published before", times[yr], "\n", sep = " ")
			out <- runSQS(data = histoccs,
				      quorum.level = quorum.levels[q], # quorum level
				      ntrials = ntrials, # number of subsampling trials
				      rank = rank, # taxonomic rank (required, options genus or species)
				      method = "SQS", # options CR, UW, or O2W (default SQS)
				      exact = "yes", # draw all occurrences and count taxa seen at the quorum level(s) (option no, default yes)
				      bycollection = "yes", # draw entire collections (option no, default yes)
				      singletons = "occ", # forced to use single-occurrence definition of singletons when EXACT option is used
				      exclude = "no", # exclude most common (dominant) genus (option yes, default no)
				      biggest = "no", # exclude taxa in most diverse collection from singleton count (option yes, default no)
				      sequential = "yes", # draw all collections in a reference before going to the next one (optional and recommended)
				      collsperref = 3, # fixed number of collections drawn per reference (default no limit, but recommended)
				      refquota = 0, # minimum number of references per bin (default none, and not recommended)
				      disperse = "no", #disabled anyway when EXACT option is used
				      interpolate = "no", # use interpolated boundary estimates in PaleoDB download file (default no, but recommended by Alroy)
				      usefailed = "no", # print data for bins sometimes under quota (default no)
				      trim = "no", # print data only for bins in the range of those that meet the subsampling target (default no)
				      mincollbybin = 1, # number of collections a bin must include to be analyzed
				      matchmax = "no", # intervals must include as many collections as the maximum drawn (default no, and not recommended)
				      deorphan = "no" # assign collections spanning multiple bins to the ones including more than half of their age estimate limits (default no; options yes or a fraction)
			)
			out <- out[out$Bin.name %in% intervals$bin, ] #trim output to interval of interest
			real.sqs.1[[rgn]][[yr]][[q]] <- out
			rownames(real.sqs.1[[rgn]][[yr]][[q]]) <- real.sqs.1[[rgn]][[yr]][[q]]$Bin.name
			real.sqs.1[[rgn]][[yr]][[q]]$Bin.name <- NULL
		}; names(real.sqs.1[[rgn]][[yr]]) <- names(quorum.levels); real.sqs.1[[rgn]][[yr]] <- abind(real.sqs.1[[rgn]][[yr]], along = 3)
	}
}

#### Run SQS Perl script version 2 ####
if (run_SQS_Perl_2 == TRUE) {
	for (yr in 1:length(times)) {
		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
		real.sqs.2[[rgn]][[yr]] <- list()
		for (q in 1:length(quorum.levels)) {
			cat("Running SQS with no quota at quorum", quorum.levels[q], "for ", names(regions)[rgn], " occurrences published before", times[yr], "\n", sep = " ")
			out <- runSQS(data = histoccs,
				      quorum.level = quorum.levels[q], # quorum level
				      ntrials = ntrials, # number of subsampling trials
				      rank = rank, # taxonomic rank (required, options genus or species)
				      method = "SQS", # options CR, UW, or O2W (default SQS)
				      exact = "yes", # draw all occurrences and count taxa seen at the quorum level(s) (option no, default yes)
				      bycollection = "no", # draw entire collections (option no, default yes)
				      singletons = "occ", # definition of singletons (option one reference and not recommended, default one occurrence); forced to use single-occurrence definition of singletons when EXACT option is used
				      exclude = "no", # exclude most common (dominant) genus (option yes, default no)
				      biggest = "no", # exclude taxa in most diverse collection from singleton count (option yes, default no)
				      sequential = "no", # draw all collections in a reference before going to the next one (optional and recommended)
				      collsperref = 0, # fixed number of collections drawn per reference (default no limit, but recommended)
				      refquota = 0, # minimum number of references per bin (default none, and not recommended)
				      disperse = "no", # disperse sampling among references with a throwback algorithm (optional and not recommended); disabled when EXACT option is used
				      interpolate = "no", # use interpolated boundary estimates in PaleoDB download file (default no, but recommended by Alroy)
				      usefailed = "no", # print data for bins sometimes under quota (default no)
				      trim = "no", # print data only for bins in the range of those that meet the subsampling target (default no)
				      mincollbybin = 1, # number of collections a bin must include to be analyzed
				      matchmax = "no", # intervals must include as many collections as the maximum drawn (default no, and not recommended)
				      deorphan = "no" # assign collections spanning multiple bins to the ones including more than half of their age estimate limits (default no; options yes or a fraction)
			)
			out <- out[out$Bin.name %in% intervals$bin, ] #trim output to interval of interest
			real.sqs.2[[rgn]][[yr]][[q]] <- out
			rownames(real.sqs.2[[rgn]][[yr]][[q]]) <- real.sqs.2[[rgn]][[yr]][[q]]$Bin.name
			real.sqs.2[[rgn]][[yr]][[q]]$Bin.name <- NULL
		}; names(real.sqs.2[[rgn]][[yr]]) <- names(quorum.levels); real.sqs.2[[rgn]][[yr]] <- abind(real.sqs.2[[rgn]][[yr]], along = 3)
	}
}



#### run iNEXT ####
if (run_iNEXT == TRUE) {
	for (yr in 1:length(times)) {
		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
		real.inext[[rgn]][[yr]] <- list()
		for (int in 1:nrow(intervals)) {
			real.inext[[rgn]][[yr]][[int]] <- list()
			dt <- intervals[int,]$FAD - intervals[int,]$LAD
			int.occs <- histoccs[histoccs$ma_min >= intervals[int,"LAD"] & histoccs$ma_max <= intervals[int,"FAD"], ]
			if (rank == "species") {
				int.occs <- int.occs[!grepl("sp\\.", int.occs$occurrence.species_name), ]
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
			} else if (rank == "genus") {
				obs <- sort(table(int.occs$occurrence.genus_name), decreasing = TRUE) #occurrence counts per genus
			}
			obs <- c(obs, use.names = TRUE) #convert to vector, because it's an array
			n_raw <- length(obs) # Observed richness
			for (q in 1:length(quorum.levels)) {
				cat(paste("iNEXT for ref_pubyr =", times[yr], "interval =", intervals$bin[int], "quorum =", quorum.levels[q], "\n", sep = " "))
				if (exists("out")) {rm(list = "out")} #iNEXT gets confused if there's an object called "out" that already exists
				if (n_raw > 2 && any(obs > 1)) {
					out.iNEXT <- estimateD(x = unname(obs), datatype = "abundance", level = quorum.levels[q], base = "coverage") #run analytical coverage-based interpolation/extrapolation with confidence intervals
				} else if (n_raw <= 2 || all(obs == 1) || nrow(out.iNEXT) > 3) {
					cat("iNEXT output malformed because of quorum issue; substituting blank output \n")
					out.iNEXT <- as.data.frame(matrix(nrow = 3, ncol = 7, dimnames = list(c(1,2,3), c("m","method","order","SC","qD","qD.LCL","qD.UCL"))))
					out.iNEXT[,"order"] <- c(0,1,2)
					out.iNEXT[,"m"] <- c(0,0,0)
					out.iNEXT[,"SC"] <- rep(unname(quorum.levels[q]), 3)
					out.iNEXT[,"method"] <- rep("interpolated", 3)
				}
				out.iNEXT$quorum <- rep(names(quorum.levels)[q], 3)
				real.inext[[rgn]][[yr]][[int]][[q]] <- out.iNEXT
			}; real.inext[[rgn]][[yr]][[int]] <- do.call("rbind", real.inext[[rgn]][[yr]][[int]])
		}; names(real.inext[[rgn]][[yr]]) <- intervals$bin
	}
}



#### run iNEXT CR ####
if (run_CR == TRUE) {
	for (yr in 1:length(times)) {
		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
		real.CR[[rgn]][[yr]][[int]] <- list()
		for (int in 1:nrow(intervals)) {
			dt <- intervals[int,]$FAD - intervals[int,]$LAD
			int.occs <- histoccs[histoccs$ma_min >= intervals[int,"LAD"] & histoccs$ma_max <= intervals[int,"FAD"], ]
			if (rank == "species") {
				int.occs <- int.occs[!grepl("sp\\.", int.occs$occurrence.species_name), ]
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
			} else if (rank == "genus") {
				obs <- sort(table(int.occs$occurrence.genus_name), decreasing = TRUE) #occurrence counts per genus
			}
			obs <- c(obs, use.names = TRUE) #convert to vector, because it's an array for some reason
			n_raw <- length(obs) # Observed richness
			for (q in 1:length(quota.levels)) {
				cat(paste("iNEXT CR for ref_pubyr =", times[yr], "interval =", intervals$bin[int], "quota =", quota.levels[q], "\n", sep = " "))
				if (exists("out")) {rm(list = "out")} #iNEXT gets confused if there's an object called "out" that already exists
				if (n_raw > 2 && any(obs > 1)) {
					out.CR <- estimateD(x = unname(obs), datatype = "abundance", level = quota.levels[q], base = "size") #run analytical coverage-based interpolation/extrapolation with confidence intervals
				} else if (n_raw <= 2 || all(obs == 1) || nrow(out.iNEXT) > 3) {
					cat("iNEXT output malformed because of quota issue; substituting blank output :'( \n")
					out.CR <- as.data.frame(matrix(nrow = 3, ncol = 7, dimnames = list(c(1,2,3), c("m","method","order","SC","qD","qD.LCL","qD.UCL"))))
					out.CR[,"order"] <- c(0,1,2)
					out.CR[,"m"] <- c(0,0,0)
					out.CR[,"SC"] <- rep(NA, 3)
					out.CR[,"method"] <- rep("interpolated", 3)
				}
				out.CR$quota <- rep(names(quota.levels)[q], 3)
				real.CR[[rgn]][[yr]][[int]][[q]] <- out.CR
			}; real.CR[[rgn]][[yr]][[int]] <- do.call("rbind", real.CR[[rgn]][[yr]][[int]])
		}; names(real.CR[[rgn]][[yr]]) <- intervals$bin
	}
}

#### run lambda-5 ####
if (run_lambda_5 == TRUE) {
	for (yr in 1:length(times)) {
		cat(paste("Running lambda-5 for ref_pubyr =", times[yr], "...\n", sep = " "))
		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
		real.lambda5[[rgn]][[yr]] <- list()
		real.lambda5[[rgn]][[yr]] <- mclapply(1:nrow(intervals), function(int) {
			int.occs <- histoccs[histoccs$ma_min >= intervals[int,"LAD"] & histoccs$ma_max <= intervals[int,"FAD"], ]
			if (rank == "species") {
				int.occs <- int.occs[!grepl("sp\\.", int.occs$occurrence.species_name), ]
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
			} else if (rank == "genus") {
				obs <- sort(table(int.occs$occurrence.genus_name), decreasing = TRUE) #occurrence counts per genus
			}
			obs <- c(obs, use.names = TRUE) #convert to vector, because it's an array for some reason
			out <- bootstrapLambda(obs)
			names(out) <- c("lambda5.LCL","lambda5","lambda5.UCL")
			out
		}, mc.cores = detectCores(), mc.preschedule = T)
		real.lambda5[[rgn]][[yr]] <- do.call("rbind", real.lambda5[[rgn]][[yr]]); rownames(real.lambda5[[rgn]][[yr]]) <- intervals$bin
	}
}


#### run Chao ####
if (run_chao == TRUE) {
	for (yr in 1:length(times)) {
		cat(paste("Running Chao for ref_pubyr =", times[yr], "...\n", sep = " "))
		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
		real.chao[[rgn]][[yr]] <- list()
		for (int in 1:nrow(intervals)) {
			dt <- intervals[int,]$FAD - intervals[int,]$LAD
			int.occs <- histoccs[histoccs$ma_min >= intervals[int,"LAD"] & histoccs$ma_max <= intervals[int,"FAD"], ]
			if (rank == "species") {
				int.occs <- int.occs[!grepl("sp\\.", int.occs$occurrence.species_name), ]
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
			} else if (rank == "genus") {
				obs <- sort(table(int.occs$occurrence.genus_name), decreasing = TRUE) #occurrence counts per genus
			}
			obs <- c(obs, use.names = TRUE) #convert to vector, because it's an array for some reason
			n_raw <- length(obs) # Observed richness
			real.chao[[rgn]][[yr]][[int]] <- bootstrapChao(obs)
			names(real.chao[[rgn]][[yr]][[int]]) <- c("chao.LCL","chao","chao.UCL")
		}; real.chao[[rgn]][[yr]] <- do.call("rbind", real.chao[[rgn]][[yr]]); rownames(real.chao[[rgn]][[yr]]) <- intervals$bin
	}
}


#### calculate sample coverage (Chao and Jost (2012)/Chao and Shen (2010) formula ####
for (yr in 1:length(times)) {
	histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
	for (int in 1:nrow(intervals)) {
		occvec <- histoccs %>% filter(ma_max <= intervals[int, "FAD"] & ma_min >= intervals[int, "LAD"]) %>% select(occurrence.binomial) %>% table
		n <- sum(occvec)
		f1 <- sum(occvec == 1); f2 <- sum(occvec == 2)
		# out <- 1 - (f1 / sum(occvec)) #singleton-only coverage estimator
		out <- 1 - (f1/n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2)) #Chao and Jost (2012) coverage estimator using singletons and doubletons
		out[is.nan(out)] <- NA
		out[is.infinite(out)] <- NA
		out[out < 0.01] <- NA
		goodsu[[rgn]][int, yr] <- out
	}
}


#### run TRiPS ####
if (run_TRiPS == TRUE) {
	for (yr in 1:length(times)) {
		histoccs <- filter(occs[[rgn]], rowID %in% histoccRowIDs[[rgn]][[yr]])
		real.trips[[rgn]][[yr]] <- list()
		for (int in 1:nrow(intervals)) {
			dt <- intervals[int,]$FAD - intervals[int,]$LAD
			cat(paste("TRiPS for ref_pubyr", times[yr], "interval", intervals$bin[int], "\n", sep = " "))
			int.occs <- histoccs[histoccs$ma_min >= intervals[int,"LAD"] & histoccs$ma_max <= intervals[int,"FAD"], ]
			if (rank == "species") {
				int.occs <- int.occs[!grepl("sp\\.", int.occs$occurrence.species_name), ]
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
			} else if (rank == "genus") {
				obs <- sort(table(int.occs$occurrence.genus_name), decreasing = TRUE) #occurrence counts per genus
			}
			obs <- c(obs, use.names = TRUE) #convert to vector, because it's an array
			n_raw <- length(obs) # Observed richness
			multitonprop <- 1 - (length(obs[obs == 1])/length(obs)) # Good's u (proportion of non-singletons), used to determine if TRiPS should be run.
			if (any(obs > 1) && multitonprop > 0.1 && length(obs) > 1) {
				if (exists("out")) {rm(list = "out")}
				real.trips[[rgn]][[yr]][[int]] <- doTRiPS_abs(obs, t = dt)
			} else {
				real.trips[[rgn]][[yr]][[int]] <- matrix(nrow = 3, ncol = 3, dimnames = list(c("Sampling rate","Sampling probability","Estimated richness"),c("MLE","lower CI","upper CI")))
			}
		}; names(real.trips[[rgn]][[yr]]) <- intervals$bin; real.trips[[rgn]][[yr]] <- abind(real.trips[[rgn]][[yr]], along = 3)
	}
}


if (run_PGS_calcs == TRUE) {pgs[[rgn]] <- abind(pgs[[rgn]], along = 3)}
if (run_SQS_Perl_1 == TRUE) {real.sqs.1[[rgn]] <- abind(real.sqs.1[[rgn]], along = 4)}
if (run_SQS_Perl_2 == TRUE) {real.sqs.2[[rgn]] <- abind(real.sqs.2[[rgn]], along = 4)}
if (run_TRiPS == TRUE) {real.trips[[rgn]] <- abind(real.trips[[rgn]], along = 4)}


#### generate null distribution for TRiPS ####
if (nulldist == TRUE && run_TRiPS == TRUE) {
	null.trips[[rgn]] <- list()
	for (i in 1:nshuff) {
		shuffled <- occs[[rgn]]
		shuffled$ref_pubyr <- sample(shuffled$ref_pubyr, replace = FALSE)
		null.trips[[rgn]][[i]] <- list()
		for (j in 1:length(times)) {
			cat("TRiPS null rep", i, " for occurrences published before", times[j], "\n", sep = " ")
			tmpoccs <- subset(shuffled, ref_pubyr <= times[j])
			null.trips[[rgn]][[i]][[j]] <- list()
			for (k in 1:nrow(intervals)) {
				dt <- intervals[k,]$FAD - intervals[k,]$LAD
				int.occs <- tmpoccs[tmpoccs$ma_min >= intervals[k,"LAD"] & tmpoccs$ma_max <= intervals[k,"FAD"], ] # & !is.na(match(all_data$country,regions[[i]]))
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
				n_raw <- length(obs) # Observed richness
				if (any(obs > 1)) {
					null.trips[[rgn]][[i]][[j]][[k]] <- doTRiPS_abs(obs, t = dt)
				} else {null.trips[[rgn]][[i]][[j]][[k]] <- matrix(nrow = 3, ncol = 3, dimnames = list(c("Sampling rate","Sampling probability","Estimated richness"),c("MLE","lower CI","upper CI")))}
			}; names(null.trips[[rgn]][[i]][[j]]) <- intervals$bin; null.trips[[rgn]][[i]][[j]] <- abind(null.trips[[rgn]][[i]][[j]], along = 3)
			rm(tmpoccs)
		}; names(null.trips[[rgn]][[i]]) <- paste("ref_pubyr", times, sep = ""); null.trips[[rgn]][[i]] <- abind(null.trips[[rgn]][[i]], along = 4)
	}; names(null.trips[[rgn]]) <- paste("rep", 1:nshuff, sep = ""); null.trips[[rgn]] <- abind(null.trips[[rgn]], along = 5)
}


#### generate null distribution for lambda-5 ####
if (nulldist == TRUE && run_lambda_5 == TRUE) {
	null.lambda5[[rgn]] <- list()
	for (i in 1:nshuff) {
		shuffled <- occs[[rgn]]
		shuffled$ref_pubyr <- sample(shuffled$ref_pubyr, replace = FALSE)
		null.lambda5[[rgn]][[i]] <- list()
		for (j in 1:length(times)) {
			cat("lambda-5 null rep", i, " for occurrences published before", times[j], "\n", sep = " ")
			tmpoccs <- subset(shuffled, ref_pubyr <= times[j])
			null.lambda5[[rgn]][[i]][[j]] <- list()
			for (k in 1:nrow(intervals)) {
				int.occs <- tmpoccs[tmpoccs$ma_min >= intervals[k,"LAD"] & tmpoccs$ma_max <= intervals[k,"FAD"], ]
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
				if (length(obs >= 1)) {
					null.lambda5[[rgn]][[i]][[j]][[k]] <- lambda5(obs)
				} else {null.lambda5[[rgn]][[i]][[j]][[k]] <- NA}
			}; names(null.lambda5[[rgn]][[i]][[j]]) <- intervals$bin; null.lambda5[[rgn]][[i]][[j]] <- abind(null.lambda5[[rgn]][[i]][[j]], along = 2)
			rm(tmpoccs)
		}; names(null.lambda5[[rgn]][[i]]) <- paste("ref_pubyr", times, sep = ""); null.lambda5[[rgn]][[i]] <- abind(null.lambda5[[rgn]][[i]], along = 3)
	}; names(null.lambda5[[rgn]]) <- paste("rep", 1:nshuff, sep = ""); null.lambda5[[rgn]] <- abind(null.lambda5[[rgn]], along = 4)
}


#### generate null distribution for Chao1 ####
if (nulldist == TRUE && run_chao == TRUE) {
	null.chao[[rgn]] <- list()
	for (i in 1:nshuff) {
		shuffled <- occs[[rgn]]
		shuffled$ref_pubyr <- sample(shuffled$ref_pubyr, replace = FALSE)
		null.chao[[rgn]][[i]] <- list()
		for (j in 1:length(times)) {
			cat("Chao1 null rep", i, " for occurrences published before", times[j], "\n", sep = " ")
			tmpoccs <- subset(shuffled, ref_pubyr <= times[j])
			null.chao[[rgn]][[i]][[j]] <- list()
			for (k in 1:nrow(intervals)) {
				int.occs <- tmpoccs[tmpoccs$ma_min >= intervals[k,"LAD"] & tmpoccs$ma_max <= intervals[k,"FAD"], ] # & !is.na(match(all_data$country,regions[[i]]))
				obs <- sort(table(int.occs$occurrence.binomial), decreasing = TRUE) #occurrence counts per species
				if (length(obs >= 1)) {
					null.chao[[rgn]][[i]][[j]][[k]] <- Chao1(obs)
				} else {null.chao[[rgn]][[i]][[j]][[k]] <- NA}
			}; names(null.chao[[rgn]][[i]][[j]]) <- intervals$bin; null.chao[[rgn]][[i]][[j]] <- abind(null.chao[[rgn]][[i]][[j]], along = 2)
			rm(tmpoccs)
		}; names(null.chao[[rgn]][[i]]) <- paste("ref_pubyr", times, sep = ""); null.chao[[rgn]][[i]] <- abind(null.chao[[rgn]][[i]], along = 3)
	}; names(null.chao[[rgn]]) <- paste("rep", 1:nshuff, sep = ""); null.chao[[rgn]] <- abind(null.chao[[rgn]], along = 4)
}


#### generate null distribution for SQS Perl Script Version 1 ####
if (nulldist == TRUE && run_SQS_Perl_1 == TRUE) {
	null.sqs.1[[rgn]] <- list()
	for (i in 1:nshuff) {
		shuffled <- occs[[rgn]]
		shuffled$ref_pubyr <- sample(shuffled$ref_pubyr, replace = FALSE)
		null.sqs.1[[rgn]][[i]] <- list()
		for (j in 1:length(times)) {
			cat("SQS Perl Script Version 1 null rep", i, "for occurrences published before", times[j], "\n", sep = " ")
			tmpoccs <- subset(shuffled, ref_pubyr <= times[j])
			null.sqs.1[[rgn]][[i]][[j]] <- list()
			for (k in 1:length(quorum.levels)) {
				out <- runSQS(data = tmpoccs,
					      quorum.level = quorum.levels[k], # quorum level
					      ntrials = nulltrials, # number of subsampling trials
					      rank = rank, # taxonomic rank (required, options genus or species)
					      method = "SQS", # options CR, UW, or O2W (default SQS)
					      exact = "yes", # draw all occurrences and count taxa seen at the quorum level(s) (option no, default yes)
					      bycollection = "yes", # draw entire collections (option no, default yes)
					      singletons = "occ", # forced to use single-occurrence definition of singletons when EXACT option is used
					      exclude = "no", # exclude most common (dominant) genus (option yes, default no)
					      biggest = "no", # exclude taxa in most diverse collection from singleton count (option yes, default no)
					      sequential = "yes", # draw all collections in a reference before going to the next one (optional and recommended)
					      collsperref = 3, # fixed number of collections drawn per reference (default no limit, but recommended)
					      refquota = 0, # minimum number of references per bin (default none, and not recommended)
					      disperse = "no", #disabled anyway when EXACT option is used
					      interpolate = "no", # use interpolated boundary estimates in PaleoDB download file (default no, but recommended by Alroy--except interpolated boundaries aren't present in PaleoBioDB API downloads, so useless)
					      usefailed = "no", # print data for bins sometimes under quota (default no)
					      trim = "no", # print data only for bins in the range of those that meet the subsampling target (default no)
					      mincollbybin = 1, # number of collections a bin must include to be analyzed
					      matchmax = "no", # intervals must include as many collections as the maximum drawn (default no, and not recommended)
					      deorphan = "no", # assign collections spanning multiple bins to the ones including more than half of their age estimate limits (default no; options yes or a fraction)
					      ignore.stderr = T, #don't print errors to the console
					      ignore.stdout = T #don't print standard output to console (may stop function slowing down over time?)
				)
				out <- out[out$Bin.name %in% intervals$bin, ] #trim output to interval of interest
				null.sqs.1[[rgn]][[i]][[j]][[k]] <- out
				rownames(null.sqs.1[[rgn]][[i]][[j]][[k]]) <- null.sqs.1[[rgn]][[i]][[j]][[k]]$Bin.name
				null.sqs.1[[rgn]][[i]][[j]][[k]]$Bin.name <- NULL
			}; names(null.sqs.1[[rgn]][[i]][[j]]) <- names(quorum.levels); null.sqs.1[[rgn]][[i]][[j]] <- abind(null.sqs.1[[rgn]][[i]][[j]], along = 3)
			rm(tmpoccs)
		}; names(null.sqs.1[[rgn]][[i]]) <- paste("ref_pubyr", times, sep = ""); null.sqs.1[[rgn]][[i]] <- abind(null.sqs.1[[rgn]][[i]], along = 4)
	}; names(null.sqs.1[[rgn]]) <- paste("rep", 1:nshuff, sep = ""); null.sqs.1[[rgn]] <- abind(null.sqs.1[[rgn]], along = 5)
}


#### generate null distribution for SQS Perl Script Version 2 ####
if (nulldist == TRUE && run_SQS_Perl_2 == TRUE) {
	null.sqs.2[[rgn]] <- list()
	for (i in 1:nshuff) {
		shuffled <- occs[[rgn]]
		shuffled$ref_pubyr <- sample(shuffled$ref_pubyr, replace = FALSE)
		null.sqs.2[[rgn]][[i]] <- list()
		for (j in 1:length(times)) {
			cat("SQS Perl Script Version 2 null rep", i, "for occurrences published before", times[j], "\n", sep = " ")
			tmpoccs <- subset(shuffled, ref_pubyr <= times[j])
			null.sqs.2[[rgn]][[i]][[j]] <- list()
			for (k in 1:length(quorum.levels)) {
				out <- runSQS(data = tmpoccs,
					      quorum.level = quorum.levels[k], # quorum level
					      ntrials = nulltrials, # number of subsampling trials
					      rank = rank, # taxonomic rank (required, options genus or species)
					      method = "SQS", # options CR, UW, or O2W (default SQS)
					      exact = "yes", # draw all occurrences and count taxa seen at the quorum level(s) (option no, default yes)
					      bycollection = "no", # draw entire collections (option no, default yes)
					      singletons = "occ", # forced to use single-occurrence definition of singletons when EXACT option is used
					      exclude = "no", # exclude most common (dominant) genus (option yes, default no)
					      biggest = "no", # exclude taxa in most diverse collection from singleton count (option yes, default no)
					      sequential = "yes", # draw all collections in a reference before going to the next one (optional and recommended)
					      collsperref = 0, # fixed number of collections drawn per reference (default no limit, but recommended)
					      refquota = 0, # minimum number of references per bin (default none, and not recommended)
					      disperse = "no", #disabled anyway when EXACT option is used
					      interpolate = "no", # use interpolated boundary estimates in PaleoDB download file (default no, but recommended)
					      usefailed = "no", # print data for bins sometimes under quota (default no)
					      trim = "no", # print data only for bins in the range of those that meet the subsampling target (default no)
					      mincollbybin = 1, # number of collections a bin must include to be analyzed
					      matchmax = "no", # intervals must include as many collections as the maximum drawn (default no, and not recommended)
					      deorphan = "no", # assign collections spanning multiple bins to the ones including more than half of their age estimate limits (default no; options yes or a fraction)
					      ignore.stderr = T, #don't print errors to the console (may stop function slowing down over time?)
					      ignore.stdout = T #don't print standard output to console (may stop function slowing down over time?)
				)
				out <- out[out$Bin.name %in% intervals$bin, ] #trim output to interval of interest
				null.sqs.2[[rgn]][[i]][[j]][[k]] <- out
				rownames(null.sqs.2[[rgn]][[i]][[j]][[k]]) <- null.sqs.2[[rgn]][[i]][[j]][[k]]$Bin.name
				null.sqs.2[[rgn]][[i]][[j]][[k]]$Bin.name <- NULL
			}; names(null.sqs.2[[rgn]][[i]][[j]]) <- names(quorum.levels); null.sqs.2[[rgn]][[i]][[j]] <- abind(null.sqs.2[[rgn]][[i]][[j]], along = 3)
			rm(tmpoccs)
		}; names(null.sqs.2[[rgn]][[i]]) <- paste("ref_pubyr", times, sep = ""); null.sqs.2[[rgn]][[i]] <- abind(null.sqs.2[[rgn]][[i]], along = 4)
	}; names(null.sqs.2[[rgn]]) <- paste("rep", 1:nshuff, sep = ""); null.sqs.2[[rgn]] <- abind(null.sqs.2[[rgn]], along = 5)
}
