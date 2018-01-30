#### DEFINE TIME INTERVALS ####
# time_int1 <- read.csv("https://paleobiodb.org/data1.2/intervals/list.txt?scale=1&limit=all") # 1.2 timescale containing more intervals
# write.csv(time_int1, file = "./input-data/pbdb_time_int.csv", row.names = FALSE)
time_int <- read.csv(file = "./input-data/pbdb_time_int.csv")

time_int <- filter(time_int, scale_no == 1) #only use Scale 1

if (use.stages == TRUE) {
	stages <- as.character(filter(time_int, scale_level == 5)$interval_name)
	bins <- setNames(as.list(stages), stages)

} else if (use.stages == FALSE ) {
	bins <- c(
			setNames(list(c("Calabrian", "Middle Pleistocene", "Late Pleistocene", "Holocene")), "Ng4"),
			setNames(list(c("Messinian", "Zanclean", "Piacenzian", "Gelasian")), "Ng3"),
			setNames(list(c("Langhian", "Serravallian", "Tortonian")), "Ng2"),
			setNames(list(c("Aquitanian", "Burdigalian")), "Ng1"),
			setNames(list(c("Chattian", "Rupelian")), "Pg5"),
			setNames(list(c("Bartonian", "Priabonian")), "Pg4"),
			setNames(list(c("Lutetian")), "Pg3"),
			setNames(list(c("Ypresian")), "Pg2"),
			setNames(list(c("Selandian", "Thanetian")), "Pg1"),
			setNames(list(c("Danian")), "Pg0"),
			setNames(list(c("Maastrichtian")), "K8"),
			setNames(list(c("Campanian")), "K7"),
			setNames(list(c("Coniacian", "Santonian", "Turonian")), "K6"),
			setNames(list(c("Cenomanian")), "K5"),
			setNames(list(c("Albian" )), "K4"),
			setNames(list(c("Aptian")), "K3"),
			setNames(list(c("Barremian", "Hauterivian")), "K2"),
			setNames(list(c("Berriasian", "Valanginian")), "K1"),
			setNames(list(c("Kimmeridgian", "Tithonian")), "J6"),
			setNames(list(c("Callovian", "Oxfordian")), "J5"),
			setNames(list(c("Bajocian", "Bathonian")), "J4"),
			setNames(list(c("Aalenian", "Toarcian")), "J3"),
			setNames(list(c("Pliensbachian")), "J2"),
			setNames(list(c("Hettangian", "Sinemurian")), "J1"),
			setNames(list(c("Rhaetian")), "Tr5"),
			setNames(list(c("Norian")), "Tr4"),
			setNames(list(c("Carnian")), "Tr3"),
			setNames(list(c("Ladinian")), "Tr2"),
			setNames(list(c("Anisian","Olenekian","Induan")), "Tr1"),
			setNames(list(c("Changhsingian", "Wuchiapingian")), "P5"),
			setNames(list(c("Capitanian", "Wordian")), "P4"),
			setNames(list(c("Kungurian", "Roadian")), "P3"),
			setNames(list(c("Artinskian")), "P2"),
			setNames(list(c("Asselian", "Sakmarian")), "P1"),
			setNames(list(c("Gzhelian", "Kasimovian")), "C5"),
			setNames(list(c("Moscovian")), "C4"),
			setNames(list(c("Bashkirian", "Serpukhovian")), "C3"),
			setNames(list(c("Visean")), "C2"),
			setNames(list(c("Tournaisian")), "C1"),
			setNames(list(c("Famennian")), "D5"),
			setNames(list(c("Frasnian")), "D4"),
			setNames(list(c("Eifelian", "Givetian")), "D3"),
			setNames(list(c("Emsian")), "D2"),
			setNames(list(c("Lochkovian", "Pragian")), "D1"),
			setNames(list(c("Gorstian", "Ludfordian", "Pridoli")), "S3"),
			setNames(list(c("Homerian", "Sheinwoodian", "Telychian")), "S2"),
			setNames(list(c("Aeronian", "Rhuddanian")), "S1"),
			setNames(list(c("Hirnantian", "Katian")), "Or5"),
			setNames(list(c("Sandbian")), "Or4"),
			setNames(list(c("Darriwilian")), "Or3"),
			setNames(list(c("Floian","Dapingian")), "Or2"),
			setNames(list(c("Tremadocian")), "Or1"),
			setNames(list(c("Paibian","Jiangshanian","Stage 10")), "Cm5"),
			setNames(list(c("Guzhangian","Drumian","Stage 5")), "Cm4"),
			setNames(list(c("Stage 4","Stage 3")), "Cm3"),
			setNames(list(c("Stage 2")), "Cm2"),
			setNames(list(c("Fortunian")), "Cm1")
		)
}


#construct intervals table from PBDB timescale data and bin definitions and write it to text file
pbdb_intervals <- as.data.frame(t(sapply(bins, function(x) range(filter(time_int, interval_name %in% x)[c("min_ma","max_ma")])))); colnames(pbdb_intervals) <- c("LAD","FAD") # extract FAD and LAD dates for bins from PBDB interval data
pbdb_intervals$stages <- do.call("rbind", lapply(bins, function(x) paste(x, collapse = ", ")))
pbdb_intervals$bin <- rownames(pbdb_intervals)
intervals <- pbdb_intervals[c("bin","stages","LAD","FAD")]; rownames(intervals) <- NULL
write.table(intervals[c("FAD","bin","stages")], file = "./input-data/timebins.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

intervals$midpoint <- (intervals$LAD + intervals$FAD) / 2 #calculate bin midpoints
intervals$duration <- intervals$FAD - intervals$LAD #calculate interval durations

lengths <- intervals$FAD - intervals$LAD #calculate lengths of intervals
names(lengths) <- intervals$bin #name the interval lengths

intervals.untrimmed <- intervals #make a backup of the untrimmed intervals

intervals <- intervals[which(intervals$bin %in% last.int):which(intervals$bin %in% first.int), ] #subset to the intervals of interest for this study

rownames(intervals) <- NULL
