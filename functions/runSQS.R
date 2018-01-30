#### runSQS: run SQS Perl script ####
runSQS <- function(
	data = data,
	collections = NULL,
	rank = c("species","genus"),
	method = c("SQS","CR","UW","O2W"),
	quorum.level = 0.4,
	ntrials = 1000,
	scale = "./input-data/timebins.txt",
	maxma = "650",
	exact = c("yes","no"),
	abund = c("no","yes"),
	bycollection = c("yes","no"),
	singletons = c("occ","ref"),
	exclude = c("no","yes"),
	biggest = c("no","yes"),
	sequential = c("yes","no"),
	collsperref = 3,
	refquota = 0,
	disperse = c("no","yes"),
	interpolate = c("no","yes"),
	usefailed = c("no","yes"),
	trim = c("no","yes"),
	mincollbybin = 1,
	matchmax = c("no","yes"),
	deorphan = c("no","yes"),
	ignore.stdout = F,
	ignore.stderr = F
	) {

	#match arguments
	rank <- match.arg(rank, c("species", "genus"))
	method <- match.arg(method, c("SQS","CR","UW","O2W"))
	exact <- match.arg(exact, c("yes","no"))
	abund <- match.arg(exact, c("no","yes"))
	exact <- match.arg(exact, c("yes","no"))
	bycollection <- match.arg(bycollection, c("yes","no"))
	singletons <- match.arg(singletons, c("occ","ref"))
	exclude <- match.arg(exclude, c("no","yes"))
	biggest <- match.arg(biggest, c("no","yes"))
	sequential <- match.arg(sequential, c("yes","no"))
	disperse <- match.arg(disperse, c("no","yes"))
	interpolate <- match.arg(interpolate, c("no","yes"))
	usefailed <- match.arg(usefailed, c("no","yes"))
	trim <- match.arg(trim, c("no","yes"))
	matchmax <- match.arg(matchmax, c("no","yes"))
	deorphan <- match.arg(deorphan, c("no","yes"))

	sqs.input.path <- "./sqs/input/"; dir.create(file.path(sqs.input.path), showWarnings = FALSE)
	sqs.output.path <- "./sqs/output/"; dir.create(file.path(sqs.output.path), showWarnings = FALSE)
	do.call(file.remove,list(list.files(sqs.input.path, full.names = TRUE)))
	do.call(file.remove,list(list.files(sqs.output.path, full.names = TRUE)))

	#names of variables needed by SQS Perl script
	sqs_variables <- c("collection_no","occurrence.genus_name","occurrence.species_reso","occurrence.species_name","occurrence.abund_value","stage","subepoch","epoch","FR2_bin","Peters.interval","X10_my_bin","ma_max","ma_min","max_ma","min_ma","collection.reference_no","species_extant","genus_extant") #"interpolated_base","interpolated_top","interval_base","interval_top",
	data <- data[,colnames(data) %in% sqs_variables] #dropping all variables except those needed ensures that the data isn't corrupted by special characters like ", / or \

	if (is.null(collections)) {
		write.csv(data, file = paste(sqs.input.path, "sqs_input_data.csv", sep = ""), quote = F)
	}
	if (!is.null(collections)) {
		write.csv(data[data$collection_no %in% collections, ], file = paste(sqs.input.path, "sqs_input_data.csv", sep = ""), quote = F)
	}

	sqs.job <- list.files(sqs.input.path) #get list of all .csv input files from input directory
	sqs.job <- sub("[.][^.]*$", "", sqs.job, perl = TRUE) #remove .csv extensions

	# cat(paste("RUNNING SQS FOR ", sqs.job, " at quorum level ", quorum.level, "\n", sep = ""))
	system(noquote(paste('perl ./functions/quorumSubsample4-2.pl --PATH "', sqs.input.path, '" ',
			     '--FILES "', sqs.job, '.csv" ',
			     '--TRIALS ', ntrials, ' ',
			     '--SCALE "', scale, '" ',
			     '--EXACT "', exact, '" ',
			     '--BYCOLLECTION "', bycollection, '" ',
			     '--SINGLETONS "', singletons, '" ',
			     '--COLLSPERREF ', collsperref, ' ',
			     '--REFQUOTA ', refquota, ' ',
			     '--DEORPHAN "',deorphan, '" ',
			     '--SEQUENTIAL "', sequential, '" ',
			     '--DISPERSE "', disperse, '" ',
			     '--METHOD "', method, '" ',
			     '--QUORUM ', quorum.level, ' ',
			     '--RANK "', rank, '" ',
			     '--INTERPOLATE "', interpolate, '" ',
			     '--USEFAILED "', usefailed, '" ',
			     '--TRIM "', trim, '" ',
			     '--MINCOLLBYBIN ', mincollbybin, ' ',
			     '--MATCHMAX "', matchmax, '"',
			     sep = "")), intern = FALSE, ignore.stdout = ignore.stdout, ignore.stderr = ignore.stderr)
	#read and format .csv file with results from perl script
	sqs.out <- read.csv("./sqs/output/quorumSubsample4.txt", skip = 0, header = T, sep = "	", stringsAsFactors = FALSE)

	return(sqs.out)

}




