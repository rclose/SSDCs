#### clear workspace ####
rm(list = ls())

#### load packages ####
source("./subscripts/01. load-packages.R", echo = F)

#### analysis settings ####
#put results in a new date-stamped output folder?
new.output.folder <- TRUE

#clear out old plots if using old folder?
clear_plots <- FALSE

#choose your clade (tetrapoda, dinosauromorpha)
clade <- "dinosauromorpha"
# load("./saved_settings.RData") #uncomment when using master script to loop through clades

#use stages for intervals, or custom timebins?
use.stages <- T

#set first and last intervals for each clade
if (clade == "tetrapoda") {
	if (use.stages == FALSE) {
		# first.int <- "D4"; last.int <- "Ng4"
		first.int <- "Tr1"; last.int <- "Pg2"
	} else {
		# first.int <- "Frasnian"; last.int <- "Holocene"
		first.int <- "Induan"; last.int <- "Maastrichtian"
	}
} else if (clade == "dinosauromorpha") {
	if (use.stages == FALSE) {
		first.int <- "Tr1"; last.int <- "K8"
	} else {
		first.int <- "Induan"; last.int <- "Maastrichtian"
	}
}

#choose beginning and end of DTTTT analysis
beginning <- 1866; end <- 2016

#choose period of DTTTT analysis (i.e., every n years)
period <- 10

#include global occurrences as a 'region'?
include.global <- T

#choose rank at which to run analyses
rank <- "species"

#turn calculation of richness estimators on and off
run_PGS_calcs <- T
run_SQS_Perl_1 <- T #SQS with quota of 3 collections per reference
run_SQS_Perl_2 <- T #SQS with no collection-per-reference quota
run_iNEXT <- T
run_lambda_5 <- T
run_chao <- T
run_CR <- T
run_TRiPS <- T

#choose quorum (SQS) and quota (CR) levels at which to run analyses
quorum.levels <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8); names(quorum.levels) <- paste("Quorum", quorum.levels, sep = " ")
quota.levels <- c(100,200,300,400,800,1600); names(quota.levels) <- paste("Quota", quota.levels, sep = " ")

#choose quorum level for summary plots, etc.
quorum <- "Quorum 0.4"

#set number of subsampling trials for real data
ntrials <- 1000


#### null distribution settings ####
#include null distribution of randomly-ordered discoveries?
nulldist <- T

#set number of subsampling trials for null distribution
nulltrials <- 1000

#set number of shuffled replicates in null
nshuff <- 100

#drop pseudo-references in PaleoDB that report faulty ref_pubyr data, such as UCMP database and Sepkoski compendia?
dropdodgyrefs <- T


#### load and clean occurrence data ####
source("./subscripts/02. load-occs.R", echo = T)


#### define intervals ####
source("./subscripts/03. make_intervals.R", echo = T)


#### define regions ####
source("./subscripts/04. define_regions.R", echo = T)


#### wipe old plots and/or create directory structure ####
source("./subscripts/05. create-directories.R", echo = T)


#### plotting settings ####
#set layout of multipanel plots depending on number of intervals used
multipanel.arr <- c(ceiling(nrow(intervals)/4),4)
multipanel.width <- multipanel.arr[2] * 4; multipanel.height <- multipanel.arr[1] * 4
if (length(quorum.levels) > 8) {quorum.mfrows <- c(6,3)} else if (length(quorum.levels) <= 8) {quorum.mfrows <- c(4,4)}
multiquorum.width <- quorum.mfrows[1]*4; multiquorum.height <- quorum.mfrows[2]*4
uniplot.width <- 8; uniplot.height <- 8
legend.cex <- 0.55
black.lwd <- 4; colour.lwd <- 3
logged <- log_xy <- ""
log_x <- ""
dtttt.xlim <- rev(range(intervals$midpoint))
#set transparency of null distribution
null.alpha <- 0.1
#set line width of null distribution (vertical lines showing min/max)
null.lwd <- 5
#define boundaries for plotting lines
period.boundaries <- c(541.0, 485.4, 443.8, 419.2, 358.9, 298.9, 252.17, 201.3, 145, 66)


#### trim occurrence data to desired interval ####
pbdb_data <- filter(pbdb_data, ma_max <= intervals[intervals$bin == first.int, "FAD"] & ma_min >= intervals[intervals$bin == last.int, "LAD"])


#### subset data by region ####
#create a unique rowID to help with filtering obsolete occurrence IDs
pbdb_data$rowID <- paste("rowID", 1:nrow(pbdb_data), sep = "")

#create a backup of the regions object
regions.bak <- regions

#make a 'global' region with all countries
global <- unique(pbdb_data$country)

#add global to end of regions if specified
if (include.global == TRUE) {
	regions[["global"]] <- global
}

#subset to just some regions
regions <- regions[c("NAm","global")]

if (clade == "dinosauromorpha") {
	regions <- regions[c("global")]
} else if (clade == "tetrapoda") {
	regions <- regions[c("NAm")]
}

#subset data by region and combine in list
occs <- lapply(regions, function(x) pbdb_data[!is.na(match(pbdb_data$country,x)),])


#### define timescale colours ####
#get names of first stage in each interval
stages <- sapply(1:nrow(intervals), function(x) strsplit(intervals$stages, split = ", ", fixed = T)[[x]][1]); names(stages) <- intervals$bin
#subset timescale data from package geoscale, keeping ages only
timescale <- read.csv("./input-data/timescale_cols.csv", stringsAsFactors = F)
timescale <- filter(timescale, Type == "Age"); timescale <- timescale[!is.na(timescale$Name),]
#get rows that correspond to stages that were drawn from interval definitions
timescale <- filter(timescale, Name %in% gsub("Stage ","Stage",stages))
#use RGB values to make colours that correspond to each bin
interval.cols <- timescalecols <- rgb(red = timescale$Col_R, green = timescale$Col_G, blue = timescale$Col_B, maxColorValue = 255, names = intervals$bin)


#### make all the empty lists we need ####
if (run_SQS_Perl_1 == TRUE) {null.sqs.1 <- real.sqs.1 <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
if (run_SQS_Perl_2 == TRUE) {null.sqs.2 <- real.sqs.2 <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
if (run_iNEXT == TRUE) {null.inext <- real.inext <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
if (run_lambda_5 == TRUE) {null.lambda5 <- real.lambda5 <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
if (run_chao == TRUE) {null.chao <- real.chao <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
if (run_CR == TRUE) {null.CR <- real.CR <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
if (run_TRiPS == TRUE) {null.trips <- real.trips <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
if (run_PGS_calcs == TRUE) {pgs <- setNames(vector(mode = "list", length = length(regions)), names(regions))}
histoccRowIDs <- goodsu <- multiton_count <- singleton_count <- newoccs <- newtax <- numoccs <- setNames(vector(mode = "list", length = length(regions)), names(regions))


#### make vector of years ####
times <- rev(seq(from = beginning, to = end, by = period))
names(times) <- paste("year", times, sep = "_")
decades <- seq(from = min(times), to = max(times), by = 10)
names(decades) <- paste("year", decades, sep = "_")
tentimes <- times[seq(from = 1, to = length(times), length.out = 10)]


##### make colour palettes #####
region.cols <- setNames(adjustcolor(brewer.pal(length(regions), "Set3"), alpha.f = 0.9), names(regions)) #colour palette for plotting subMST trees onto global MST @UNDO
region.cols <- region.cols[1:length(regions)]
timeseries.cols <- rainbow(n = length(times))
quorum.cols <- rainbow(n = length(quorum.levels))
black2 <- adjustcolor("black", alpha.f = 0.5)
interval.cols <- adjustcolor(timescalecols, alpha.f = 1); names(interval.cols) <- intervals$bin
nullint.cols <- adjustcolor(darken(interval.cols), alpha.f = 0.4); names(nullint.cols) <- intervals$bin


##### do the analyses #####
for (rgn in 1:length(regions)) {
	if (nrow(occs[[rgn]]) > 2) {
		source("./subscripts/06. runrunrun.R", echo = T)
	}
}

if (run_SQS_Perl_1 == T | run_SQS_Perl_2 == T) {
	real.sqs <- real.sqs.2 #choose your SQS flavour for plotting
}

##### flatten lists and arrays into tidy data frames #####
source("./subscripts/06.5 flatten-lists.R", echo = T)


#### save the workspace ####
save.image(paste("./saved_workspaces/saved_workspace_", clade, "_stages_equals_", use.stages,".RData", sep = ""))
# load(paste("./saved_workspaces/saved_workspace_", clade, "_stages_equals_", use.stages, ".RData", sep = ""))


##### plot the manuscript figures #####
source("./subscripts/10. figure-code.R", echo = T)


##### run richness estimator simulations ####
source("./subscripts/11.0 Richness Estimator Simulations.R", echo = T)
