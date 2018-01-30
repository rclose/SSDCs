# coll_data <- tbl_df(read.csv("https://paleobiodb.org/data1.2/colls/list.csv?base_name=Metazoa&show=full,loc,paleoloc,stratext,lithext,geo,methods,ref,refattr,secref,crmod", header = T, stringsAsFactors = FALSE)); save(coll_data, file = "./input-data/coll_data.RData")
load("./input-data/coll_data.RData")
if (clade == "tetrapoda") {
	# tetrapoda <- tbl_df(read.csv("https://paleobiodb.org/data1.2/occs/list.csv?base_name=Tetrapoda^Aves^Chiroptera^Pterosauria&ident=all&pres=regular&show=attr,class,classext,genus,subgenus,ident,abund,ecospace,taphonomy,ctaph,etbasis,pres,coll,coords,loc,paleoloc,prot,stratext,lithext,env,geo,methods,resgroup,ref,refattr,ent,entname,crmod", header = T, stringsAsFactors = FALSE)); save(tetrapoda, file = "./input-data/tetrapoda.RData")
	load("./input-data/tetrapoda.RData")
} else if (clade == "dinosauromorpha") {
	# dinosauromorpha <- tbl_df(read.csv("https://paleobiodb.org/data1.2/occs/list.csv?base_name=Dinosauromorpha^Aves&ident=all&pres=regular&show=attr,class,classext,genus,subgenus,ident,abund,ecospace,taphonomy,ctaph,etbasis,pres,coll,coords,loc,paleoloc,prot,stratext,lithext,env,geo,methods,resgroup,ref,refattr,ent,entname,crmod", header = T, stringsAsFactors = FALSE)); save(dinosauromorpha, file = "./input-data/dinosauromorpha.RData")
	load("./input-data/dinosauromorpha.RData")
}

#copy desired occurrence data to new object
pbdb_data <- get(clade); rm(list = clade)

if (rank == "species") {
pbdb_data <- filter(pbdb_data, (identified_rank %in% c("species")))
} else if (rank == "genus") {
	pbdb_data <- filter(pbdb_data, (identified_rank %in% c("species","genus")))
}

#look at most prolific references
# paste(pbdb_data$ref_author, " ", pbdb_data$ref_pubyr, ", ", pbdb_data$reference_no, sep = "") %>% table %>% sort(decreasing = F) %>% View

#drop refs that don't really correspond to discovery years for occurrences
if (dropdodgyrefs == TRUE) {
	droprefs <- c(
		"13" #J. J. Sepkoski Jr. 1998. Rates of speciation in the fossil record.
		,"6294" #J. Alroy. 2002. Synonymies and reidentifications of North American fossil mammals.
		,"13103" #UCMP Database. 2005. UCMP collections database.
	)
	pbdb_data <- filter(pbdb_data, !(reference_no %in% droprefs))
}


#merge in variables from collection data
pbdb_data$collection.reference_no <- coll_data[match(pbdb_data$collection_no, coll_data$collection_no), ]$reference_no #slot in reference numbers for collections
pbdb_data$collection.ref_pubyr <- coll_data[match(pbdb_data$collection_no, coll_data$collection_no), ]$ref_pubyr #slot in ref_pubyrs for collections
rm(list = "coll_data") #remove from env to save memory

#rename variables to match those used by Alroy's SQS Perl script v. 4.3
pbdb_data$occurrence.species_name <- pbdb_data$species_name; pbdb_data$species_name <- NULL
pbdb_data$occurrence.genus_name <- pbdb_data$primary_name; pbdb_data$genus <- NULL #@TODO verify that $primary_name is better than $genus
pbdb_data$occurrence.reference_no <- pbdb_data$reference_no; pbdb_data$reference_no <- NULL
pbdb_data$occurrence.species_reso <- pbdb_data$species_reso; pbdb_data$species_reso <- NULL
pbdb_data$ma_max <- pbdb_data$max_ma; pbdb_data$max_ma <- NULL
pbdb_data$ma_min <- pbdb_data$min_ma; pbdb_data$min_ma <- NULL
pbdb_data$max_interval <- pbdb_data$early_interval; pbdb_data$early_interval <- NULL
pbdb_data$min_interval <- pbdb_data$late_interval; pbdb_data$late_interval <- NULL
pbdb_data$occurrence.abund_value <- pbdb_data$abund_value; pbdb_data$abund_value <- NULL
pbdb_data$country <- pbdb_data$cc
pbdb_data$occurrence.binomial <- paste(pbdb_data$occurrence.genus_name, pbdb_data$occurrence.species_name, sep = " ")
pbdb_data <- subset(pbdb_data, pbdb_data$occurrence.genus_name != "")
pbdb_data$paleolatdec <- pbdb_data$paleolat; pbdb_data$paleolngdec <- pbdb_data$paleolng

#Then use lists of trace, egg and marine names stored in text files to remove those occurrences
trace.terms <- scan("./input-data/exclude-terms/trace-terms.txt", what = "character"); trace.terms <- trace.terms[trace.terms != ""]
egg.terms <- scan("./input-data/exclude-terms/egg-terms.txt", what = "character"); egg.terms <- egg.terms[egg.terms != ""]
marine.terms <- scan("./input-data/exclude-terms/marine-terms.txt", what = "character"); marine.terms <- marine.terms[marine.terms != ""]

#By this point we have a list of taxon names that we'd like to include
exclude.terms <- c(
	trace.terms,
	egg.terms,
	# "Jeholornis", #dinobird
	"Agelaius", #random "dinosaur" in the Cenozoic
	"Crocodylus", "Crocodilus", "Alligator", "Testudo", "Lacerta", #wastebasket genera
	"Fenestrosaurus","Ovoraptor","Ornithoides" #names from popular article by Osborne (1925) that should not be in the database
)
exclude.terms <- exclude.terms[exclude.terms != ""] #remove any blank entries that may have crept in

#strip out everything that needs to be excluded
pbdb_data <- pbdb_data[!(pbdb_data$order %in% exclude.terms | pbdb_data$family %in% exclude.terms | pbdb_data$occurrence.genus_name %in% exclude.terms), ]

#remove occurrences with soft-tissue preservation or traces
# pbdb_data <- pbdb_data[!grepl("soft|trace|coprolite",pbdb_data$pres_mode), ]
pbdb_data <- pbdb_data[!grepl("trace|coprolite",pbdb_data$pres_mode), ]

#remove informal identifications (e.g. "informal indet. 1", etc.)
pbdb_data <- filter(pbdb_data, occurrence.species_reso != "informal")

#remove marine tetrapods
pbdb_data <- pbdb_data[!(pbdb_data$order %in% marine.terms | pbdb_data$family %in% marine.terms | pbdb_data$occurrence.genus_name %in% marine.terms), ]


