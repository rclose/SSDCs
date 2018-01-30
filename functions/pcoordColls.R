#find all collections with identical palaeocoordinates to those listed
pcoordColls <- function(dat, coll_nos) {

	dat <- as.data.frame(dat)
	dat$pcoords <- paste(dat$paleolatdec,dat$paleolngdec)
	coll_pcoords <- unique(dat[dat$collection_no %in% as.integer(coll_nos), ]$pcoords)
	coincident_colls <- unique(subset(dat, pcoords %in% coll_pcoords)$collection_no)
	return(coincident_colls)

}

