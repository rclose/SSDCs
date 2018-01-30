collRegions <- function(dat = NULL, colls = NULL, regions = NULL) {

	countries <- unique(dat[dat$collection_no %in% colls,]$country) #get countries present in current fixed-spread MST
	region <- sapply(regions, function(x) any(countries %in% x))
	region <- names(region[region == TRUE])
	if (length(region) == 0) {region <- "Undefined"}
	if (length(region) > 1) {region <- paste(region, collapse = ", ")}
	return(region)

}
