goodsU <- function(occvec) {
	n <- sum(occvec)
	f1 <- sum(occvec == 1); f2 <- sum(occvec == 2)
	# out <- 1 - (f1 / sum(occvec)) #singleton-only coverage estimator
	out <- 1 - (f1/n) * (((n - 1) * f1) / ((n - 1) * f1 + 2 * f2)) #Chao and Jost (2012) coverage estimator using singletons and doubletons
	out[is.nan(out) | is.infinite(out)] <- NA
	return(out)
}
