Chao1 <- function(ab)	{
	return(length(ab) + length(which(ab == 1)) * (length(which(ab == 1)) - 1) / (2 *  (length(which(ab == 2)) + 1)))
}

otherEstimators <- function(s1) {
	raw = length(s1)
	chao = Chao1(s1)
	lambda = lambda5(s1)
	# multi = multiton(s1, target = 10)
	if (exists("out")) {rm(list = "out")}
	trips <- tryCatch(doTRiPS_abs(s1, t = 1), error = function(err) matrix(nrow = 3, ncol = 3, dimnames = list(c("Sampling rate","Sampling probability","Estimated richness"),c("MLE","lower CI","upper CI"))))
	trips <- trips[3,1]
	return(data.frame(raw = raw, chao = chao, lambda = lambda, trips = trips))
}

coverageEstimators <- function(s1, q) {
	sqs_inext = estimateD(x = s1, datatype = "abundance", base = "coverage", level = q, conf = NULL)[1,4]
	sqs_inext = as.numeric(sqs_inext)
	sqs_alroy = simpleSQS(n = s1, quorum = q)
	return(data.frame(sqs_inext = sqs_inext, sqs_alroy = sqs_alroy))
}

multitonRatio <- function(n) {

	S <- length(n)
	s1 <- sum(n == 1)
	(S - s1) / S

}

rescaledMultitonRatio <- function(n) {
	S <- length(n)
	singletons <- length(which(n == 1))
	N <- sum(n)
	N / (singletons + 1) * (S - singletons) / S
}

sampleCoverage <- goodsU <- function(n) {
	N <- sum(n)
	f1 <- sum(n == 1); f2 <- sum(n == 2)
	# out <- 1 - (f1 / sum(n)) #singleton-only coverage estimator
	if (f1 == 0 && f2 == 0) {
		return(1)
	} else {
		out <- 1 - (f1/N) * (((N - 1) * f1) / ((N - 1) * f1 + 2 * f2)) #Chao and Jost (2012) coverage estimator using singletons and doubletons
	}
	out[is.nan(out) | is.infinite(out)] <- NA
	return(out)
}
