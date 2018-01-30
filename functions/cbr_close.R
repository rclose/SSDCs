rescaledMultitonRatio <- function(n) {
	S <- length(n)
	singletons <- length(which(n == 1))
	N <- sum(n)
	N / (singletons + 1) * (S - singletons) / S
}

multitonRatio <- function(n) {

	S <- length(n)
	s1 <- sum(n == 1)
	(S - s1) / S

}


cbr <- function(occurrences, quora = seq(from = 0, to = 1, by = 0.1), trials = 100) {

	#modified versions of TRiPS code that doesn't return CIs to speed things up

	doTRiPS_abs_noCI <- function(abs,t=1){
		# Performs TRiPS on the count of observations in abs. t=1 (default) is assumed duration
		p_lam = estimatePoiss_noCI(rep(t,length(abs)),abs);
		p_bino = 1-exp(-p_lam*t);
		estimatetrue_noCI(length(abs),p_bino[1])[1]

	}

	estimatePoiss_noCI <- function(dTs,Occs){
		# Negative log likelihood of the observations.
		nllnow <- function(lambda){ loglikepoissint(lambda,dTs,Occs)}
		# Getting maximum likelihood estimate of the poisson rate.
		fit1 <- mle(nllnow,start = list(lambda=1),nobs=length(Occs),method="Brent",lower = 1e-8,upper=30)
		return(coef(fit1))
	}

	estimatetrue_noCI <- function(nobs,binomprob) {
		if (!is.na(binomprob)){
			n <- seq(0,nobs/binomprob * 4+10)
			liks <- log(dbinom(nobs,size=n,prob=binomprob))

			return(n[which.max(liks)])
		} else {
			return(NA)
		}
	}

	TRiPS_noCI <- function(n) {
		if (!any(n > 1)) {
			return(NA)
		} else {
			tryCatch(doTRiPS_abs_noCI(n, t = 1), error = function(err) NA)
		}
	}

	# convert taxon names to numbers to speed things up
	occurrences <- xtfrm(occurrences)

	doSQS <- function(i) {
		scrambled <- sample(occurrences)
		cross_chao <- cross_lambda <- cross_raw <- cross_trips <- rep(list(array()), length(quora))
		crosses <- rep(0, length(quora)) #vector with one element for each quorum level
		lastu <- 0
		lastr_chao <- lastr_lambda <- lastr_raw <- lastr_trips <- 0
		for (j in 1:length(scrambled))	{
			# compute the abundance distribution
			# tabulate is faster than table
			n <- tabulate(scrambled[1:j])
			# tabulate returns zero counts, unlike table
			n <- n[which(n > 0)]
			# compute Good's u
			u <- 1 - (sum(n == 1) / sum(n))
			r_raw <- length(n)
			# r_trips <- TRiPS(n)
			r_chao <- Chao1(n)
			r_lambda <- lambda5(n)
			# record richness at each step that surpasses the quorum
			for (q in 1:length(quora)) {
				if (lastu < quora[q] && u > quora[q]) {
					crosses[q] <- crosses[q] + 1
					cross_raw[[q]][crosses[q]] <- (r_raw + lastr_raw) / 2
					# cross_trips[[q]][crosses[q]] <- (r_trips + lastr_trips) / 2
					cross_chao[[q]][crosses[q]] <- (r_chao + lastr_chao) / 2
					cross_lambda[[q]][crosses[q]] <- (r_lambda + lastr_lambda) / 2
				} else if (lastu < quora[q] && u == quora[q]) {
					crosses[q] <- crosses[q] + 1
					cross_raw[[q]][crosses[q]] <- r_raw
					# cross_trips[[q]][crosses[q]] <- r_trips
					cross_chao[[q]][crosses[q]] <- r_chao
					cross_lambda[[q]][crosses[q]] <- r_lambda
				}
			}
			lastu <- u
			lastr_raw <- r_raw
			# lastr_trips <- r_trips
			lastr_chao <- r_chao
			lastr_lambda <- r_lambda
		}
		# this trial's richness is the median of all values at crossing points
		output <- cbind(
			sapply(cross_raw, median, na.rm = T),
			# sapply(cross_trips, median, na.rm = T),
			sapply(cross_chao, median, na.rm = T),
			sapply(cross_lambda, median, na.rm = T)
		)
		output
	}
	test <- doSQS(i = 1)
	if (all(is.na(test))) {
		# return(c(NA,NA,NA))
		matrix(data = NA, nrow = length(quora), ncol = ncol(test))
	} else {
		richness <- mclapply(1:trials, function(i) doSQS(i), mc.cores = detectCores())
		richness <- abind(richness, along = 3)
		richness1 <- matrix(data = NA, nrow = length(quora), ncol = dim(richness)[2])
		for (i in 1:dim(richness)[2]) {
			richness1[,i] <- apply(richness[,i,], MARGIN = 1, function(x) {
				if (any(is.na(x))) {
					# NA values caused by failure to cross are unacceptable
					# c(NA,NA,NA)
					NA
				} else 	{
					# return the median and 95% confidence interval
					# quantile(richness,probs=c(0.025,0.5,0.975))
					median(x)
				}
			}
			)
		}
		colnames(richness1) <- c(
			"raw",
			# "trips",
			"chao",
			"lambda"
		)
		richness1
	}
}

#rarefy by multiton ratio
mbr <- function(occurrences, mr = seq(from = 0, to = 1, by = 0.1), trials = 100) {
	# TRiPS <- function(n) {
	# 	if (!any(n > 1)) {
	# 		return(NA)
	# 	} else {
	# 		tryCatch(doTRiPS_abs(n, t = 1), error = function(err) matrix(nrow = 3, ncol = 3, dimnames = list(c("Sampling rate","Sampling probability","Estimated richness"),c("MLE","lower CI","upper CI"))))[3,1]
	# 	}
	# }
	# convert taxon names to numbers to speed things up
	occurrences <- xtfrm(occurrences)

	doSQS <- function(i) {
		scrambled <- sample(occurrences)
		cross_chao <- cross_lambda <- cross_raw <- cross_trips <- rep(list(array()), length(mr))
		crosses <- rep(0, length(mr)) #vector with one element for each quorum level
		lastu <- 0
		lastr_chao <- lastr_lambda <- lastr_raw <- lastr_trips <- 0
		for (j in 1:length(scrambled))	{
			# compute the abundance distribution
			# tabulate is faster than table
			n <- tabulate(scrambled[1:j])
			# tabulate returns zero counts, unlike table
			n <- n[which(n > 0)]
			# compute rescaled multiton ratio
			u <- multitonRatio(n)
			r_raw <- length(n)
			# r_trips <- TRiPS(n)
			r_chao <- Chao1(n)
			r_lambda <- lambda5(n)
			# record richness at each step that surpasses the quorum
			for (q in 1:length(mr)) {
				if (lastu < mr[q] && u > mr[q]) {
					crosses[q] <- crosses[q] + 1
					cross_raw[[q]][crosses[q]] <- (r_raw + lastr_raw) / 2
					# cross_trips[[q]][crosses[q]] <- (r_trips + lastr_trips) / 2
					cross_chao[[q]][crosses[q]] <- (r_chao + lastr_chao) / 2
					cross_lambda[[q]][crosses[q]] <- (r_lambda + lastr_lambda) / 2
				} else if (lastu < mr[q] && u == mr[q]) {
					crosses[q] <- crosses[q] + 1
					cross_raw[[q]][crosses[q]] <- r_raw
					# cross_trips[[q]][crosses[q]] <- r_trips
					cross_chao[[q]][crosses[q]] <- r_chao
					cross_lambda[[q]][crosses[q]] <- r_lambda
				}
			}
			lastu <- u
			lastr_raw <- r_raw
			# lastr_trips <- r_trips
			lastr_chao <- r_chao
			lastr_lambda <- r_lambda
		}
		# this trial's richness is the median of all values at crossing points
		output <- cbind(
			sapply(cross_raw, median, na.rm = T),
			# sapply(cross_trips, median, na.rm = T),
			sapply(cross_chao, median, na.rm = T),
			sapply(cross_lambda, median, na.rm = T)
		)
		output
	}
	test <- doSQS(i = 1)
	if (all(is.na(test))) {
		# return(c(NA,NA,NA))
		matrix(data = NA, nrow = length(mr), ncol = ncol(test))
	} else {
		richness <- mclapply(1:trials, function(i) doSQS(i), mc.cores = detectCores())
		richness <- abind(richness, along = 3)
		richness1 <- matrix(data = NA, nrow = length(mr), ncol = dim(richness)[2])
		for (i in 1:dim(richness)[2]) {
			richness1[,i] <- apply(richness[,i,], MARGIN = 1, function(x) {
				if (any(is.na(x))) {
					# NA values caused by failure to cross are unacceptable
					# c(NA,NA,NA)
					NA
				} else 	{
					# return the median and 95% confidence interval
					# quantile(richness,probs=c(0.025,0.5,0.975))
					median(x)
				}
			}
			)
		}
		colnames(richness1) <- c(
			"raw",
			# "trips"
			"chao",
			"lambda"
		)
		richness1
	}
}
