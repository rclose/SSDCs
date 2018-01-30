# written by John Alroy on 29 April 2017

sqsbyref<-function(occurrences,samples,references,quorum=0.5,trials=100,sample.quota=1)	{
	# convert taxon names to numbers to speed things up
	occurrences <- xtfrm(occurrences)
	richness <- rep(NA,trials)
	# create a lookup of reference numbers corresponding to sample numbers
	lookup <- array()
	for (i in 1:length(samples))
		lookup[samples[i]] <- references[i]
	# create a list of unique references
	scrambled <- unique(references)
	for (i in 1:trials)	{
		# randomly order the list
		scrambled <- sample(scrambled)
		subsamples <- array()
		cross <- array()
		crosses <- 0
		lastu <- 0
		lastr <- 0
		for (j in 1:length(scrambled))	{
			# extract sample indices associated with this reference
			s <- which(lookup == scrambled[j])
			# impose the sample quota
			if (length(s) < sample.quota)
				next
			# draw sample.quota samples and append to a growing list
			subsamples[((j - 1)*sample.quota+1):(j*sample.quota)] <- s[sample(length(s),sample.quota)]
			# compute the abundance distribution
			# tabulate is faster than table
			n <- tabulate(occurrences[which(samples %in% subsamples)])
			# tabulate returns zero counts, unlike table
			n <- n[which(n > 0)]
			# compute Good's u
			u <- 1 - length(which(n == 1)) / sum(n)
			# record richness at each step that surpasses the quorum
			if (lastu < quorum && u > quorum)	{
				crosses <- crosses + 1
				cross[crosses] <- (length(n) + lastr) / 2
			} else if (lastu < quorum && u == quorum)	{
				crosses <- crosses + 1
				cross[crosses] <- length(n)
			}
			lastu <- u
			lastr <- length(n)
		}
		# this trial's richness is the median of all values at crossing points
		richness[i] <- median(cross)
		# NA values caused by failure to cross are unacceptable
		if (is.na(richness[i]))
			return(c(NA,NA,NA))
	}
	# return the median and 95% confidence interval
	return(quantile(richness,probs=c(0.025,0.5,0.975)))
}
