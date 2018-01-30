# performs shareholder quorum subsampling on a count vector
# written by John Alroy on 29 April 2017

sqs<-function(n,quorum)	{
	S <- length(n)
	N <- sum(n)
	if (1 - length(which(n == 1)) / N < quorum)
		return(NA)
	lastS <- 1
	for (i in 1:N)	{
		not <- 0
		thisSingletons <- 0
		all <- lchoose(N , i)
		for (j in 1:length(n))	{
			not <- not + exp(lchoose(N - n[j] , i) - all)
			thisSingletons <- thisSingletons + exp(log(n[j]) + lchoose(N - n[j] , i - 1) - all)
		}
		thisS <- S - not
		if (1 - thisSingletons / i > quorum)
			break
		lastS <- thisS
	}
	return(as.numeric((lastS + thisS) / 2))
}
