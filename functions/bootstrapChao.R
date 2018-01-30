bootstrapChao <- function(n)	{
	l <- array()
	for (i in 1:10000)	{
		l[i] <- Chao1(sample(n,replace=T))
	}
	quantile(l,probs=c(0.025,0.5,0.975))
}

