lambda5<-function(n)    {
        S <- length(n)
        N <- sum(n)
        s1 <- length(which(n == 1))
        if (s1 == 0)
                return(S)
        if (s1 == S)
                return(S + s1^2 / 2)
	l <- 1
	lastl <- l * 2
	target <- log(N/(S - s1)) * s1/S
	estimate <- target * 2
	i <- 1
	while (i < 100 && abs(log(l / lastl)) > 0.0001)	{
		i <- i + 1
		lastl <- l
		estimate <- log(l / (1 - exp(-l) - l * exp(-l))) * l * exp(-l)/(1 - exp(-l))
		l <- l * (estimate / target)^(1/3)
	}
	return(S / (1 - exp(-l)))
}

