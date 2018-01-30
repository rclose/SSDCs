## Functions:
doTRiPS_abs <- function(abs,t=1){
  # Performs TRiPS on the count of observations in abs. t=1 (default) is assumed duration
  p_lam = estimatePoiss(rep(t,length(abs)),abs);
  p_bino = 1-exp(-p_lam*t);
  N_true = c(estimatetrue(length(abs),p_bino[1])[1],
             min(estimatetrue(length(abs),p_bino[3])),
             max(estimatetrue(length(abs),p_bino[2])))
  out = array(NA,c(3,3))

  out[1,]=p_lam
  out[2,]=p_bino
  out[3,]=N_true
  rownames(out)<-c("Sampling rate","Sampling probability","Estimated richness")
  colnames(out)<-c("MLE","lower CI","upper CI")
  return(out)
}


get_ltt <- function(Out){
  # Function to tally lineages through time from the output of a simulateBDF.
  # Out is an array of size (no_lineages) by (3)
  # where [,1] is start of lineage
  # [,2] end of lineages
  # [,3] number of fossils left (not used for this function)
  spectimes=sort(unique(as.vector(Out[,1:2])))
  N_tt = array(NA,c(length(spectimes)));
  # Need to think here; spectimes is when events happen. So we want N_tt[ii] to be
  # the number of lineages between spectimes[ii] and spectimes[ii+1]. I.e. if a lineage
  # started at spectimes[ii] we includ it. If it dies out at spectimes[ii] we do not.
  for (ii in 1:(length(spectimes)-1)){
    N_tt[ii] = sum(Out[,2]>spectimes[ii] & Out[,1]<=spectimes[ii])
  }
  N_tt[ii+1]=N_tt[ii];

  return(list(spectimes,N_tt))
}

doTRiPS <- function(Out){
  # doing the whole TRiPS estimation with the Out array from a simulation from simulateBDF
  #
  # Occs = Out[Out[,3]>0,3];
  # dTs = max(Out[,2]); # assuming this is the duration.
  p_lam = estimatePoiss(rep(max(Out[,2]),sum(Out[,3]>0)),Out[Out[,3]>0,3])
  p_bino = 1-exp(-max(Out[,2])*p_lam);
  N_true = c(estimatetrue(sum(Out[,3]>0),p_bino[1])[1],
             min(estimatetrue(sum(Out[,3]>0),p_bino[3])),
             max(estimatetrue(sum(Out[,3]>0),p_bino[2])))
}



# Drawing a number/index from an empirical probability distribution.
Emprand <- function(x) {
  j <- runif(1) ## drawing random uniform number
  which(cumsum(x)/sum(x)>j)[1] #returning first entry in the cumsum/sum [empirical density function]
}

## Likelihood function of the poisson intensity conditioned on more than one occurrence.
likepoissint <- function(lambda,dt,nobs) {
  (((lambda*dt)^(nobs))/(factorial(nobs))*(exp(-lambda*dt)))/(1- exp(-lambda*dt))
  # These factorials blow up rather quickly and become computationally demanding.
  # Below the code generates a look-up table for relevant factorials, and
  # a function using the look-up table and calculating log likelihoods is generated.
  #
}

# Can we utilize the much simpler log likelihood for Poisson, even when we want to condition on >0?
# These are the same, but perhaps with some numerical stability
# JOS;: 030915
# Calculating and storing the log(factorial(n!)) needed for Poisson estimation,
tab_lfactorials = array(NA,c(1000,1));
tab_lfactorials[1] = log(1);
for (ii in 1:1000){
  tab_lfactorials[ii+1] = tab_lfactorials[ii] + log(ii+1)
}

likepoissint2 <-function(lambda,dt,nobs){
  nobs*log(lambda*dt) - length(nobs)*(lambda*dt) - log(1-exp(-lambda*dt)) - tab_lfactorials[nobs];
}


## Negative Log likelihood for many observations
loglikepoissint <- function(lambda,dtin,nobs) {
  if (length(nobs)!=length(dtin)){
    dt = rep(dt,length(nobs))
  }else{
    dt = dtin
  } # if not equal lengths, assume dt has 1 value and repeat

  ll <- 0
  for (ii in (1:length(dt))){
    # ll <- ll - log(likepoissint(lambda,dt[ii],nobs[ii]))
    # testing for likepoissint2 (using log likelihod for stability in calcs.)
    ll <- ll - likepoissint2(lambda,dt[ii],nobs[ii]);
  }
  return(ll)
}

# Function to return phat and pci
estimatePoiss <- function(dTs,Occs){
  # Negative log likelihood of the observations.
  nllnow <- function(lambda){ loglikepoissint(lambda,dTs,Occs)}
  # Getting maximum likelihood estimate of the poisson rate.
  fit1 <- mle(nllnow,start = list(lambda=1),nobs=length(Occs),method="Brent",lower = 1e-8,upper=30)
  # Defining the likelihood ratio function to estimate confidence intervals
  mycis2 <- function(l_x) ((2*(-nllnow(coef(fit1))+nllnow(l_x)))-qchisq(0.95,1))^2

  # lower CI
  ci1<-optim(0.9*coef(fit1),mycis2,method="Brent",lower=1e-6,upper=coef(fit1))
  # upper CI
  ci2<-optim(1.1*coef(fit1),mycis2,method="Brent",lower=coef(fit1),upper=40)
  return(c(coef(fit1),ci1$par,ci2$par))
}

# To estimate maximum likelihood of the true number of species given observed number and
# binomial sampling probability.
# pdf =
estimatetrue <- function(nobs,binomprob) {
  if (!is.na(binomprob)){
    n <- seq(0,nobs/binomprob * 4+10)
    liks <- log(dbinom(nobs,size=n,prob=binomprob))
    tmp <- n[which((2*(max(liks)-liks))<qchisq(0.95,1))]

    return(c(n[which.max(liks)],min(tmp),max(tmp)))
  } else {
    return(c(NA,NA,NA))
  }
}


# Code for simulating a BDF process. This code was developed for other uses as well,
# and is therefore slightly more complicated than needed for the PhilTrans paper.

SimulateBDF <- function(spec,mu,lambda,lambdavar,tmax,n_init,dt){
  ## Simulating a birth-death-fossilize model.
  # Written by Jostein Starrfelt (jostein.starrfelt[AT]ibv.uio.no)
  # Last checked 09.09.15
  # Should be updated to only use 'functional' approach to rates actually defined as function. Now
  # if any are defined as functions, all are treated as functions. This is inefficient. [JOS 051015]
  # Also there must be an error here;

  #
  # Updated to have potentially variable rates over time. These rates are inputted as functions taking the argument t-time.
  # This now works but then without lambdavar. Implementing lambdavar as a function that outputs a scaling factor for of sampling.
  # Lambdavar is then essentially a scaling distribution, which should have mean 1. The sampling rate for the whole clade is
  # lambda*lambdavar.
  # v3:
  # In this way one can implement whatever kind of distribution of differential sampling rates across lineages. The sampling rate for an
  # individual lineage at time t will then be lambda(t) * lambdavar(1). Future implementations might include also this as a function of time.
  # If not species do not vary in sampling rate, set lambdavar to 1. This is different from previous applications. Both the 'fixed' rate and
  # the functional code below needs to be augmented.

  # Implementing in fixed rate code; DONE.

  # Implementing in variable rate code; DONE, not tested.

  # Calling functions
  # When calling a function you can specify arguments by position, by complete name, or by partial name. Arguments are matched first by exact name (perfect matching), then by prefix matching, and finally by position.

  # How to implement
  #   spec = 0.01;
  #   mu = 0;
  #   lambda = 0.4;
  #   lambdavar = 0;
  #   tmax = 10;
  #   n_init = 100;
  # Splitting code into fixed and variable rates
  # Assuming fixed speciation and extinction rates for whole interval.
  if (any(list(typeof(spec),typeof(mu),typeof(lambda)) == "closure")){
    # Then one or more of these arguments are functions and needs to be treated accordingly.
    # If one or more are given as constants; make into functions;
    if (typeof(spec)=="double"){spold = spec; spec <- function(t){spold}; spec <- Vectorize(spec)}
    if (typeof(mu)=="double"){muold = mu; mu <- function(t){muold}; mu <- Vectorize(mu)}
    if (typeof(lambda)=="double"){lold = lambda; lambda <- function(t){lold}; lambda <- Vectorize(lambda)}
    # The following lines calculates the estimated mean number of cumulative species richness according to e BD process;
    # cf eq 41 in Kendall 1948
    myrhotmp <-  function(t){(mu(t)-spec(t))};
    myrho    <- function(t) {integrate(myrhotmp,0,t)}
    tmpfun <- function(t){exp(-myrho(t)$value)}
    mytmp  <- function(t){exp(-myrho(t)$value*spec(t))}
    nosp <- round(n_init*(1 + integrate(mytmp,0,tmax)$value))
    Out = array(NA,c(nosp*100,3));
    # Should put in dynamic checker for this array too, since it might be too small sometimes.
    Lamds = array(NA,c(nosp*100,1)); # for testing the lambdavalues
    done = 0; #checker for done calculating or not.
    Out[1:n_init,1] = 0; # initial lineages have start at t=0;
    tix = 1; # current lineage being simulated
    ntix = n_init+1; # current index into first entry in Out not populated with species yet.
    while (done==0){
      # while not done.
      # Here I will insert a switch to redefine the functions spec, mu and lambda, in cases where they will vary between individual lineages.
      # Exactly how this will be simulated remains to be seen.
      # if (doindfun==1){
      #redefine functions so they are lineage specific here. Most relevant for sampling function
      # }
      # Now simulating lineage tix.
      dead = 0; # switches to 1 when this lineage goes extinct;
      t_x = Out[tix+1-1,1]; # simulated time, this starts at time of speciation for lineage tix
      # dt = 1e-2; # Traveling along time with this resolution. The code below will integrate the relevant rates over many short intervals
      # of this duration.
      while (dead==0){
        if (runif(1)<integrate(mu,t_x,t_x+dt)$value){
          # if this is true it has died
          dead = t_x;
          Out[tix,2] = t_x;

        } else{
          t_x = t_x+dt; # just update time
        }
        if (t_x>tmax){
          # If not dead by end of simulation it has not died, but has duration until time tmax
          dead=-1; # To stop loop
          Out[tix,2] = tmax;  # This lineage has not gone extinct before tmax.
        }

      }
      # Now we have the full duration for this lineage.
      # Getting the speciation event for this lineage;
      # List of speciation times;
      S = array(NA,10*round(integrate(spec,0,tmax)$value));
      stix = 1; # Entry into the S array
      # dt = 1e-3;
      t_x = Out[tix+2-2,1]; # WHen this lineage first appears
      while (t_x<=Out[tix,2]){
        if (runif(1)<integrate(spec,t_x,t_x+dt)$value){
          # if this is true a speciation events happens now.
          # Need to check if S is big enough;
          if (stix>length(S)){
            S_tmp = S;
            S = array(NA,dim=c(round(length(S_tmp)*1.2),1)); # increasing size by 20 %
            S[1:stix-1] = S_tmp[1:stix-1];
          }
          S[stix] = t_x; # speciation time
          stix = stix+1; # increase counter
          t_x = t_x+dt;  # increase time
        } else {
          # nothing happens, time is updated
          t_x = t_x + dt;
        }
      }
      S = S[!is.na(S)]; # keeping only true speciation events from preallocated array.
      # show(S)
      # show(tix)

      if ((ntix+length(S)-1)>nrow(Out)){
        ## need to enlarge OUT
        oldout= Out;
        # Generating new array ~50 % bigger.
        Out = array(NA,c(round(nrow(Out)*1.5),3));
        Out[1:(ntix-1),]=oldout[1:(ntix-1),];
      }

      if (length(S)>0){
        Out[ntix:(ntix+length(S)-1)]=S;
        ntix =ntix +length(S); # updating index into first non-populated part of Out.
      }



      # Fossilization/Sampling;
      nofos = 0; # counting number of fossils for this lineage
      t_x = Out[tix+3-3,1]; #starting at first appearance of this lineage
      # Lambda scaling factor for this lineage
      lambdatmp = lambdavar(1);
      Lamds[tix] = lambdatmp;
      while (t_x<=Out[tix,2]){
        # show(t_x)
        # For the duration of this lineage
        if (runif(1)<lambdatmp*integrate(lambda,t_x,t_x+dt)$value){
          # If fossilization event inside the time t_x ... t_x + dt, adda fossil to count
          nofos = nofos+1;
          t_x = t_x + dt; # update time
        } else {
          # No fossilization event
          t_x = t_x + dt; # update time.
        }
      }
      Out[tix,3] = nofos;

      tix = tix+1; #now done simulating lineage tix, moving to tix+1
      # show(tix/ntix)
      if (tix==(ntix)){
        # if this was the last lineage in the array, we are done.
        done=1;
      }
    }


  } else {
    # If rates are fixed over time.
    rho <- function(t,spec,mu){mu*t - spec*t}
    rho2 <- function(t){mu*t - spec*t}
    tmpfun <- function(t){exp(-rho2(t))*spec}
    nosp <- function(t,spec,mu,n_init){ n_init*(1 + integrate(tmpfun,0,t)$value)}
    # Rough estimate of total species richness expected by Kendall
    minspec <- round(nosp(tmax,spec,mu,n_init))
    Out = array(NA,c(minspec*10,3)); # preallocating the array for species
    Lamds = array(NA,c(minspec*10,1)); # for testing the lambdavalues
    done = 0;
    Out[1:n_init,1]=0;
    tix = 1;
    ntix = n_init+1; #index into first row in Out with no entries yet.
    while (done==0) {
      # When does this lineage go extinct
      if (mu>0){ # if extinction rate is nonzero
        Out[tix,2] = min(Out[tix,1]+rexp(1,mu),tmax)
      } else {
        Out[tix,2] = tmax;
      }
      # Drawing number of fossils for this lineage
      # draw one number from Poisson distribution with mean drawn from a normal distribution with mean lambda
      # and st.dev lambdavar times the duration of the taxon.
      # Out[tix,3] = rpois(1,max(0,rnorm(1,lambda,lambdavar))*(Out[tix,2]-Out[tix,1]))
      # This was before lambdavar = 1 if no variability.
      indsamp = lambda*lambdavar(1);
      Out[tix,3] = rpois(1,max(0,indsamp*(Out[tix,2]-Out[tix,1]))); #
      Lamds[tix] = indsamp/lambda; # Now this is also a scaling factor
      # Drawing waiting times to possible speciation events.
      if (spec>0){
        tmp <- rexp(1e1,spec)
        tmptix = 1;
        # The above draw might be too small;
        while (sum(tmp)<(Out[tix,2]-Out[tix,1])){
          tmp = rexp(10^(tmptix),spec);
          tmptix = tmptix+1; # increase number of draws until end of
          # lineage is traversed
        }
        spectimes = cumsum(tmp)<(Out[tix,2]-Out[tix,1]);
        if (sum(cumsum(tmp)<(Out[tix,2]-Out[tix,1]))>0){

          if ((ntix+sum(spectimes)-1)>nrow(Out)){
            ## need to enlarge OUT
            oldout= Out;
            # Generating new array ~50 % bigger.
            Out = array(NA,c(round(nrow(Out)*1.5),3));
            Out[1:(ntix-1),]=oldout[1:(ntix-1),];
          }
          # if some of these 'speciation times' are inside the actual duration of taxon tix
          Out[ntix:(ntix+sum(spectimes)-1),1]=Out[tix,1]+cumsum(tmp)[spectimes]
          ntix = ntix+sum(spectimes)

        }
      } else {
        # tmp = 0;
        # Do nothing if speciation rate is 0
      }

      tix = tix+1;
      if (tix==(ntix)){
        done=1;
      }

    }
  }
  Out = Out[!is.na(Out[,1]),]; # removin excess NA's
  Lamds = Lamds[!is.na(Lamds)];
  return(list(Out,Lamds))
}



## Functions below are specific for treating the dinosaur data downloaded from PBDB.


## Function to collect occurrances and durations from PBDB output
# This is custom made for dinosaur data, it will only extract occurrences from
# pbdb intervals 112:138
createDataArrs <- function(dinos){

  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Data has [species by interval] with number of occurrences per species in each interval.
  Data = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Times are the durations in a matrix of same size for ease of computation.
  Times = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Times[,ii] <- Bins[ii,3]

  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        # OR after Ladinian (WE include Ladinian here)

        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Data[tix,binow-(Bins[1,4]-1)] <- Data[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  results <- list(Data = Data, Times=Times)
  return(results)
}


createDataArrs_v2 <- function(dinos,removedoubles=FALSE){

  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  # This function adds a third list [Isbird] which is true/false if the taxon is a bird.
  # This is to 'unify' the datamatrices with dinos and theropods with and without birds.
  # if removedoubles is set to TRUE, each occurrence that spans more than one interval is ignored.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Data has [species by interval] with number of occurrences per species in each interval.
  Data = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  Isbird = array(NA,c(length(uniqspec),1))
  # Times are the durations in a matrix of same size for ease of computation.
  Times = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Times[,ii] <- Bins[ii,3]

  }
  ## species by interval matrix
  tix = 1 # loop counter this is for each unique species
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # TESTING TO INCLUDE LADINIAN_ 040615
        if (j1[jj]>j2[jj]) {
          if (removedoubles==TRUE){
            binow <- -100 #make it not count.
          }else{
            countdoubles=countdoubles+1
            bix = seq(j1[jj],j2[jj])
            x = Bins[bix-(Bins[1,4]-1),3]
            binow <- bix[Emprand(x)]
          }
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Data[tix,binow-(Bins[1,4]-1)] <- Data[tix,binow-(Bins[1,4]-1)]+1


        }
      }
    }
    Isbird[tix] = any(dinos[dinos$mid==ii,]$cln==36616)
    Isbird[tix]
    tix <-tix+1
  }

  results <- list(Data = Data, Times=Times,Isbird=Isbird)
  return(results)
}


createDataArrs_v3 <- function(dinos){

  ## Counting occurrences inside each bin for each uniqe species
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$mid[dinos$mra==3])
  # Dats has [species by interval] with number of occurrences per species in each interval.
  Dats = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  uniqnam<-unique(dinos$mna[dinos$mra==3])
  rownames(Dats)<-uniqnam;
  colnames(Dats)<-interval.names;
  # Trying to not defined Data, but do it iteratively to accurately get the ones included. earlier (bf 070815)
  # this code included rows with no occurrences used.
  # The prealloc is needed for computation, how about removing the zeros
  # Times are the durations in a matrix of same size for ease of computation.
  Tims = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Tims[,ii] <- Bins[ii,3]

  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Dats[tix,binow-(Bins[1,4]-1)] <- Dats[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  tmp<-rowSums(Dats)>0;
  Data = Dats[tmp,]
  Times = Tims[tmp,]
  results <- list(Data = Data, Times=Times)
  return(results)
}


createDataArrs_v4 <- function(dinos){

  ## Counting occurrences inside each bin for each uniqe GENUS
  # A uniqe indexlist of occurrences matched to species rank.
  uniqspec <- unique(dinos$gnn)
  uniqspec = uniqspec[!is.na(uniqspec)] # if there's a NA

  # Dats has [species by interval] with number of occurrences per species in each interval.
  Dats = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  dinos[which(dinos$gnn==uniqspec[2]),]$gnl
  uniqnam = array(NA,dim=c(length(uniqspec),1))
  for (ii in 1:length(uniqspec)){
    uniqnam[ii] = unique(dinos[which(dinos$gnn==uniqspec[ii]),]$gnl)
  }
  # uniqnam<-unique(dinos$gnl)
  rownames(Dats)<-uniqnam;
  colnames(Dats)<-interval.names;
  # Trying to not defined Data, but do it iteratively to accurately get the ones included. earlier (bf 070815)
  # this code included rows with no occurrences used.
  # The prealloc is needed for computation, how about removing the zeros
  # Times are the durations in a matrix of same size for ease of computation.
  Tims = matrix(data=0,nrow=length(uniqspec),ncol=nrow(Bins))
  # Inputting the durations in the Times matrix
  for (ii in 1:nrow(Bins)){
    Tims[,ii] <- Bins[ii,3]

  }
  ## species by interval matrix
  tix = 1 # loop counter
  countdoubles = 0; # counting how many occurrences span more than 1 bin
  for (ii in uniqspec) {
    ##  dinos$ein[dinos$mid==ii]-dinos$lin[dinos$mid==ii]
    j1<- dinos$ein[which(dinos$gnn==ii)]#dinos$ein[dinos$mid==ii]
    j2<- dinos$lin[which(dinos$gnn==ii)]#dinos$lin[dinos$mid==ii]
    for (jj in (1:length(j1))) {
      # For each occurrence of this species
      # if (j2[jj]>(Bins[1,4]-1)&j2[jj]<Bins[27,4]){
      if (j2[jj]>(Bins[1,4]-1)&j2[jj]<(Bins[27,4]+1)){
        # if the last interval of this occurrence in this species is before Maastricthian-1
        if (j1[jj]>j2[jj]) {
          countdoubles=countdoubles+1
          bix = seq(j1[jj],j2[jj])
          x = Bins[bix-(Bins[1,4]-1),3]
          binow <- bix[Emprand(x)]
        } else {
          binow <- j1[jj]
        }
        if (binow<Bins[1,4] | binow>Bins[nrow(Bins),4]){
          ignor<- 1
        } else {
          Dats[tix,binow-(Bins[1,4]-1)] <- Dats[tix,binow-(Bins[1,4]-1)]+1
        }
      }
    }
    tix <-tix+1
  }
  tmp<-rowSums(Dats)>0;
  Data = Dats[tmp,]
  Times = Tims[tmp,]
  results <- list(Data = Data, Times=Times)
  return(results)
}

