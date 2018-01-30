# FUNCTIONS:

makedeep<-function(x,y,y_label,quartz=TRUE,main="")
{
  y.max=max(sort(y))
  y.min=0
  young.bin<--8 # Offset for youngest bin label...allows for age < 0 to center label
  plot.min=y.min-.03*(y.max-y.min)
  seg.min=y.min-.07*(y.max-y.min)
  text.min=y.min-.036*(y.max-y.min)
  time.boundaries=c(5.332,23.03,33.9,55.8,65.5,99.6,145.5,161.2,175.6,199.6)
  interval.names=c("Pl","M","O","E","P","UK","LK","UJ","MJ","LJ")
  interval.midpoint=time.boundaries-diff(c(0,time.boundaries))/2
  plot(1,1,xlim=c(max(time.boundaries),0),ylim=c(plot.min,y.max),type="n",xlab="Geologic time (Ma)", ylab=y_label,main=main)
  abline(h=y.min)
  segments(c(time.boundaries,0),y.min,c(time.boundaries,0),seg.min)
  text(interval.midpoint,text.min,labels=interval.names)
  lines(x,y)
}

makedeeppoly<-function(x,y,y_label,upper95,lower95,quartz=TRUE,main="")
{
  y.max=max(c(sort(y),sort(upper95)))
  y.min=0
  young.bin<--8 # Offset for youngest bin label...allows for age < 0 to center label
  plot.min=y.min-.03*(y.max-y.min)
  seg.min=y.min-.07*(y.max-y.min)
  text.min=y.min-.036*(y.max-y.min)
  time.boundaries=c(5.332,23.03,33.9,55.8,65.5,99.6,145.5,161.2,175.6,199.6)
  interval.names=c("Pl","M","O","E","P","UK","LK","UJ","MJ","LJ")
  interval.midpoint=time.boundaries-diff(c(0,time.boundaries))/2
  plot(1,1,xlim=c(max(time.boundaries),0),ylim=c(plot.min,y.max),type="n",xlab="Geologic time (Ma)", ylab=y_label,main=main)
  abline(h=y.min)
  segments(c(time.boundaries,0),y.min,c(time.boundaries,0),seg.min)
  text(interval.midpoint,text.min,labels=interval.names)
  lines(x,y)
  polygon(x=c(x,rev(x)),y=c(upper95,rev(lower95)),col="grey",border=NA)
  points(x,y,type="l")
}

makedeepz<-function(x,y,z,y_label,z_label,quartz=TRUE,main="")
{
	y.max=max(sort(y))
	y.min=0
	z.max=max(sort(z))
	z.min=0
	young.bin<--8 # Offset for youngest bin label...allows for age < 0 to center label
	yplot.min=y.min-.03*(y.max-y.min)
	yseg.min=y.min-.07*(y.max-y.min)
	ytext.min=y.min-.036*(y.max-y.min)
	zplot.min=z.min-.03*(z.max-z.min)
	time.boundaries=c(5.332,23.03,33.9,55.8,65.5,99.6,145.5,161.2,175.6,199.6)
	interval.names=c("Pl","M","O","E","P","UK","LK","UJ","MJ","LJ")
	interval.midpoint=time.boundaries-diff(c(0,time.boundaries))/2
	par(mar=c(5,4,4,4) + 0.1) # Leave space for z axis
	plot(1,1,xlim=c(max(time.boundaries),0),ylim=c(yplot.min,y.max),type="n",xlab="Geologic time (Ma)", ylab=y_label,main=main)
	abline(h=y.min)
	segments(c(time.boundaries,0),y.min,c(time.boundaries,0),yseg.min)
	text(interval.midpoint,ytext.min,labels=interval.names)
	lines(x,y)
	par(new=T)
	plot(1,1,xlim=c(max(time.boundaries),0),ylim=c(zplot.min,z.max), axes=F, type="n",xlab="", ylab="")
	lines(x,z,col="dark grey")
	axis(4, at=pretty(range(z)))
	mtext(z_label, 4, 3)
}

makephan<-function(x,y,y_label,quartz=TRUE,main="")
{
	y.max=max(sort(y))
	y.min=0
	young.bin<--8 # Offset for youngest bin label...allows for age < 0 to center label
	plot.min=y.min-.03*(y.max-y.min)
	seg.min=y.min-.07*(y.max-y.min)
	text.min=y.min-.036*(y.max-y.min)
	time.boundaries=c(23.03,65.5,145.5,199.6,251,299,359.2,416,443.7,488.3,542)
	interval.names=c("Ng","Pg","K","J","Tr","P","C","D","S","O","Cm")
	interval.midpoint=time.boundaries-diff(c(0,time.boundaries))/2
	plot(1,1,xlim=c(max(time.boundaries),0),ylim=c(plot.min,y.max),type="n",xlab="Geologic time (Ma)", ylab=y_label,main=main)
	abline(h=y.min)
	segments(c(time.boundaries,0),y.min,c(time.boundaries,0),seg.min)
	text(interval.midpoint,text.min,labels=interval.names)
	lines(x,y)
}

AICc<-function(model,n)
{
	require(stats)
	require(nlme)
	p<-length(coef(model))
	lk<-AIC(model)-(2*p)
	lk+(2*p*(n/(n-p-1)))
}

linear.model<-function(x,y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1<-nls(y~a*x,start=list(a=(max(y)/max(x))),algorithm="port",lower=list(a=0),control=list(warnOnly=TRUE))
	modelout2<-nls(y~(a*x)+b,start=list(a=(max(y)/max(x)),b=1),algorithm="port",lower=list(a=0,b=0),control=list(warnOnly=TRUE))
	modelout3<-lm(y~x)
	wts<-akaike.wts(c(AICc(modelout1,n=length(x)),AICc(modelout2,n=length(x)),AICc(modelout3,n=length(x))))
	best<-grep(TRUE,max(wts) == wts)
	if(best == 1) model<-modelout1
	if(best == 2) model<-modelout2
	if(best == 3) model<-modelout3
	sefit.model<-std.error(y-predict(model))*1.96
	sdfit.model<-sd(y-predict(model))*1.96
	result<-list(model,sefit.model,sdfit.model)
	names(result)<-c("model","sefit.model","sdfit.model")
	return(result)
}

hyperbolic.model<-function(x,y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1<-NA; modelout2<-NA
	try(modelout1<-nls(y~(a*x)/(b+x),start=list(a=max(y),b=1),algorithm="port",control=list(warnOnly=TRUE)),silent=TRUE)
	try(modelout2<-nls(y~a+((c*x)/(b+x)),start=list(a=1,b=1,c=max(y)),algorithm="port",lower=list(a=0),control=list(warnOnly=TRUE)),silent=TRUE)
	if (is.na(modelout1[1]) == TRUE) modelout1<-nls(y~(max(y)*x)/(b+x),start=list(b=1),algorithm="port",control=list(warnOnly=TRUE))
	if (is.na(modelout2[1]) == TRUE) modelout2<-nls(y~a+((max(y)*x)/(b+x)),start=list(a=1,b=1),algorithm="port",lower=list(a=0),control=list(warnOnly=TRUE))
	wts<-akaike.wts(c(AICc(modelout1,n=length(x)),AICc(modelout2,n=length(x))))
	ifelse(wts[1] >= wts[2],model<-modelout1,model<-modelout2)
	sefit.model<-std.error(y-predict(model))*1.96
	sdfit.model<-sd(y-predict(model))*1.96
	result<-list(model,sefit.model,sdfit.model)
	names(result)<-c("model","sefit.model","sdfit.model")
	return(result)
}

logarithmic.model<-function(x,y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	model<-NA
	try(model<-nls(y~a+log(b+x),start=list(a=1,b=1),algorithm="port",lower=list(a=0,b=0),control=list(warnOnly=TRUE)),silent=TRUE)
	if (is.na(model[1]) == TRUE) model<-nls(y~a+log(1+x),start=list(a=1),algorithm="port",lower=list(a=0),control=list(warnOnly=TRUE))
	sefit.model<-std.error(y-predict(model))*1.96
	sdfit.model<-sd(y-predict(model))*1.96
	result<-list(model,sefit.model,sdfit.model)
	names(result)<-c("model","sefit.model","sdfit.model")
	return(result)
}

exponential.model<-function(x,y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1<-NA; modelout2<-NA
	try(modelout1<-nls(y~a*(1-exp(-c*x)),start=list(a=max(y),c=0.1),algorithm="port",lower=list(a=0,c=0),control=list(warnOnly=TRUE)),silent=TRUE)
	try(modelout2<-nls(y~a-(b*exp(-c*x)),start=list(a=max(y),b=max(y),c=0),algorithm="port",control=list(warnOnly=TRUE)),silent=TRUE)
	if (is.na(modelout1[1]) == TRUE) modelout1<-nls(y~max(y)*(1-exp(-c*x)),start=list(c=0.000001),algorithm="port",lower=list(c=0),control=list(warnOnly=TRUE))
	if (is.na(modelout2[1]) == TRUE) modelout2<-nls(y~max(y)-(max(y)*exp(-c*x)),start=list(c=0.),algorithm="port",control=list(warnOnly=TRUE))
	wts<-akaike.wts(c(AICc(modelout1,n=length(x)),AICc(modelout2,n=length(x))))
	ifelse(wts[1] >= wts[2],model<-modelout1,model<-modelout2)
	sefit.model<-std.error(y-predict(model))*1.96
	sdfit.model<-sd(y-predict(model))*1.96
	result<-list(model,sefit.model,sdfit.model)
	names(result)<-c("model","sefit.model","sdfit.model")
	return(result)
}

sigmoidal.model<-function(x,y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout1<-NA; modelout2<-NA
	try(modelout1<-nls(y~b/(c+exp(-x)),start=list(b=max(y)/100,c=max(y)/10000),algorithm="port",control=list(warnOnly=TRUE)),silent=TRUE)
	try(modelout2<-nls(y~a+(b/(c+exp(-x))),start=list(a=1,b=0.1,c=0.001),algorithm="port",lower=list(a=0),control=list(warnOnly=TRUE)),silent=TRUE)
	if (is.na(modelout1[1]) == TRUE) modelout1<-nls(y~b/((1/max(y))+exp(-x)),start=list(b=1),algorithm="port",control=list(warnOnly=TRUE))
	if (is.na(modelout2[1]) == TRUE) modelout2<-nls(y~a+(max(y)/100)/((max(y)/100)+exp(-x)),start=list(a=1),algorithm="port",lower=list(a=0),control=list(warnOnly=TRUE))
	wts<-akaike.wts(c(AICc(modelout1,n=length(x)),AICc(modelout2,n=length(x))))
	ifelse(wts[1] >= wts[2],model<-modelout1,model<-modelout2)
	sefit.model<-std.error(y-predict(model))*1.96
	sdfit.model<-sd(y-predict(model))*1.96
	result<-list(model,sefit.model,sdfit.model)
	names(result)<-c("model","sefit.model","sdfit.model")
	return(result)
}

polynomial.model<-function(x,y)
{
	require(nlme)
	require(paleoTS)
	require(plotrix)
	modelout2<-lm(y~x+I(x^2))
	modelout3<-lm(y~x+I(x^2)+I(x^3))
	modelout4<-lm(y~x+I(x^2)+I(x^3)+I(x^4))
	wts<-akaike.wts(c(AICc(modelout2,n=length(x)),AICc(modelout3,n=length(x)),AICc(modelout4,n=length(x))))
	best<-grep(TRUE,max(wts) == wts)
	if(best == 1) model<-modelout2
	if(best == 2) model<-modelout3
	if(best == 3) model<-modelout4
	sefit.model<-std.error(y-predict(model))*1.96
	sdfit.model<-sd(y-predict(model))*1.96
	result<-list(model,sefit.model,sdfit.model)
	names(result)<-c("model","sefit.model","sdfit.model")
	return(result)
}

best.model<-function(x,y)
{
	linmod<-linear.model(x,y)$model # Fit linear model
	hypmod<-hyperbolic.model(x,y)$model # Fit hyperbolic model
	logmod<-logarithmic.model(x,y)$model # Fit logarithmic model
	expmod<-exponential.model(x,y)$model # Fit exponential model
	sigmod<-sigmoidal.model(x,y)$model # Fit sigmoidal model
	polmod<-polynomial.model(x,y)$model # Fit polynomial model
	best<-min(c(AICc(linmod,n=length(x)),AICc(hypmod,n=length(x)),AICc(logmod,n=length(x)),AICc(expmod,n=length(x)),AICc(sigmod,n=length(x)),AICc(polmod,n=length(x))))
	if (AICc(linmod,n=length(x)) == best) model<-linmod; sefit.model<-linear.model(x,y)$sefit.model; sdfit.model<-linear.model(x,y)$sdfit.model
	if (AICc(hypmod,n=length(x)) == best) model<-hypmod; sefit.model<-hyperbolic.model(x,y)$sefit.model; sdfit.model<-hyperbolic.model(x,y)$sdfit.model
	if (AICc(logmod,n=length(x)) == best) model<-logmod; sefit.model<-logarithmic.model(x,y)$sefit.model; sdfit.model<-logarithmic.model(x,y)$sdfit.model
	if (AICc(expmod,n=length(x)) == best) model<-expmod; sefit.model<-exponential.model(x,y)$sefit.model; sdfit.model<-exponential.model(x,y)$sdfit.model
	if (AICc(sigmod,n=length(x)) == best) model<-sigmod; sefit.model<-sigmoidal.model(x,y)$sefit.model; sdfit.model<-sigmoidal.model(x,y)$sdfit.model
	if (AICc(polmod,n=length(x)) == best) model<-polmod; sefit.model<-polynomial.model(x,y)$sefit.model; sdfit.model<-polynomial.model(x,y)$sdfit.model
	result<-list(model,sefit.model,sdfit.model)
	names(result)<-c("model","sefit.model","sdfit.model")
	return(result)
}

rockmodel.predictCI<-function(rockmeasure,diversitymeasure,CI=0.95)
{
	x<-sort(rockmeasure)
	y<-sort(diversitymeasure)
	model<-best.model(x,y)$model
	sefit.model<-best.model(x,y)$sefit.model
	sdfit.model<-best.model(x,y)$sdfit.model
	predicted<-predict(model,list(x=rockmeasure))
	selowerCI<-predicted-sefit.model
	seupperCI<-predicted+sefit.model
	sdlowerCI<-predicted-sdfit.model
	sdupperCI<-predicted+sdfit.model
	result<-list(predicted,selowerCI,seupperCI,sdlowerCI,sdupperCI,model)
	names(result)<-c("predicted","selowerCI","seupperCI","sdlowerCI","sdupperCI","model")
	return(result)
}

sphpolyarea<-function(vlat,vlon)
{
	require(fossil)
	splits<-strsplit(unique(paste(vlon,vlat,sep="_")),"_")
	vlat<-vlon<-vector(mode="numeric")
	for (i in 1:length(splits)) {
		vlon[i]<-splits[[i]][1]
		vlat[i]<-splits[[i]][2]
	}
	vlat<-as.numeric(vlat)
	vlon<-as.numeric(vlon)
	points<-chull(vlon,vlat)
	vlon<-vlon[points]
	vlat<-vlat[points]
	sum<-0
	nv<-length(vlat)
	R <- 40041.47/(2 * pi)
	while (nv >= 3) {
		lat1<-vlat[1]
		lat2<-vlat[2]
		lat3<-vlat[3]
		long1<-vlon[1]
		long2<-vlon[2]
		long3<-vlon[3]
		cx <- deg.dist(lat1, long1, lat2, long2)/R
		bx <- deg.dist(lat1, long1, lat3, long3)/R
		ax <- deg.dist(lat2, long2, lat3, long3)/R
		A <- acos((cos(ax) - cos(bx) * cos(cx))/(sin(bx) * sin(cx)))
		B <- acos((cos(bx) - cos(cx) * cos(ax))/(sin(cx) * sin(ax)))
		C <- acos((cos(cx) - cos(ax) * cos(bx))/(sin(ax) * sin(bx)))
		SA <- R^2 * ((A + B + C) - pi)
		sum<-sum+SA
		vlat<-vlat[-2]
		vlon<-vlon[-2]
		nv<-length(vlat)
	}
	return(sum)
}

subsample<-function(pool, ntrials, CI = 0.95)
{
  spgout<-gout<-sout<-matrix(nrow=ntrials,ncol=length(pool))
  stopat<-length(unique(pool)) # No point continuing to pull out species if they have all been sampled
  poolholder<-pool
  topCI<-ceiling((1-((1-CI)/2))*ntrials)
  bottomCI<-max(floor(((1-CI)/2)*ntrials),1)
  for (i in 1:ntrials) {
    pool<-poolholder
    poolsize<-length(pool)
    picked<-vector(mode="character") # Binomial
    pickedgen<-vector(mode="character") # Genus only
    counter<-1
    while (poolsize > 0 && length(picked) < stopat) {
      pickno<-ceiling(runif(1,0,length(pool)))
      picked<-unique(c(picked,pool[pickno]))
      pickedgen<-unique(c(pickedgen,strsplit(pool[pickno]," ")[[1]][1]))
      sout[i,counter]<-length(picked)
      gout[i,counter]<-length(pickedgen)
      spgout[i,counter]<-length(picked)/length(pickedgen)
      pool<-pool[-pickno]
      poolsize<-length(pool)
      counter<-counter+1
    }
    sout[i,grep(TRUE,is.na(sout[i,]))]<-length(picked) # Fill in remaining
    gout[i,grep(TRUE,is.na(gout[i,]))]<-length(pickedgen) # Fill in remaining
    spgout[i,grep(TRUE,is.na(spgout[i,]))]<-length(picked)/length(pickedgen) # Fill in remaining
  }
  smeanCI<-supperCI<-slowerCI<-apply(sout,2,mean)
  gmeanCI<-gupperCI<-glowerCI<-apply(gout,2,mean)
  spgmeanCI<-spgupperCI<-spglowerCI<-apply(spgout,2,mean)
  for (i in 1:length(sout[1,])) {
    supperCI[i]<-sort(sout[,i])[topCI]
    slowerCI[i]<-sort(sout[,i])[bottomCI]
    gupperCI[i]<-sort(gout[,i])[topCI]
    glowerCI[i]<-sort(gout[,i])[bottomCI]
    spgupperCI[i]<-sort(spgout[,i])[topCI]
    spglowerCI[i]<-sort(spgout[,i])[bottomCI]
  }
  sout<-rbind(supperCI,smeanCI,slowerCI)
  gout<-rbind(gupperCI,gmeanCI,glowerCI)
  spgout<-rbind(spgupperCI,spgmeanCI,spglowerCI)
  out<-list(sout,gout,spgout)
  names(out)<-c("sout","gout","spgout")
  return(out)
}

five.models<-function(intable, time, log.input=FALSE)
{
  tstable<-matrix(ncol=4,nrow=length(time))
  if (log.input == FALSE) {
    for (i in 1:length(time)) {
      tstable[i,1]<-ifelse(is.nan(mean(sort(intable[,i]))),0,mean(sort(intable[,i])))
      tstable[i,2]<-ifelse(is.nan(var(sort(intable[,i]))),0,var(sort(intable[,i])))
      tstable[i,3]<-length(sort(intable[,i]))
      tstable[i,4]<-time[i]
    }
  }
  if (log.input == TRUE) {
    for (i in 1:length(time)) {
      tstable[i,1]<-mean(sort(log(intable[,i]))[grep(TRUE,sort(log(intable[,i])) > 0)])
      tstable[i,2]<-var(sort(log(intable[,i]))[grep(TRUE,sort(log(intable[,i])) > 0)])
      tstable[i,3]<-length(sort(log(intable[,i]))[grep(TRUE,sort(log(intable[,i])) > 0)])
      tstable[i,4]<-time[i]
    }
  }
  for (i in length(tstable[,1]):1) {
    kill<-0; a<-tstable[i,1]; b<-tstable[i,2]; c<-tstable[i,3]
    if(is.nan(a) || is.na(a) || a == Inf || a == -Inf || a == 0) kill<-1
    if(is.nan(b) || is.na(b) || b == Inf || b == -Inf || b == 0) kill<-1
    if(is.nan(c) || is.na(c) || c == Inf || c == -Inf || c == 0) kill<-1
    if(kill == 1) tstable<-tstable[-i,]
  }
  tstable<-as.paleoTS(tstable[,1],tstable[,2],tstable[,3],tstable[,4],start.age=max(tstable[,4]))
  threemodels<-fit3models(tstable,silent=T)
  twostases<-fitGpunc(tstable,ng=2,oshare=F,silent=T)
  threestases<-fitGpunc(tstable,ng=3,oshare=F,silent=T)
  results<-cbind(c(threemodels$aic,twostases$AIC,threestases$AIC),c(threemodels$aicc,twostases$AICc,threestases$AICc),round(akaike.wts(c(threemodels$aicc,twostases$AICc,threestases$AICc)),digits=5))
  colnames(results)<-c("AIC","AICc","Akaike Weights"); rownames(results)<-c("Directional trend","Random walk","Stasis","Punctuation with two stases","Punctuation with three stases")
  result<-list(tstable,threemodels,twostases,threestases,results,log.input)
  names(result)<-c("tstable","threemodels","twostases","threestases","results","log.input")
  return(result)
}

plot.stases<-function(fivemodel, rawtable, y_label)
{
  makedeep(x=fivemodel$tstable$tt,y=fivemodel$tstable$mm,y_label=y_label)
  fivemodel$results
  bestmodel<-rownames(fivemodel$results)[grep(TRUE,fivemodel$results[,"Akaike Weights"] == max(fivemodel$results[,"Akaike Weights"]))]
  text(max(fivemodel$tstable$tt),max(fivemodel$tstable$mm),paste("Best model: ",bestmodel,sep=""),pos=4)
  if (bestmodel == "Punctuation with three stases") {
    firststasis<-1:(fivemodel$threestases$shift.start[1]-1)
    secondstasis<-fivemodel$threestases$shift.start[1]:(fivemodel$threestases$shift.start[2]-1)
    thirdstasis<-fivemodel$threestases$shift.start[2]:length(fivemodel$tstable$mm)
  }
  if (bestmodel == "Punctuation with two stases") {
    firststasis<-1:(fivemodel$twostases$shift.start[1]-1)
    secondstasis<-fivemodel$twostases$shift.start[1]:length(fivemodel$tstable$mm)
  }
  SD<-rawtable
  if (fivemodel$log.input == FALSE) {
    for (i in 1:length(SD[1,])) SD[1,i]<-sd(sort(SD[,i]))
  }
  if (fivemodel$log.input == TRUE) {
    for (i in 1:length(SD[1,])) SD[1,i]<-sd(log(sort(SD[,i])))
  }
  SD<-1.96*SD[1,1:length(fivemodel$tstable$mm)]
  if (bestmodel == "Punctuation with three stases") {
    polygon(c(fivemodel$tstable$tt[firststasis],rev(fivemodel$tstable$tt[firststasis])),c(rep(fivemodel$threestases$par[1]+fivemodel$threestases$par[4],length(firststasis)),rep(fivemodel$threestases$par[1]-fivemodel$threestases$par[4],length(firststasis))),col="grey",border=NA)
    polygon(c(fivemodel$tstable$tt[secondstasis],rev(fivemodel$tstable$tt[secondstasis])),c(rep(fivemodel$threestases$par[2]+fivemodel$threestases$par[5],length(secondstasis)),rep(fivemodel$threestases$par[2]-fivemodel$threestases$par[5],length(secondstasis))),col="grey",border=NA)
    polygon(c(fivemodel$tstable$tt[thirdstasis],rev(fivemodel$tstable$tt[thirdstasis])),c(rep(fivemodel$threestases$par[3]+fivemodel$threestases$par[6],length(thirdstasis)),rep(fivemodel$threestases$par[3]-fivemodel$threestases$par[6],length(thirdstasis))),col="grey",border=NA)
  }
  if (bestmodel == "Punctuation with two stases") {
    polygon(c(fivemodel$tstable$tt[firststasis],rev(fivemodel$tstable$tt[firststasis])),c(rep(fivemodel$twostases$par[1]+fivemodel$twostases$par[4],length(firststasis)),rep(fivemodel$twostases$par[1]-fivemodel$twostases$par[4],length(firststasis))),col="grey",border=NA)
    polygon(c(fivemodel$tstable$tt[secondstasis],rev(fivemodel$tstable$tt[secondstasis])),c(rep(fivemodel$twostases$par[2]+fivemodel$twostases$par[5],length(secondstasis)),rep(fivemodel$twostases$par[2]-fivemodel$twostases$par[5],length(secondstasis))),col="grey",border=NA)
  }
  points(x=fivemodel$tstable$tt,y=fivemodel$tstable$mm,type="l")
  points(fivemodel$tstable$tt,fivemodel$tstable$mm,cex=0.5)
  for (i in 1:length(fivemodel$tstable$mm)) lines(c(fivemodel$tstable$tt[i],fivemodel$tstable$tt[i]),c(fivemodel$tstable$mm[i]+SD[i],fivemodel$tstable$mm[i]-SD[i]))
  if (bestmodel == "Punctuation with three stases") {
    points(fivemodel$tstable$tt[firststasis],rep(fivemodel$threestases$par[1],length(firststasis)),type="l")
    points(fivemodel$tstable$tt[secondstasis],rep(fivemodel$threestases$par[2],length(secondstasis)),type="l")
    points(fivemodel$tstable$tt[thirdstasis],rep(fivemodel$threestases$par[3],length(thirdstasis)),type="l")
  }
  if (bestmodel == "Punctuation with two stases") {
    points(fivemodel$tstable$tt[firststasis],rep(fivemodel$twostases$par[1],length(firststasis)),type="l")
    points(fivemodel$tstable$tt[secondstasis],rep(fivemodel$twostases$par[2],length(secondstasis)),type="l")
  }
}

ts.corr<-function(var1, var2)
{
  shared.time<-intersect(var1$tstable$tt,var2$tstable$tt)
  var1.comp<-var1$tstable$mm[match(shared.time,var1$tstable$tt)]
  var2.comp<-var2$tstable$mm[match(shared.time,var2$tstable$tt)]
  raw.corr<-cor.test(var1.comp, var2.comp, method="spearman")
  firstdiff.corr<-cor.test(diff(var1.comp), diff(var2.comp), method="spearman")
  result<-list(raw.corr, firstdiff.corr)
  names(result)<-c("raw.corr", "firstdiff.corr")
  return(result)
}

multi.ts.3comp<-function(var1, var2, var3, var2name, var3name, cuts=NA) # First variable is dependent, second and third are explanatory
{
	require(paleoTS)
	# Create comparable vectors (i.e. where all variables are sampled in bin)
	shared.time<-intersect(intersect(var1$tstable$tt,var2$tstable$tt),var3$tstable$tt)
	if(!is.na(cuts)[1]) {
	  bottom<-max(cuts)
	  top<-min(cuts)
	  shared.time<-shared.time[intersect(grep(TRUE,shared.time >= top),grep(TRUE,shared.time <= bottom))]
	}
	var1.comp<-var1$tstable$mm[match(shared.time,var1$tstable$tt)]
	var2.comp<-var2$tstable$mm[match(shared.time,var2$tstable$tt)]
	var3.comp<-var3$tstable$mm[match(shared.time,var3$tstable$tt)]
	fd.var1.comp<-diff(var1.comp)
	fd.var2.comp<-diff(var2.comp)
	fd.var3.comp<-diff(var3.comp)
	gd.var1.comp<-gen.diff(var1.comp, shared.time)
	gd.var2.comp<-gen.diff(var2.comp, shared.time)
	gd.var3.comp<-gen.diff(var3.comp, shared.time)
	# Fit linear models:
	lm.var2<-lm(var1.comp ~ var2.comp)
	lm.var3<-lm(var1.comp ~ var3.comp)
	lm.var2var3<-lm(var1.comp ~ var2.comp + var3.comp)
	fd.lm.var2<-lm(fd.var1.comp ~ fd.var2.comp)
	fd.lm.var3<-lm(fd.var1.comp ~ fd.var3.comp)
	fd.lm.var2var3<-lm(fd.var1.comp ~ fd.var2.comp + fd.var3.comp)
	gd.lm.var2<-lm(gd.var1.comp ~ gd.var2.comp)
	gd.lm.var3<-lm(gd.var1.comp ~ gd.var3.comp)
	gd.lm.var2var3<-lm(gd.var1.comp ~ gd.var2.comp + gd.var3.comp)
	# Get proportion explained by each model, plus unexplained:
	var2.rsq<-summary(lm.var2var3)$r.squared-summary(lm.var3)$r.squared
	var3.rsq<-summary(lm.var2var3)$r.squared-summary(lm.var2)$r.squared
	var2var3.rsq<-summary(lm.var2)$r.squared-var2.rsq
	ue.rsq<-1-(var2.rsq+var3.rsq+var2var3.rsq)
	fd.var2.rsq<-summary(fd.lm.var2var3)$r.squared-summary(fd.lm.var3)$r.squared
	fd.var3.rsq<-summary(fd.lm.var2var3)$r.squared-summary(fd.lm.var2)$r.squared
	fd.var2var3.rsq<-summary(fd.lm.var2)$r.squared-fd.var2.rsq
	fd.ue.rsq<-1-(fd.var2.rsq+fd.var3.rsq+fd.var2var3.rsq)
	gd.var2.rsq<-summary(gd.lm.var2var3)$r.squared-summary(gd.lm.var3)$r.squared
	gd.var3.rsq<-summary(gd.lm.var2var3)$r.squared-summary(gd.lm.var2)$r.squared
	gd.var2var3.rsq<-summary(gd.lm.var2)$r.squared-gd.var2.rsq
	gd.ue.rsq<-1-(gd.var2.rsq+gd.var3.rsq+gd.var2var3.rsq)
	# Get AIC and AICc for each model plus Akaike weights:
	AIC.var2<-AIC(lm.var2)
	AIC.var3<-AIC(lm.var3)
	AIC.var2var3<-AIC(lm.var2var3)
	fd.AIC.var2<-AIC(fd.lm.var2)
	fd.AIC.var3<-AIC(fd.lm.var3)
	fd.AIC.var2var3<-AIC(fd.lm.var2var3)
	gd.AIC.var2<-AIC(gd.lm.var2)
	gd.AIC.var3<-AIC(gd.lm.var3)
	gd.AIC.var2var3<-AIC(gd.lm.var2var3)
	AICc.var2<-AICc(lm.var2, length(var1.comp))
	AICc.var3<-AICc(lm.var3, length(var1.comp))
	AICc.var2var3<-AICc(lm.var2var3, length(var1.comp))
	fd.AICc.var2<-AICc(fd.lm.var2, length(fd.var1.comp))
	fd.AICc.var3<-AICc(fd.lm.var3, length(fd.var1.comp))
	fd.AICc.var2var3<-AICc(fd.lm.var2var3, length(fd.var1.comp))
	gd.AICc.var2<-AICc(gd.lm.var2, length(gd.var1.comp))
	gd.AICc.var3<-AICc(gd.lm.var3, length(gd.var1.comp))
	gd.AICc.var2var3<-AICc(gd.lm.var2var3, length(gd.var1.comp))
	wts<-akaike.wts(c(AICc.var2, AICc.var3, AICc.var2var3))
	fd.wts<-akaike.wts(c(fd.AICc.var2, fd.AICc.var3, fd.AICc.var2var3))
	gd.wts<-akaike.wts(c(gd.AICc.var2, gd.AICc.var3, gd.AICc.var2var3))
	# Write formulae for models:
	fm.var2<-paste("(",round(coef(lm.var2)[2],2)," x ",var2name,") + ",round(coef(lm.var2)[1],2),sep="")
	fm.var3<-paste("(",round(coef(lm.var3)[2],2)," x ",var3name,") + ",round(coef(lm.var3)[1],2),sep="")
	fm.var2var3<-paste("(",round(coef(lm.var2var3)[2],2)," x ",var2name,") + (",round(coef(lm.var2var3)[3],2)," x ",var3name,") + ",round(coef(lm.var2var3)[1],2),sep="")
	fd.fm.var2<-paste("(",round(coef(fd.lm.var2)[2],2)," x ",var2name,") + ",round(coef(fd.lm.var2)[1],2),sep="")
	fd.fm.var3<-paste("(",round(coef(fd.lm.var3)[2],2)," x ",var3name,") + ",round(coef(fd.lm.var3)[1],2),sep="")
	fd.fm.var2var3<-paste("(",round(coef(fd.lm.var2var3)[2],2)," x ",var2name,") + (",round(coef(fd.lm.var2var3)[3],2)," x ",var3name,") + ",round(coef(fd.lm.var2var3)[1],2),sep="")
	gd.fm.var2<-paste("(",round(coef(gd.lm.var2)[2],2)," x ",var2name,") + ",round(coef(gd.lm.var2)[1],2),sep="")
	gd.fm.var3<-paste("(",round(coef(gd.lm.var3)[2],2)," x ",var3name,") + ",round(coef(gd.lm.var3)[1],2),sep="")
	gd.fm.var2var3<-paste("(",round(coef(gd.lm.var2var3)[2],2)," x ",var2name,") + (",round(coef(gd.lm.var2var3)[3],2)," x ",var3name,") + ",round(coef(gd.lm.var2var3)[1],2),sep="")
	# Make table:
	models<-c(var2name, var3name, paste(var2name," + ",var3name,sep=""), "Unexplained")
	results<-cbind(models, c(fm.var2, fm.var3, fm.var2var3, NA), c(fd.fm.var2, fd.fm.var3, fd.fm.var2var3, NA), c(gd.fm.var2, gd.fm.var3, gd.fm.var2var3, NA), c(var2.rsq, var3.rsq, var2var3.rsq, ue.rsq), c(fd.var2.rsq, fd.var3.rsq, fd.var2var3.rsq, fd.ue.rsq), c(gd.var2.rsq, gd.var3.rsq, gd.var2var3.rsq, gd.ue.rsq), c(AIC.var2, AIC.var3, AIC.var2var3, NA), c(fd.AIC.var2, fd.AIC.var3, fd.AIC.var2var3, NA), c(gd.AIC.var2, gd.AIC.var3, gd.AIC.var2var3, NA), c(AICc.var2, AICc.var3, AICc.var2var3, NA), c(fd.AICc.var2, fd.AICc.var3, fd.AICc.var2var3, NA), c(gd.AICc.var2, gd.AICc.var3, gd.AICc.var2var3, NA), c(wts, NA), c(fd.wts, NA), c(gd.wts, NA))
	colnames(results)<-c("Explanatory variable", "Model formula (raw data)", "Model formula (first differences)", "Model formula (generalised differences)", "Proportion variance explained (raw data)", "Proportion variance explained (first differences)", "Proportion variance explained (generalised differences)", "AIC (raw data)", "AIC (first differences)", "AIC (generalised differences)", "AICc (raw data)", "AICc (first differences)", "AICc (generalised differences)", "Akaike weight (raw data)", "Akaike weight (first differences)", "Akaike weight (generalised differences)")
	result<-list(shared.time, var1.comp, var2.comp, var3.comp, fd.var1.comp, fd.var2.comp, fd.var3.comp, gd.var1.comp, gd.var2.comp, gd.var3.comp, lm.var2, lm.var3, lm.var2var3, fd.lm.var2, fd.lm.var3, fd.lm.var2var3, gd.lm.var2, gd.lm.var3, gd.lm.var2var3, results)
	names(result)<-c("shared.time", "var1.comp", "var2.comp", "var3.comp", "fd.var1.comp", "fd.var2.comp", "fd.var3.comp", "gd.var1.comp", "gd.var2.comp", "gd.var3.comp", "lm.var2", "lm.var3", "lm.var2var3", "fd.lm.var2", "fd.lm.var3", "fd.lm.var2var3", "gd.lm.var2", "gd.lm.var3", "gd.lm.var2var3", "results")
	return(result)
}

moving.average<-function(ts.holder, mav)
{
  ts<-vector(mode="numeric")
  for (i in 1:length(ts.holder[1,])) ts[i]<-mean(sort(ts.holder[,i]))
  trimtb<-floor(mav/2)
  mav.vector<-vector(mode="numeric",length=length(ts)-(2*trimtb))
  for (i in (trimtb+1):(length(ts)-trimtb)) mav.vector[(i-trimtb)]<-mean(sort(ts[(i-trimtb):(i+trimtb)]))
  mav.vector<-as.numeric(gsub(NaN,NA,mav.vector))
  mav.vector<-c(rep(NA,trimtb),mav.vector,rep(NA,trimtb))
  mav.vector
}

moving.averagets<-function(ts, mav)
{
	trimtb<-floor(mav/2)
	mav.vector<-vector(mode="numeric",length=length(ts)-(2*trimtb))
	for (i in (trimtb+1):(length(ts)-trimtb)) mav.vector[(i-trimtb)]<-mean(sort(ts[(i-trimtb):(i+trimtb)]))
	mav.vector<-as.numeric(gsub(NaN,NA,mav.vector))
	mav.vector<-c(rep(NA,trimtb),mav.vector,rep(NA,trimtb))
	mav.vector
}

corr.plot<-function(x.holder, y.holder, time, log.x=FALSE, log.y=FALSE, fd=FALSE, gd=FALSE, do.mav=FALSE, mav, xlab="", ylab="", main="")
{
	y<-x<-vector(mode="numeric")
	for (i in 1:length(x.holder[1,])) x[i]<-mean(sort(x.holder[,i]))
	for (i in 1:length(y.holder[1,])) y[i]<-mean(sort(y.holder[,i]))
	if(log.x == FALSE) x1<-as.numeric(gsub(NaN,NA,x))
	if(log.x == TRUE) x1<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,x)))))
	if(log.y == FALSE) y1<-as.numeric(gsub(NaN,NA,y))
	if(log.y == TRUE) y1<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,y)))))
	if(fd == TRUE) {
		x1<-diff(x1)
		y1<-diff(y1)
	}
	if(gd == TRUE) {
		x1<-gen.diff(x1, time)
		y1<-gen.diff(y1, time)
	}	
	if(do.mav == TRUE) {
		trimtb<-floor(mav/2)
		mav.x1<-vector(mode="numeric",length=length(x1)-(2*trimtb))
		mav.y1<-vector(mode="numeric",length=length(y1)-(2*trimtb))
		for (i in (trimtb+1):(length(x1)-trimtb)) mav.x1[(i-trimtb)]<-mean(sort(x1[(i-trimtb):(i+trimtb)]))
		for (i in (trimtb+1):(length(y1)-trimtb)) mav.y1[(i-trimtb)]<-mean(sort(y1[(i-trimtb):(i+trimtb)]))
		mav.x1<-as.numeric(gsub(NaN,NA,mav.x1))
		mav.y1<-as.numeric(gsub(NaN,NA,mav.y1))
		mav.x1<-c(rep(NA,trimtb),mav.x1,rep(NA,trimtb))
		mav.y1<-c(rep(NA,trimtb),mav.y1,rep(NA,trimtb))
		x1<-x1-mav.x1
		y1<-y1-mav.y1
	}
	plot(x1,y1,xlab=xlab,ylab=ylab,main=main)
	text(x=min(sort(x1)),y=max(sort(y1)),labels=paste("Spearman rho: ",round(cor.test(x1,y1,method="spearman")$estimate,2)," (p = ",format(round(cor.test(x1,y1,method="spearman")$p.value, 4), nsmall = 4),")",sep="",collapse=""),adj=c(0,1))
	abline(lsfit(x1,y1)$coefficients[1],lsfit(x1,y1)$coefficients[2])
}

ts.corr.plot<-function(var1.holder, var2.holder, var3.holder, log.var1=FALSE, log.var2=FALSE, log.var3=FALSE, var1.name, var2.name, var3.name, time)
{
  var1<-var2<-var3<-vector(mode="numeric")
  for (i in 1:length(var1.holder[1,])) var1[i]<-mean(sort(var1.holder[,i]))
  for (i in 1:length(var2.holder[1,])) var2[i]<-mean(sort(var2.holder[,i]))
  for (i in 1:length(var3.holder[1,])) var3[i]<-mean(sort(var3.holder[,i]))
  if(log.var1 == FALSE) var1<-as.numeric(gsub(NaN,NA,var1))
  if(log.var1 == TRUE) var1<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,var1)))))
  if(log.var2 == FALSE) var2<-as.numeric(gsub(NaN,NA,var2))
  if(log.var2 == TRUE) var2<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,var2)))))
  if(log.var3 == FALSE) var3<-as.numeric(gsub(NaN,NA,var3))
  if(log.var3 == TRUE) var3<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,var3)))))
  par(mfrow=c(2,2))
  plot(var1, var2, cex=0.5, xlab=var1.name, ylab=var2.name)
  abline(a=lsfit(var1,var2)$coefficients[1],b=lsfit(var1,var2)$coefficients[2])
  text(x=min(sort(var1)), y=max(sort(var2)), paste("rho = ",round(cor.test(var1, var2, method="spearman")$estimate,2)), adj=c(0,1))
  plot(gen.diff(var1,time), gen.diff(var2,time), cex=0.5, xlab=var1.name, ylab=var2.name)
  abline(a=lsfit(gen.diff(var1,time),gen.diff(var2,time))$coefficients[1],b=lsfit(gen.diff(var1,time),gen.diff(var2,time))$coefficients[2])
  lines(x=c(0,0),y=c(-10,10),col="grey",lty=2)
  lines(x=c(-10,10),y=c(0,0),col="grey",lty=2)
  text(x=min(sort(gen.diff(var1,time))), y=max(sort(gen.diff(var2,time))), paste("rho = ",round(cor.test(gen.diff(var1,time), gen.diff(var2,time), method="spearman")$estimate,2)), adj=c(0,1))
  plot(var1, var3, cex=0.5, xlab=var1.name, ylab=var3.name)
  abline(a=lsfit(var1,var3)$coefficients[1],b=lsfit(var1,var3)$coefficients[2])
  text(x=min(sort(var1)), y=max(sort(var3)), paste("rho = ",round(cor.test(var1, var3, method="spearman")$estimate,2)), adj=c(0,1))
  plot(gen.diff(var1,time), gen.diff(var3,time), cex=0.5, xlab=var1.name, ylab=var3.name)
  abline(a=lsfit(gen.diff(var1,time),gen.diff(var3,time))$coefficients[1],b=lsfit(gen.diff(var1,time),gen.diff(var3,time))$coefficients[2])
  lines(x=c(0,0),y=c(-10,10),col="grey",lty=2)
  lines(x=c(-10,10),y=c(0,0),col="grey",lty=2)
  text(x=min(sort(gen.diff(var1,time))), y=max(sort(gen.diff(var3,time))), paste("rho = ",round(cor.test(gen.diff(var1,time), gen.diff(var3,time), method="spearman")$estimate,2)), adj=c(0,1))
}

ts.diffs.plot<-function(var1.holder, var2.holder, var3.holder, log.var1=FALSE, log.var2=FALSE, log.var3=FALSE, diff.type="fd", mav, time)
{
	var1<-var2<-var3<-vector(mode="numeric")
	for (i in 1:length(var1.holder[1,])) var1[i]<-mean(sort(var1.holder[,i]))
	for (i in 1:length(var2.holder[1,])) var2[i]<-mean(sort(var2.holder[,i]))
	for (i in 1:length(var3.holder[1,])) var3[i]<-mean(sort(var3.holder[,i]))
	if(log.var1 == FALSE) var1<-as.numeric(gsub(NaN,NA,var1))
	if(log.var1 == TRUE) var1<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,var1)))))
	if(log.var2 == FALSE) var2<-as.numeric(gsub(NaN,NA,var2))
	if(log.var2 == TRUE) var2<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,var2)))))
	if(log.var3 == FALSE) var3<-as.numeric(gsub(NaN,NA,var3))
	if(log.var3 == TRUE) var3<-as.numeric(gsub(-Inf,NA,log(as.numeric(gsub(NaN,NA,var3)))))
	if(diff.type == "fd") {
		var1<-diff(var1)
		var2<-diff(var2)
		var3<-diff(var3)
	}
	if(diff.type == "gd") {
		var1<-gen.diff(var1, time)
		var2<-gen.diff(var2, time)
		var3<-gen.diff(var3, time)
	}
	if(diff.type == "mav") {
		trimtb<-floor(mav/2)
		mav.var1<-vector(mode="numeric",length=length(var1)-(2*trimtb))
		mav.var2<-vector(mode="numeric",length=length(var2)-(2*trimtb))
		mav.var3<-vector(mode="numeric",length=length(var3)-(2*trimtb))
		for (i in (trimtb+1):(length(var1)-trimtb)) mav.var1[(i-trimtb)]<-mean(sort(var1[(i-trimtb):(i+trimtb)]))
		for (i in (trimtb+1):(length(var2)-trimtb)) mav.var2[(i-trimtb)]<-mean(sort(var2[(i-trimtb):(i+trimtb)]))
		for (i in (trimtb+1):(length(var3)-trimtb)) mav.var3[(i-trimtb)]<-mean(sort(var3[(i-trimtb):(i+trimtb)]))
		mav.var1<-as.numeric(gsub(NaN,NA,mav.var1))
		mav.var2<-as.numeric(gsub(NaN,NA,mav.var2))
		mav.var3<-as.numeric(gsub(NaN,NA,mav.var3))
		mav.var1<-c(rep(NA,trimtb),mav.var1,rep(NA,trimtb))
		mav.var2<-c(rep(NA,trimtb),mav.var2,rep(NA,trimtb))
		mav.var3<-c(rep(NA,trimtb),mav.var3,rep(NA,trimtb))
		var1<-var1-mav.var1
		var2<-var2-mav.var2
		var3<-var3-mav.var3
	}
	par(mfrow=c(2,1))
	if(diff.type == "fd") {
		time<-time[1:(length(time)-1)]+diff(time)/2
		plot(time, var1, xlim=c(max(time),0), ylim=c(-max(sqrt(sort(c(var1,var2))^2)),max(sqrt(sort(c(var1,var2))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="First diff.")
		points(time, var2, type="l")
	}
	if(diff.type == "gd") {
		time<-time[1:(length(time)-1)]+diff(time)/2
		plot(time, var1, xlim=c(max(time),0), ylim=c(-max(sqrt(sort(c(var1,var2))^2)),max(sqrt(sort(c(var1,var2))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Gen. diff.")
		points(time, var2, type="l")
	}
	if(diff.type == "mav") {
		plot(time, var1, xlim=c(max(time),0), ylim=c(-max(sqrt(sort(c(var1,var2))^2)),max(sqrt(sort(c(var1,var2))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Moving average detrended")
		points(time, var2, type="l")
	}
	lines(x=c(max(time)+20,-20), y=c(0,0), col="grey", lty=2)
	if(diff.type == "fd") plot(time, var1, xlim=c(max(time),0), ylim=c(-max(sqrt(sort(c(var1,var3))^2)),max(sqrt(sort(c(var1,var3))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="First diff.")
	if(diff.type == "gd") plot(time, var1, xlim=c(max(time),0), ylim=c(-max(sqrt(sort(c(var1,var3))^2)),max(sqrt(sort(c(var1,var3))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Gen. diff.")
	if(diff.type == "mav") plot(time, var1, xlim=c(max(time),0), ylim=c(-max(sqrt(sort(c(var1,var3))^2)),max(sqrt(sort(c(var1,var3))^2))), type="l", col="grey", xlab="Time (Ma)", ylab="Moving average detrended")
	points(time, var3, type="l")
	lines(x=c(max(time)+20,-20), y=c(0,0), col="grey", lty=2)
}

parse.nexustotnt<-function(file)
{
	X<-scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
	# Remove trees:
	deleteline<-grep("BEGIN TREES",X)-1
	X<-X[1:deleteline]
	# Replace #NEXUS with xread:
	nexusline<-grep("NEXUS",X)
	X[nexusline]<-"xread"
	# Replace [!...] with '...':
	textline<-grep("[",X,fixed=TRUE)
	X[textline]<-gsub("[!","'",X[textline],fixed=TRUE)
	X[textline]<-gsub("]","'",X[textline],fixed=TRUE)
	# Replace DATA block with char, ntax:
	ntaxline<-grep("NTAX",X)
	X[ntaxline]<-gsub("\tDIMENSIONS  NTAX=","",X[ntaxline])
	X[ntaxline]<-gsub("NCHAR=","",X[ntaxline])
	X[ntaxline]<-gsub(";","",X[ntaxline])
	X[ntaxline]<-paste(strsplit(X[ntaxline]," ")[[1]][2],strsplit(X[ntaxline]," ")[[1]][1],sep=" ")
	# Delete now redundant starting lines:
	X<-X[-(3:(ntaxline-1))]
	formatline<-grep("FORMAT",X)
	matrixline<-grep("MATRIX",X)
	semicolonline<-grep(";",X,fixed=TRUE)[grep(TRUE,grep(";",X,fixed=TRUE) > matrixline)][1]
	if(length(grep("(",X[(matrixline+1):(semicolonline-1)],fixed=TRUE)) > 0) X[(matrixline+1):(semicolonline-1)]<-gsub("(","[",X[(matrixline+1):(semicolonline-1)],fixed=TRUE) # Replace polymorphic ( with [
	if(length(grep(")",X[(matrixline+1):(semicolonline-1)],fixed=TRUE)) > 0) X[(matrixline+1):(semicolonline-1)]<-gsub(")","]",X[(matrixline+1):(semicolonline-1)],fixed=TRUE) # Replace polymorphic ) with ]
	missingchar<-strsplit(strsplit(X[formatline],"MISSING=")[[1]][2],"")[[1]][1] # Find character used for missing states
	gapchar<-strsplit(strsplit(X[formatline],"GAP=")[[1]][2],"")[[1]][1] # Find character used for gaps
	if (missingchar != "?") X[(matrixline+1):(semicolonline-1)]<-gsub(missingchar,"?",X[(matrixline+1):(semicolonline-1)],fixed=TRUE) # Replace with ? if not ?
	if (gapchar != "-") X[(matrixline+1):(semicolonline-1)]<-gsub(gapchar,"-",X[(matrixline+1):(semicolonline-1)],fixed=TRUE) # Replace with - if not -
	X<-X[-matrixline]
	X<-X[-formatline]
	# Re-find end of matrix:
	semicolonline<-grep(";",X,fixed=TRUE)[1]
	# Search for ordered character list:
	if (length(grep("TYPESET",X)) > 0) {
		orderingline<-grep("TYPESET",X)
		X[orderingline]<-gsub("\tTYPESET * UNTITLED  = ","",X[orderingline],fixed=TRUE) # Clean up
		X[orderingline]<-gsub(";","",X[orderingline],fixed=TRUE) # Clean up
		unordered<-strsplit(X[orderingline],", ")[[1]][1] # Get unordered list
		ordered<-strsplit(X[orderingline],", ")[[1]][2] # Get ordered list
		unordered<-gsub("unord: ","",unordered,fixed=TRUE) # Clean up
		ordered<-gsub("ord: ","",ordered,fixed=TRUE) # Clean up
		if(length(grep("-",unordered)) > 0) {
			while(length(grep("-",unordered)) > 0) {
				unordered<-strsplit(unordered," ")[[1]]
				unordered[grep("-",unordered)[1]]<-paste(c(strsplit(unordered[grep("-",unordered)[1]],"-")[[1]][1]:strsplit(unordered[grep("-",unordered)[1]],"-")[[1]][2]),collapse=" ")
				unordered<-paste(unordered,collapse=" ")
			}
		}
		if(length(grep("-",ordered)) > 0) {
			while(length(grep("-",ordered)) > 0) {
				ordered<-strsplit(ordered," ")[[1]]
				ordered[grep("-",ordered)[1]]<-paste(c(strsplit(ordered[grep("-",ordered)[1]],"-")[[1]][1]:strsplit(ordered[grep("-",ordered)[1]],"-")[[1]][2]),collapse=" ")
				ordered<-paste(ordered,collapse=" ")
			}
		}
		orderline<-paste("ccode +",ordered," -",unordered,";",sep="",collapse="")
		endline<-c(orderline,"proc/;")
	} else {
		orderline<-paste("ccode -",paste(1:strsplit(X[3]," ")[[1]][1],collapse=" "),";",sep="",collapse="")
		endline<-c(orderline,"proc/;")
	}
	X<-X[-((semicolonline+1):length(X))]
	X<-c(X,endline)
	return(X)
}

path.lengths<-function(phy)
{
	ntips<-length(phy$tip.label)
	pathlengths<-vector(mode="numeric")
	for (i in 1:ntips) {
		taxon<-i
		pathedges<-vector(mode="numeric")
		pathedges[1]<-grep(TRUE,phy$edge[,2] == i)
		while (phy$edge[pathedges[length(pathedges)],1] != (ntips+1)) {
			i<-grep(TRUE,phy$edge[,2] == phy$edge[pathedges[length(pathedges)],1])
			pathedges<-c(pathedges,i)
		}
		pathlengths[taxon]<-sum(phy$edge.length[pathedges])
	}
	names(pathlengths)<-phy$tip.label
	return(pathlengths)
}

find.descendants<-function(n,tree)
{
	n<-as.vector(n)
	ancs<-vector(mode="numeric")
	desc<-vector(mode="numeric")
	while (length(n) > 0) {
		for (i in 1:length(n)) ancs<-c(tree$edge[grep(TRUE,tree$edge[,1] == n[i]),2],ancs)
		n<-vector(mode="numeric")
		for (i in length(ancs):1) {
			if (ancs[i] <= Ntip(tree)) desc<-c(desc,ancs[i])
			if (ancs[i] > Ntip(tree)) n<-c(n,ancs[i])
			ancs<-ancs[-i]
		}
	}
	return(desc)
}

date.nodes<-function(tree, root.age)
{
	rootnode<-Ntip(tree)+1
	paths<-matrix(NA,nrow=Ntip(tree),ncol=(Ntip(tree)+Nnode(tree)-1))
	paths[1,]<-c(1:Ntip(tree),(Ntip(tree)+2):(Ntip(tree)+Nnode(tree)))
	for (i in 1:length(paths[1,])) {
		j<-1
		currentnode<-paths[j,i]
		while (currentnode != rootnode) {
			currentnode<-paths[j+1,i]<-tree$edge[match(currentnode,tree$edge[,2]),1]
			j<-j+1
		}
	}
	nodeages<-vector(mode="numeric",length=Ntip(tree)+Nnode(tree))
	for (i in 1:length(paths[1,])) {
		nodeages[paths[1,i]]<-sum(tree$edge.length[match(paths[1:(match(NA,paths[,i])-2),i],tree$edge[,2])])
	}
	nodeages<-root.age-nodeages
	return(nodeages)
}

pat.dist.phylo<-function(tree, comp, nchar)
{
	require(ape)
	nt<-Ntip(tree)
	cchar<-nchar-comp[tree$tip.label,1]
	for (i in 1:length(tree$edge[,1])) {
		if (tree$edge[i,2] <= max(nt))  tree$edge.length[i]<-tree$edge.length[i]/cchar[tree$edge[i,2]]
		if (tree$edge[i,2] > max(nt))  tree$edge.length[i]<-tree$edge.length[i]/nchar
	}
	return(tree)
}

date.phylo<-function(tree, ages, rlen, method="equal", ptree)
{
	require(ape)	
	nodes<-c(Ntip(tree)+1):(Nnode(tree)+Ntip(tree))
	nodeages<-nodes
	for (i in 1:length(nodes)) nodeages[i]<-max(ages[tree$tip.label[find.descendants(nodes[i],tree)],1])
	allages<-as.vector(c(ages[tree$tip.label,1],nodeages))
	allages[Ntip(tree)+1]<-allages[Ntip(tree)+1]+rlen
	ttree<-tree
	for (i in 1:length(ttree$edge[,1])) ttree$edge.length[i]<-allages[ttree$edge[i,1]]-allages[ttree$edge[i,2]]
	if (method != "standard") {
		btc<-grep(TRUE,ttree$edge.length == 0) # Branches to check
		bro<-cbind(ttree$edge[,2],node.depth(ttree)[ttree$edge[,2]]) # Get branch list
		rownames(bro)<-1:length(bro[,1])
		bro<-bro[btc,] # Just zero-length branches
		bro<-bro[order(bro[,2]),]
		bro<-as.numeric(rownames(bro))
		for (i in 1:length(bro)) {
			i<-bro[i]
			if (ttree$edge.length[i] == 0) {
				brs<-vector(mode="numeric")
				brs[length(brs)+1]<-ttree$edge[i,2] # Find actual branch
				while (ttree$edge.length[match(ttree$edge[i,1],ttree$edge[,2])] == 0 && match(ttree$edge[i,1],ttree$edge[,2]) != 1) { # Find other preceding branches of zero length and stop at root!
					i<-match(ttree$edge[i,1],ttree$edge[,2])
					brs[length(brs)+1]<-ttree$edge[i,2]
				}
				brs[length(brs)+1]<-ttree$edge[match(brs[length(brs)],ttree$edge[,2]),1] # Add final branch
				tt<-sum(ttree$edge.length[match(brs,ttree$edge[,2])]) # Amount of time to be shared
				if (method == "equal") { # Sharing proportions if method is equal
					ppt<-1:length(brs)
					for (j in 1:length(ppt)) ppt[j]<-1/length(brs)
				}
				if (method == "ruta") { # Sharing proportions if method is Ruta
					if (sum(ptree$edge.length[match(brs,ttree$edge[,2])]) > 0) ppt<-ptree$edge.length[match(brs,ttree$edge[,2])]/sum(ptree$edge.length[match(brs,ttree$edge[,2])])
					if (sum(ptree$edge.length[match(brs,ttree$edge[,2])]) == 0) {
						ppt<-1:length(brs)
						for (j in 1:length(ppt)) ppt[j]<-1/length(brs)
					}
				}
				tc<-ppt*tt # Actual branch lengths (proportion of total length shared)
				for (j in 2:length(tc)) tc[j]<-tc[j]+tc[j-1] # Cumulative branch ages
				allages[ttree$edge[match(brs,ttree$edge[,2]),1]][1:length(tc)-1]<-allages[ttree$edge[match(brs,ttree$edge[,2]),1]][1:length(tc)-1]+tc[1:length(tc)-1] # Update ages of nodes/tips
				ttree$edge.length<-allages[tree$edge[,1]]-allages[tree$edge[,2]] # Update actual branch lengths
			}
		}
	}
	return(ttree)
}

SCM<-function(names, ages, tbins)
{
	out<-matrix(NA, nrow=length(unique(names)), ncol=length(tbins)-1)
	rownames(out)<-sort(unique(names))
	for (i in 1:length(out[,1])) {
		dates<-sort(ages[grep(TRUE,names == rownames(out)[i])])
		for(j in 1:length(dates)) {
			out[i,intersect(grep(TRUE,dates[j] <= tbins[1:(length(tbins)-1)]),grep(TRUE,dates[j] > tbins[2:length(tbins)]))]<-1
		}
		start<-min(grep(TRUE,out[i,] == 1))
		stop<-max(grep(TRUE,out[i,] == 1))
		out[i,(start+grep(TRUE,is.na(out[i,start:stop]))-1)]<-0
	}
	out.adj<-out
	for (i in length(out.adj[,1]):1) {
		out.adj[i,min(grep(TRUE,out.adj[i,] == 1))]<-NA
		out.adj[i,max(grep(TRUE,out[i,] == 1))]<-NA
		if(length(sort(out.adj[i,])) == 0) out.adj<-out.adj[-i,]
	}
	SCM.taxa.adj<-SCM.tbins.adj<-SCM.taxa<-SCM.tbins<-vector(mode="numeric")
	for(i in 1:length(out[1,])) {
		SCM.tbins[i]<-mean(sort(out[,i]))
	}
	for(i in 1:length(out[,1])) {
		SCM.taxa[i]<-mean(sort(out[i,]))
	}
	for(i in 1:length(out.adj[1,])) {
		SCM.tbins.adj[i]<-mean(sort(out.adj[,i]))
	}
	for(i in 1:length(out.adj[,1])) {
		SCM.taxa.adj[i]<-mean(sort(out.adj[i,]))
	}
	names(SCM.taxa)<-rownames(out)
	names(SCM.taxa.adj)<-rownames(out.adj)
	tbin.midpoints<-((tbins[1:(length(tbins)-1)]-tbins[2:length(tbins)])/2)+tbins[2:length(tbins)]
	result<-list(out, SCM.tbins, SCM.taxa, SCM.tbins.adj, SCM.taxa.adj, tbin.midpoints)
	names(result)<-c("out", "SCM.tbins", "SCM.taxa", "SCM.tbins.adj", "SCM.taxa.adj", "tbin.midpoints")
	return(result)
}

make.phan<-function(x,y,y_label,quartz=TRUE,main="")
{
	y.max<-max(sort(y))
	y.min<-0
	young.bin<--8 # Offset for youngest bin label...allows for age < 0 to center label
	plot.min<-y.min-.03*(y.max-y.min)
	seg.min<-y.min-.07*(y.max-y.min)
	text.min<-y.min-.036*(y.max-y.min)
	time.boundaries<-c(23.03,65.5,145.5,199.6,251,299,359.2,416,443.7)
	interval.names<-c("Ng","Pg","K","J","Tr","P","C","D","S")
	interval.midpoint<-time.boundaries-diff(c(0,time.boundaries))/2
	plot(1,1,xlim=c(max(time.boundaries),0),ylim=c(plot.min,y.max),type="n",xlab="Geologic time (Ma)", ylab=y_label,main=main)
	abline(h=y.min)
	segments(c(time.boundaries,0),y.min,c(time.boundaries,0),seg.min)
	text(interval.midpoint,text.min,labels=interval.names)
	lines(x,y)
}

randomisation.phylo<-function(tree, permutations)
{
	nchang<-sum(tree$edge.length)
	permat<-matrix(0,nrow=length(tree$edge.length),ncol=permutations)
	brk<-vector(mode="numeric",length=length(tree$edge.length))
	brk[1]<-1/length(tree$edge.length)
	for (i in 2:length(brk)) brk[i]<-brk[i-1]+1/length(tree$edge.length)
	for (i in 1:permutations) {
		randno<-runif(nchang,min=0,max=1)
		randno<-sort(randno)
		for (j in 1:length(tree$edge.length)) {
			permat[j,i]<-length(grep(TRUE,randno <= brk[j]))
			randno<-randno[-grep(TRUE,randno <= brk[j])]
		}
	}
	probs<-vector(mode="numeric")
	for (i in 1:length(tree$edge.length)) {
		probs[i]<-wilcox.test(tree$edge.length[i],permat[i,],alternative="greater")$p.value
	}
	probs
}

randombranch.phylo<-function(tree, ttree, comp, nchar, permutations)
{
	charcorr<-c((comp[ttree$tip.label,]/nchar),rep(1,Nnode(ttree)))
	ttree$edge.length<-charcorr[tree$edge[,2]]*ttree$edge.length
	nchang<-sum(tree$edge.length)
	permat<-matrix(0,nrow=length(tree$edge.length),ncol=permutations)
	brk<-vector(mode="numeric",length=length(ttree$edge.length))
	brk[1]<-ttree$edge.length[1]/length(ttree$edge.length)
	for (i in 2:length(brk)) brk[i]<-brk[i-1]+ttree$edge.length[i]/length(ttree$edge.length)
	brk<-brk/(sum(ttree$edge.length)/length(ttree$edge.length))
	for (i in 1:permutations) {
		randno<-runif(nchang,min=0,max=1)
		randno<-sort(randno)
		for (j in 1:length(tree$edge.length)) {
			breaker<-length(grep(TRUE,randno <= brk[j]))
			permat[j,i]<-breaker
			if (breaker > 0) randno<-randno[-(1:breaker)]
		}
	}
	sigs<-vector(mode="numeric")
	for (i in 1:length(tree$edge.length)) {
		ifelse(tree$edge.length[i] > sort(permat[i,])[ceiling(0.95*length(permat[i,]))],sigs[i]<-1,sigs[i]<-0)
	}
	sigs
}

get.branches<-function(node, froms, tos)
{
	branches<-vector(mode="numeric")
	for (i in 1:length(froms)) {
		for (j in 1:length(node)) {
			branches[length(branches)+1:length(branches)+2]<-grep(node[j],froms)
			branches<-sort(unique(branches))
		}
		node<-unique(c(node,froms[sort(match(tos[branches],froms))]))
	}
	branches<-vector(mode="numeric")
	for (i in 1:length(node)) {
		branches[length(branches)+1:length(branches)+2]<-grep(node[i],froms)
		branches<-sort(unique(branches))
	}
	branches
}

broken.stick<-function(timeseries,time)
{
	require(stats)
	require(nlme)
	require(paleoTS)
	if(length(grep("1 1",paste(diff(grep(TRUE,diff(timeseries) == 0)),collapse=" "))) > 0) { # Case if 4 values in a row are the same
		print("Some bins deleted due to consecutive equal values (confounding model comparison)")
		deletes<-vector(mode="numeric")
		for(i in 2:length(timeseries)) {
			if(timeseries[i-1] == timeseries[i]) deletes[length(deletes)+1]<-i
		}
		timeseries<-timeseries[-deletes]
		time<-time[-deletes]
	}
	tslength<-length(timeseries)
	onestick<-lm(timeseries~time) # First fit single stick model
	onestickAICc<-AICc(onestick,tslength)
	if(length(timeseries) < 15 && length(timeseries) >= 10) print("Time series too short (too few bins) for three stick model") # Case if broken stick pointless
	if(length(timeseries) < 10) print("Time series too short (too few bins) for two or three stick model") # Case if broken stick pointless
	if(length(timeseries) >= 10) { # Case if two stick doable
		splits<-vector(length=tslength-9,mode="numeric")
		for (i in 1:(tslength-9)) splits[i]<-(4+i) # Find all splits of at least five bins in length
		AICcsplits<-splits # Set up vector to store AICc of two stick models
		for(i in 1:length(splits)) { # Fill the above
			firstsplit<-1:splits[i]
			secondsplit<-(splits[i]+1):tslength
			AICcsplits[i]<-AICc(lm(timeseries[firstsplit]~time[firstsplit]),length(firstsplit))+AICc(lm(timeseries[secondsplit]~time[secondsplit]),length(secondsplit))
		}
		AICcsplits<-as.numeric(AICcsplits)
		twostickAICc<-min(AICcsplits)
		besttwostick<-splits[grep(TRUE,AICcsplits == min(AICcsplits))] # Establish split with lowest combined AICc
		twostickone<-lm(timeseries[1:besttwostick]~time[1:besttwostick]) # Fit first stick of two stick model
		twosticktwo<-lm(timeseries[(besttwostick+1):tslength]~time[(besttwostick+1):tslength]) # Fit second stick of two stick model
	}
	if(length(timeseries) >= 15) { # Case if three stick doable
		splits<-vector(mode="character")
		for(i in 5:(tslength-10)) {
			for(j in (i+5):(tslength-5)) {
				splits[(length(splits)+1)]<-paste(i,":",j,sep="",collapse="")
			}
		}
		AICcsplits<-splits
		for(i in 1:length(splits)) {
			AICcsplits[i]<-AICc(lm(timeseries[1:strsplit(splits[i],":")[[1]][1]]~time[1:strsplit(splits[i],":")[[1]][1]]),length(1:strsplit(splits[i],":")[[1]][1]))+AICc(lm(timeseries[(as.numeric(strsplit(splits[i],":")[[1]][1])+1):strsplit(splits[i],":")[[1]][2]]~time[(as.numeric(strsplit(splits[i],":")[[1]][1])+1):strsplit(splits[i],":")[[1]][2]]),length((as.numeric(strsplit(splits[i],":")[[1]][1])+1):strsplit(splits[i],":")[[1]][2]))+AICc(lm(timeseries[(as.numeric(strsplit(splits[i],":")[[1]][2])+1):tslength]~time[(as.numeric(strsplit(splits[i],":")[[1]][2])+1):tslength]),length((as.numeric(strsplit(splits[i],":")[[1]][2])+1):tslength))
		}
		AICcsplits<-as.numeric(AICcsplits)
		threestickAICc<-min(AICcsplits)
		bestthreestick<-splits[grep(TRUE,AICcsplits == min(AICcsplits))] # Establish split with lowest combined AICc
		threestickone<-lm(timeseries[1:as.numeric(strsplit(bestthreestick,":")[[1]][1])]~time[1:as.numeric(strsplit(bestthreestick,":")[[1]][1])]) # Fit first stick of three stick model
		threesticktwo<-lm(timeseries[(as.numeric(strsplit(bestthreestick,":")[[1]][1])+1):as.numeric(strsplit(bestthreestick,":")[[1]][2])]~time[(as.numeric(strsplit(bestthreestick,":")[[1]][1])+1):as.numeric(strsplit(bestthreestick,":")[[1]][2])]) # Fit second stick of three stick model
		threestickthree<-lm(timeseries[(as.numeric(strsplit(bestthreestick,":")[[1]][2])+1):tslength]~time[(as.numeric(strsplit(bestthreestick,":")[[1]][2])+1):tslength]) # Fit second stick of three stick model
	}
	modelcomp<-cbind(c(onestickAICc,twostickAICc,threestickAICc),round(akaike.wts(c(onestickAICc,twostickAICc,threestickAICc)),3))
	rownames(modelcomp)<-c("One stick","Two stick","Three stick")
	colnames(modelcomp)<-c("AICc","Akaike wt")
	result<-list(timeseries,time,onestick,twostickone,twosticktwo,threestickone,threesticktwo,threestickthree,besttwostick,bestthreestick,modelcomp)
	names(result)<-c("y.values","x.values","onestick.model","twostick.model1","twostick.model2","threestick.model1","threestick.model2","threestick.model3","twostick.split","threestick.split","model.comparison")
	return(result)
}

mars<-function(x,y)
{
	require(earth)
	AICs<-vector(mode="numeric")
	for(i in 3:21) { # Find optimal number of knots using AIC
		marsout<-earth(x,y,nk=i)
		cuts<-sort(x[match(sort(unique(marsout$cuts)),x)])
		rss<-marsout$rss
		AICs[(i-2)]<-(2*(length(cuts)+1))+(length(x)*log(rss))
	}
	nk<-max(grep(TRUE,min(AICs) == AICs))+2
	marsout<-earth(x,y,nk=nk)
	prediction<-y-marsout$residuals
	cuts<-sort(x[match(sort(unique(marsout$cuts)),x)])
	residuals<-marsout$residuals
	rss<-marsout$rss
	slopes<-vector(mode="numeric",length=length(cuts)+1)
	breaks<-sort(unique(c(1,match(cuts,x),length(x))))
	for(i in 2:length(breaks)) slopes[(i-1)]<-prediction[breaks[i]]-prediction[breaks[(i-1)]]
	slopes<-round(slopes,4)
	AIC<-(2*(length(cuts)+1))+(length(x)*log(rss))
	linearRSS<-sum((y-((lsfit(x,y)$coefficients[2]*x)+lsfit(x,y)$coefficients[1]))^2)
	linearAIC<-4+(length(x)*log(linearRSS))
	ifelse(linearAIC <= AIC,best.model<-"Linear",best.model<-paste("MARS with ",length(cuts)," hinges",sep="",collapse=""))
	result<-list(prediction,cuts,residuals,rss,slopes,AIC,linearRSS,linearAIC,best.model)
	names(result)<-c("prediction","cuts","residuals","rss","slopes","AIC","linearRSS","linearAIC","best.model")
	return(result)
}




# sqs version 2.0 by John Alroy
# performs shareholder quorum subsampling on an array of specimen counts
# can be used to perform classical rarefaction instead of SQS
# written 29 July 2010; version 2.0 completed 14 February 2011
# changes in version 2.0: improved subsampling algorithm; including the dominant 
#  taxon is now the default; improved reporting of errors and basic statistics
# warning: do not use this program with taxonomic occurrence data drawn from
#  multiple published references because it is not designed to count
#  single-reference taxa or adjust for long taxonomic lists
# warning: version 1.0 yields estimates that are downwards-biased when q < 0.6
#  and abundance distributions are highly uneven
#
# Modified by GTL 22/03/11 to include single publication occurrence correction
sqs<-function(ab,q,trials,method,dominant,p.1)	{
	
	params <- array(data=NA,dim=0,dimnames=c("raw richness"))
	if (missing(trials))	{
		trials <- 100
	}
	if (missing(method))	{
		method <- ""
	} else if (method != "" && method != "rarefaction" && method != "CR")	{
		return(print('If the method is rarefaction enter method="rarefaction" or "CR"',quote=F))
	}
	if ((q <= 0 || q >= 1) && method != "rarefaction" && method != "CR")	{
		return(print("If the method is SQS the quota must be greater than zero and less than one",quote=F))
	} else if (q < 1 && (method == "rarefaction" || method == "CR"))	{
		return(print("If the method is rarefaction the quota must be an integer",quote=F))
	}
	if (missing(dominant))	{
		dominant <- 0
	} else if (dominant != "" && dominant != "exclude" && dominant != "no")	{
		return(print('To exclude the dominant taxon, enter dominant="exclude" or "no"',quote=F))
	}
	
	# compute basic statistics
	specimens <- sum(ab)
	singletons <- 0
	doubletons <- 0
	highest <- 0
	for (i in 1:length(ab))	{
		if (ab[i] == 1)	{
			singletons <- singletons + 1
		} else if (ab[i] == 2)	{
			doubletons <- doubletons + 1
		}
		if (ab[i] > highest)	{
			highest <- ab[i]
			mostfrequent <- i
		}
	}
	
	u <- 1 - singletons / specimens
	
	# GTL modification starts
	if (missing(p.1) && dominant == "exclude") {
		u <- 1 - singletons/(specimens - highest)
	}
	if (missing(p.1) && dominant == "no") {
		u <- 1 - singletons/(specimens - highest)
	}
	if (!missing(p.1)) {
		u <- (sum(ab) - p.1)/sum(ab)
	}
	# GTL modification ends

	if (u == 0)	{
		return(print("Coverage is zero because all taxa are singletons",quote=F))
	}
	
	# compute raw taxon frequencies (temporarily)
	freq <- ab / specimens
	
	# standard recursive equation for Fishers alpha
	alpha <- 10
	oldalpha <- 0
	while (abs(alpha - oldalpha) > 0.0000001)	{
		oldalpha <- alpha
		alpha <- length(ab) / log(1 + specimens/alpha)
	}

	params["raw richness"] <- length(ab)
	params["Good's u"] <- u
	params["subsampled richness"] <- NA
	params["subsampled u"] <- NA
	params["Chao 1"] <- length(ab) + singletons**2/(2* doubletons)
	params["subsampled Chao 1"] <- NA
	# governing parameter of the geometric series distribution
	params["k"] <- abs(lm(log(sort(freq)) ~ c(1:length(freq)))$coefficients[2])
	params["Fisher's alpha"] <- alpha
	params["Shannon's H"] <- -1 * sum(freq * log(freq))
	params["Hurlbert's PIE"] <- (1 - sum(freq**2)) * length(ab) / (length(ab) - 1)
	params["dominance"] <- highest / specimens
	params["specimens"] <- specimens
	params["singletons"] <- singletons
	params["doubletons"] <- doubletons
	params["specimens drawn"] <- 0

	if (dominant != "exclude" && dominant != "no")	{
		highest <- 0
		mostfrequent <- 0
	}

	# return if the quorum target is higher than overall coverage
	if ((q > u && method != "rarefaction" && method != "CR") || (q >= sum(ab)))	{
		return(params)
	}
	# return if the rarefaction quota is equal to or higher than the
	#  specimen count
	if (method == "rarefaction" && q >= specimens - highest)	{
		return(params)
	}

	# compute adjusted taxon frequencies
	freq <- ab * u / (specimens - highest)

	# create an array in which each cell corresponds to one specimen
	ids <- array()
	n <- 0
	for (i in 1:length(ab))	{
		for (j in 1:ab[i])	{
			n <- n + 1
			ids[n] <- i
		}
	}

	# subsampling trial loop
	# s will be the subsampled taxon count
	s <- array(rep(0,trials))
	subsingle <- array(rep(0,trials))
	subdouble <- array(rep(0,trials))
	subchao <- array(rep(0,trials))
	mostfrequentdrawn <- 0
	for (trial in 1:trials)	{
		pool <- ids
		left <- length(pool)
		seen<- array(data=rep(0,length(ab)))
		subfreq <- array(rep(0,length(ab)))
		if (method != "rarefaction" && method != "CR")	{
			udrawn <- 0
			while (udrawn < q)	{
				# draw a specimen
				x <- floor(runif(1,min=1,max=left+1))
				# add to frequency and taxon sums if species has
				#  not been drawn previously
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if (seen[pool[x]] == 0)	{
					if (pool[x] != mostfrequent)	{
						udrawn <- udrawn + freq[pool[x]]
					}
					seen[pool[x]] <- 1
					# randomly throw back some draws that put the sum over q
					#  (improved algorithm added in version 2.0)
					plus <- 1
					if (udrawn > q)	{
						plus <-  1 - q / udrawn
					}
					if (runif(1) <= plus)	{
						s[trial] <- s[trial] + 1
					} else	{
						subfreq[pool[x]] <- subfreq[pool[x]] - 1
					}
				}
				# decrease pool of specimens not yet drawn
				pool[x] <- pool[left]
				left <- left - 1
			}
		} else	{
			i <- 0
			draws <- 0
			while (i < q)	{
				draws <- draws + 1
				x <- floor(runif(1,min=1,max=length(ids)-draws+2))
				subfreq[pool[x]] <- subfreq[pool[x]] + 1
				if (pool[x] != mostfrequent)	{
					i <- i + 1
				}
				if (seen[pool[x]] == 0)	{
					seen[pool[x]] <- 1
					s[trial] <- s[trial] + 1
				}
				pool[x] <- pool[length(ids)-draws+1]
			}
		}
		for (i in 1:length(ab))	{
			if (subfreq[i] == 1 && i != mostfrequent)	{
				subsingle[trial] <- subsingle[trial] + 1
			} else if (subfreq[i] == 2 && i != mostfrequent)	{
				subdouble[trial] <- subdouble[trial] + 1
			}
		}
		if (subsingle[trial] > 0 && subdouble[trial] > 0)	{
			subchao[trial] <- s[trial] + subsingle[trial]**2/(2*subdouble[trial])
		} else	{
			subchao[trial] <- s[trial]
		}
		params["specimens drawn"] <- params["specimens drawn"] + sum(subfreq)
		if (mostfrequent != 0)	{
			mostfrequentdrawn <- mostfrequentdrawn + subfreq[mostfrequent]
		}
	}
	params["specimens drawn"] <- params["specimens drawn"] / trials
	# compute vector of non-zero counts
	options(warn=-1)
	s2 <- sort(sqrt(s-1))^2+1
	options(warn=0)
	# compute geometric mean
	params["subsampled richness"] <- exp(mean(log(s2))) * length(s2)/length(s)
	# use of arithmetic means to compute Goods u is adequate
	mostfrequentdrawn <- mostfrequentdrawn / trials
	params["subsampled u"] <- 1 - mean(subsingle) / (params["specimens drawn"] - mostfrequentdrawn)
	params["subsampled Chao 1"] <- exp(mean(log(subchao)))
	return(params)
	
}

# sqs version 1.0 by John Alroy
# performs shareholder quorum subsampling on an array of specimen counts
# set method="rarefaction" or "CR" to perform classical rarefaction instead of SQS
# written 29 July 2010
# Modified by GTL 14/12/10 to include single publication occurrence correction
#sqs<-function(ab,q,trials,method,dominant,p.1)
#{
#	params <- array(data=NA,dim=0,dimnames=c("raw richness"))
#	if (missing(trials))  {
#		trials <- 100
#	}
#	if (missing(method))  {
#		method <- ""
#	}
#	if ((q <= 0 || q >= 1) && method != "rarefaction" && method != "CR")  {
#		print("If the method is SQS the quota must be greater than zero and less than one")
#		return(params)
#	} else if (q < 1 && (method == "rarefaction" || method == "CR"))  {
#		print("If the method is rarefaction the quota must be an integer")
#		return(params)
#	}
#	if (missing(dominant))  {
#		dominant <- 0
#	}
#	
#	# compute basic statistics
#	specimens <- sum(ab)
#	singletons <- 0
#	doubletons <- 0
#	highest <- 0
#	for (i in 1:length(ab))  {
#		if (ab[i] == 1)  {
#			singletons <- singletons + 1
#		} else if (ab[i] == 2)  {
#			doubletons <- doubletons + 1
#		}
#		if (ab[i] > highest)  {
#			highest <- ab[i]
#			mostfrequent <- i
#		}
#	}
#	# exclude dominant taxon unless told to include it
#	if (dominant == "include")  {
#		highest <- 0
#	}
#	
#	if (missing(p.1)) {
#		u <- 1 - singletons/(specimens - highest)
#	}
#	if (!missing(p.1)) {
#		u <- (sum(ab) - p.1)/sum(ab)
#	}
#	if (u == 0)  {
#		print("Coverage is zero because all taxa are singletons")
#		return(params)
#	}
#	
#	# compute raw taxon frequencies (temporarily)
#	freq <- ab / specimens
#	
#	# standard recursive equation for Fishers alpha
#	alpha <- 10
#	oldalpha <- 0
#	while (abs(alpha - oldalpha) > 0.0000001)  {
#		oldalpha <- alpha
#		alpha <- length(ab) / log(1 + specimens/alpha)
#	}
#	
#	params["raw richness"] <- length(ab)
#	params["Good's u"] <- u
#	params["subsampled richness"] <- NA
#	params["Chao 1"] <- length(ab) + singletons**2/(2* doubletons)
#	params["subsampled Chao 1"] <- NA
#	# governing parameter of the geometric series distribution
#	params["k"] <- abs(lm(log(sort(freq)) ~ c(1:length(freq)))$coefficients[2])
#	params["Fisher's alpha"] <- alpha
#	params["Shannon's H"] <- -1 * sum(freq * log(freq))
#	params["Hurlbert's PIE"] <- (1 - sum(freq**2)) * length(ab) / (length(ab) - 1)
#	params["dominance"] <- highest / specimens
#	params["specimens"] <- specimens
#	params["singletons"] <- singletons
#	params["doubletons"] <- doubletons
#	params["specimens drawn"] <- 0
#	
#	# return if the quorum target is higher than overall coverage
#	if ((q > u && method != "rarefaction" && method != "CR") || (q >= sum(ab)))  {
#		return(params)
#	}
#	
#	# compute adjusted taxon frequencies
#	freq <- ab * u / (specimens - highest)
#	
#	# create an array in which each cell corresponds to one specimen
#	ids <- array()
#	n <- 0
#	for (i in 1:length(ab))  {
#		for (j in 1:ab[i])  {
#			n <- n + 1
#			ids[n] <- i
#		}
#	}
#	
#	# subsampling trial loop
#	# s will be the subsampled taxon count
#	s <- array(rep(0,trials))
#	subchao <- array(rep(0,trials))
#	for (trial in 1:trials)  {
#		pool <- ids
#		left <- length(pool)
#		seen <- array(data=rep(0,length(ab)))
#		subfreq <- array(rep(0,length(ab)))
#		
#		if (method != "rarefaction" && method != "CR")  {
#			sumfreq <- 0
#			under <- q
#			while (sumfreq < q)  {
#				# draw a specimen
#				x <- floor(runif(1,min=1,max=left+1))
#				subfreq[pool[x]] <- subfreq[pool[x]] + 1
#				# add to frequency and taxon sums if species has
#				#  not been drawn previously
#				if (seen[pool[x]] == 0)  {
#					if (pool[x] != mostfrequent || dominant == "include")  {
#						sumfreq <- sumfreq + freq[pool[x]]
#					}
#					seen[pool[x]] <- 1
#					# count the taxon putting the sum over the target
#					#  only if the resulting overshoot would be less
#					#  than the undershoot created by not counting it
#					if ((sumfreq >= q && sumfreq - q < under) || sumfreq < q)  {
#						s[trial] <- s[trial] + 1
#					} else  {
#						subfreq[pool[x]] <- subfreq[pool[x]] - 1
#					}
#					under <- q - sumfreq
#				}
#				# decrease pool of specimens not yet drawn
#				pool[x] <- pool[left]
#				left <- left - 1
#			}
#		} else  {
#			for (i in 1:q)  {
#				x <- floor(runif(1,min=1,max=length(ids)-i+2))
#				subfreq[pool[x]] <- subfreq[pool[x]] + 1
#				if (seen[pool[x]] == 0)  {
#					seen[pool[x]] <- 1
#					s[trial] <- s[trial] + 1
#				}
#				pool[x] <- pool[length(ids)-i+1]
#			}
#		}
#		subsingle <- 0
#		subdouble <- 0
#		for (i in 1:length(ab))  {
#			if (subfreq[i] == 1)  {
#				subsingle <- subsingle + 1
#			} else if (subfreq[i] == 2)  {
#				subdouble <- subdouble + 1
#			}
#		}
#		if (subsingle > 0 && subdouble > 0)  {
#			subchao[trial] <- s[trial] + subsingle**2/(2*subdouble)
#		} else  {
#			subchao[trial] <- s[trial]
#		}
#		params["specimens drawn"] <- params["specimens drawn"] + sum(subfreq)
#	}
#	
#	params["subsampled richness"] <- exp(mean(log(s)))
#	params["subsampled Chao 1"] <- exp(mean(log(subchao)))
#	params["specimens drawn"] <- params["specimens drawn"] / trials
#	return(params)
#}

gen.diff<-function(x,time)
{
	#if(cor.test(time,x)$p.value > 0.05) print("Warning: variables not significantly correlated, generalised differencing not recommended")
	dt<-x-((lsfit(time,x)$coefficients[2]*time)+lsfit(time,x)$coefficients[1])
	m<-lsfit(dt[1:(length(dt)-1)],dt[2:length(dt)])$coefficients[2]
	gendiffs<-dt[1:(length(dt)-1)]-(dt[2:length(dt)]*m)
	gendiffs
}

getthedamndatain<-function(file, sep="\t", header=F)
{
	X<-scan(file = file, what = "", sep = "\n", quiet = TRUE)
	data<-matrix(nrow=length(X),ncol=length(strsplit(X[1],sep)[[1]]))
	for(i in 1:length(X)) {
		data[i,1:length(strsplit(X[1],sep)[[1]])]<-strsplit(X[i],sep)[[1]]
	}
	if(header == TRUE) {
		colnames(data)<-data[1,]
		data<-data[-1,]
	}
	return(data)
}

my.ccf<-function(x,y,xlab,ylab) {
	completeforboth<-intersect(grep(TRUE,!is.na(x)),grep(TRUE,!is.na(y)))
	differences<-diff(completeforboth)
	if(length(grep(TRUE,differences != 1)) > 0) { # Does gap exist?
		for(i in 1:length(grep(TRUE,differences != 1))) {
			if(i == 1) {
				split<-completeforboth[1:grep(TRUE,differences != 1)[i]]
				if(length(grep(TRUE,differences != 1)) == 1) {
					split.a<-completeforboth[(grep(TRUE,differences != 1)[i]+1):length(completeforboth)]
					if(length(split.a) > length(split)) split<-split.a
				}
			}
			if(i > 1) {
				split.a<-completeforboth[(grep(TRUE,differences != 1)[(i-1)]+1):grep(TRUE,differences != 1)[i]]
				if(length(split.a) > length(split)) split<-split.a
			}
			if(i == length(grep(TRUE,differences != 1)) && i > 1) {
				split.a<-completeforboth[(grep(TRUE,differences != 1)[i]+1):length(completeforboth)]
				if(length(split.a) > length(split)) split<-split.a
			}
		}
	}
	if(length(grep(TRUE,differences != 1)) == 0) split<-completeforboth
	x<-x[split]
	y<-y[split]
	plot.title<-paste(xlab,"vs.",ylab,"\n(Continuous from bin",split[1],"to",split[length(split)],")")
	ccf(y,x,main=plot.title)
}

find.ancestor<-function(descs,tree) {
    require(ape)
    tipnos<-match(descs,tree$tip.label)
    anc.node<-sort(unique(tree$edge[,1][match(tipnos,tree$edge[,2])]))
    while(length(anc.node) > 1) {
        highestnode<-anc.node[length(anc.node)]
        anc.node<-anc.node[-length(anc.node)]
        anc.node<-sort(unique(c(anc.node,tree$edge[match(highestnode,tree$edge[,2]),1])))
    }
    return(anc.node)
}

get.nexus<-function(file)
{
	# Read in raw file:
	X<-scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
	# Remove trees:
	if(length(grep("BEGIN TREES",X)) > 0) { # If trees are present
		deleteline<-grep("BEGIN TREES",X)-1 # Find start line
		X<-X[1:deleteline] # And remove
	}
    if(length(grep("\\[!",X)) > 0) {
        textlines<-grep("\\[!",X)
        textlines<-gsub("]","",gsub("\\[!","",X[textlines]))
    } else {
        textlines<-""
    }
	dimensionsline<-grep("DIMENSIONS",X)
	ntax<-as.numeric(gsub("NTAX=","",strsplit(X[dimensionsline]," ")[[1]][grep("NTAX",strsplit(X[dimensionsline]," ")[[1]],fixed=TRUE)]))
	nchar<-as.numeric(gsub(";","",gsub(" ","",gsub("NCHAR=","",strsplit(X[dimensionsline]," ")[[1]][grep("NCHAR",strsplit(X[dimensionsline]," ")[[1]],fixed=TRUE)]))))
	formatline<-grep("FORMAT",X)
	missing<-gsub("MISSING=","",strsplit(X[formatline]," ")[[1]][grep("MISSING=",strsplit(X[formatline]," ")[[1]],fixed=TRUE)])
	gap<-gsub("GAP=","",strsplit(X[formatline]," ")[[1]][grep("GAP=",strsplit(X[formatline]," ")[[1]],fixed=TRUE)])
	matrixline<-grep("MATRIX",X,fixed=TRUE)
	matrixblock<-X[(matrixline+1):(matrixline+ntax)]
	# Remove quotation marks and tabs from matrixblock:
	for(i in 1:length(matrixblock)) matrixblock[i]<-gsub("\t"," ",gsub("\x94","",gsub("\x93","",gsub("\xd5","",gsub("\xd4","",gsub("\xd3","",gsub("\xd2","",gsub("'","",gsub("\"","",matrixblock[i])))))))))
	MATRIX<-matrix(nrow=ntax,ncol=nchar)
	MATRIXrn<-vector(mode="character")
	for(i in 1:length(matrixblock)) {
		nameplusspaces<-paste(strsplit(matrixblock[i],"")[[1]][1:max(grep(TRUE,strsplit(matrixblock[i],"")[[1]] == " "))],collapse="")
		reversed.elements<-rev(strsplit(nameplusspaces,"")[[1]])
		just.name<-reversed.elements[min(grep(TRUE,reversed.elements != " ")):length(reversed.elements)]
		MATRIXrn[i]<-paste(rev(just.name),collapse="")
	}
	rownames(MATRIX)<-MATRIXrn
	for(i in 1:ntax) {
		matrixrowdata<-rev(rev(strsplit(matrixblock[i],"")[[1]])[1:(min(grep(TRUE,rev(strsplit(matrixblock[i],"")[[1]]) == " "))-1)])
		j<-1 # Set place in character sequence
		nc<-1 # Set actual character
		while(j <= length(matrixrowdata)) {
			# Normal character (copy as is):
			if(matrixrowdata[j] != "(" && matrixrowdata[j] != "{" && matrixrowdata[j] != "[" && matrixrowdata[j] != missing && matrixrowdata[j] != gap) MATRIX[i,nc]<-matrixrowdata[j]
			# Missing/gap character (replace with NA):
			if(matrixrowdata[j] == missing || matrixrowdata[j] == gap) MATRIX[i,nc]<-NA
			# Polymorphic character (separate states with ampersands):
			if(matrixrowdata[j] == "(" || matrixrowdata[j] == "{" || matrixrowdata[j] == "[") {
				start<-j
				while(matrixrowdata[j] != ")" && matrixrowdata[j] != "}" && matrixrowdata[j] != "]") j<-j+1
				end<-j
				MATRIX[i,nc]<-paste(matrixrowdata[(start+1):(end-1)],collapse="&")
			}
			j<-j+1 # Text character count
			nc<-nc+1 # Actual character count
		}
	}
    # Convert alphabetic characters to numerical values:
    MATRIX[grep(TRUE,MATRIX == "A")] <- 10
    MATRIX[grep(TRUE,MATRIX == "B")] <- 11
    MATRIX[grep(TRUE,MATRIX == "C")] <- 12
    MATRIX[grep(TRUE,MATRIX == "D")] <- 13
    MATRIX[grep(TRUE,MATRIX == "E")] <- 14
    MATRIX[grep(TRUE,MATRIX == "F")] <- 15
    MATRIX[grep(TRUE,MATRIX == "G")] <- 16
    MATRIX[grep(TRUE,MATRIX == "H")] <- 17
    MATRIX[grep(TRUE,MATRIX == "I")] <- 18
    MATRIX[grep(TRUE,MATRIX == "J")] <- 19
    MATRIX[grep(TRUE,MATRIX == "K")] <- 20
    MATRIX[grep(TRUE,MATRIX == "L")] <- 21
    MATRIX[grep(TRUE,MATRIX == "M")] <- 22
    MATRIX[grep(TRUE,MATRIX == "N")] <- 23
    MATRIX[grep(TRUE,MATRIX == "O")] <- 24
    MATRIX[grep(TRUE,MATRIX == "P")] <- 25
    MATRIX[grep(TRUE,MATRIX == "Q")] <- 26
    MATRIX[grep(TRUE,MATRIX == "R")] <- 27
    MATRIX[grep(TRUE,MATRIX == "S")] <- 28
    MATRIX[grep(TRUE,MATRIX == "T")] <- 29
    MATRIX[grep(TRUE,MATRIX == "U")] <- 30
    MATRIX[grep(TRUE,MATRIX == "V")] <- 31
    MATRIX[grep(TRUE,MATRIX == "W")] <- 32
    MATRIX[grep(TRUE,MATRIX == "X")] <- 33
    MATRIX[grep(TRUE,MATRIX == "Y")] <- 34
    MATRIX[grep(TRUE,MATRIX == "Z")] <- 35
	if (length(grep("TYPESET",X)) > 0) { # If ordering is specified
		orderingline<-grep("TYPESET",X)
		orderingline<-gsub(",","",gsub(";","",gsub("\t","",gsub("TYPESET","",gsub("\\*","",gsub("UNTITLED","",gsub("=","",X[orderingline],fixed=TRUE))))))) # Clean up
		unordered<-strsplit(gsub("unord:","",strsplit(orderingline," ord:")[[1]][1])," ")[[1]] # Get unordered list
		ordered<-strsplit(strsplit(orderingline," ord:")[[1]][2]," ")[[1]] # Get ordered list
		unordered<-unordered[-grep(TRUE,unordered == "")] # Remove empty cells
		ordered<-ordered[-grep(TRUE,ordered == "")] # Remove empty cells
		# Break out ranges into individual characters:
		if(length(grep("-",unordered)) > 0) {
			while(length(grep("-",unordered)) > 0) {
				unordered<-c(unordered,as.numeric(strsplit(unordered[grep("-",unordered)[1]],"-")[[1]][1]):as.numeric(strsplit(unordered[grep("-",unordered)[1]],"-")[[1]][2]))
				unordered<-unordered[-grep("-",unordered)[1]]
			}
		}
		if(length(grep("-",ordered)) > 0) {
			while(length(grep("-",ordered)) > 0) {
				ordered<-c(ordered,as.numeric(strsplit(ordered[grep("-",ordered)[1]],"-")[[1]][1]):as.numeric(strsplit(ordered[grep("-",ordered)[1]],"-")[[1]][2]))
				ordered<-ordered[-grep("-",ordered)[1]]
			}
		}
		# Place in numerical order:
		unordered<-sort(as.numeric(unordered))
		ordered<-sort(as.numeric(ordered))
		ordering<-vector(mode="character")
		ordering[unordered]<-"UNORD"
		ordering[ordered]<-"ORD"		
	} else {
		orderingline<-grep("DEFTYPE",X)
		if(strsplit(strsplit(X[orderingline],"DEFTYPE=")[[1]][2]," ")[[1]][1] == "unord") ordering<-rep("UNORD",nchar) # Case if all unordered
		if(strsplit(strsplit(X[orderingline],"DEFTYPE=")[[1]][2]," ")[[1]][1] == "ord") ordering<-rep("ORD",nchar) # Case if all ordered
	}
    if (length(grep("WTSET",X)) > 0) { # If weights are specified
		weightsline<-grep("WTSET",X)
		weightsline<-gsub(",","",gsub(";","",gsub("\t","",gsub("WTSET","",gsub("\\*","",gsub("Weights","",gsub("UNTITLED","",gsub("=","",X[weightsline],fixed=TRUE)))))))) # Clean up
        weightsline<-strsplit(weightsline,":")[[1]]
        weightvals<-as.numeric(weightsline[1]) # Weight of first weight (will add other weights to this)
        totalweightno<-length(weightsline)-1 # Number of different weights
        for (i in 2:totalweightno) weightvals<-c(weightvals,as.numeric(strsplit(weightsline[i]," ")[[1]][length(strsplit(weightsline[i]," ")[[1]])]))
        weights<-vector(mode="numeric")
        for (i in 2:length(weightsline)) {
            charnos<-strsplit(weightsline[i]," ")[[1]][grep(FALSE,strsplit(weightsline[i]," ")[[1]] == "")] # Get character numbers
            if(i < length(weightsline)) charnos<-charnos[-length(charnos)]
            if(length(grep("-",charnos)) > 0) {
                blocks<-charnos[grep("-",charnos)]
                charnos<-charnos[-grep("-",charnos)]
                for(j in blocks) charnos<-c(as.numeric(charnos),as.numeric(strsplit(j,"-")[[1]][1]):as.numeric(strsplit(j,"-")[[1]][2]))
            }
            weights[sort(as.numeric(charnos))]<-weightvals[(i-1)]
        }	
	} else {
		weights<-rep(1,nchar) # Case if no weighting (i.e. all equal weight)
	}
	min.vals<-max.vals<-vector(length=nchar)
	for(i in 1:nchar) {
		non.nas<-sort(MATRIX[,i])
		# Case if polymorphisms are available:
		while(length(grep("&",non.nas)) > 0) {
			polymorphism<-grep("&",non.nas)
			non.nas<-c(non.nas,strsplit(non.nas[polymorphism],"&")[[1]])
			non.nas<-non.nas[-polymorphism]
		}
		max.vals[i]<-max(non.nas)
		min.vals[i]<-min(non.nas)
	}
	# Remove characters not favoured by phylogenetic software (spaces, slashes, dashes and pluses):
	rownames(MATRIX)<-gsub("-","_dash_",gsub("/","_slash_",gsub("\\+","_plus_",gsub(" ","_",rownames(MATRIX)))))
	result<-list(textlines,MATRIX,ordering,weights,max.vals,min.vals)
	names(result)<-c("header","matrix","ordering","weights","max.vals","min.vals")
	return(result)
    # To do:
    # - Spit out cleaned up non-_ OTUs
}

# Function:
dist.clad.matrix<-function(clad.matrix, dist.method="euclidean") {
    ordering<-clad.matrix$ordering
    max.vals<-clad.matrix$max.vals
    min.vals<-clad.matrix$min.vals
    weights<-clad.matrix$weights
    clad.matrix<-clad.matrix$matrix
    # Distance matrices for storing:
    comp.char.matrix<-max.dist.matrix<-dist.matrix<-matrix(0, nrow=length(rownames(clad.matrix)), ncol=length(rownames(clad.matrix)))
    # Fill comparable characters diagonal:
    for(i in 1:length(clad.matrix[,1])) comp.char.matrix[i,i]<-length(clad.matrix[i,])-length(grep(TRUE,is.na(clad.matrix[i,])))
    # Main loop (finds distances):
    for(i in 1:(length(clad.matrix[,1])-1)) {
        for(j in (i+1):length(clad.matrix[,1])) {
            # Comparable characters:
            compchar<-intersect(grep(TRUE,!is.na(clad.matrix[rownames(clad.matrix)[i],])),grep(TRUE,!is.na(clad.matrix[rownames(clad.matrix)[j],])))
            # Get sets of comparable characters:
            firstrow<-clad.matrix[rownames(clad.matrix)[i],compchar]
            secondrow<-clad.matrix[rownames(clad.matrix)[j],compchar]
            # Deal with polymorphic characters (if present):
            if(length(grep("&",unique(c(firstrow,secondrow)))) > 0) {
                ampersand.elements<-sort(c(grep("&",firstrow),grep("&",secondrow)))
                for(k in 1:length(ampersand.elements)) {
                    intersection.value<-intersect(strsplit(firstrow[ampersand.elements[k]],"&")[[1]],strsplit(secondrow[ampersand.elements[k]],"&")[[1]])
                    if(length(intersection.value) > 0) { # Case if polymorphic and non-/polymorphic values overlap
                        firstrow[ampersand.elements[k]]<-0
                        secondrow[ampersand.elements[k]]<-0
                    }
                    if(length(intersection.value) == 0) { # Case if polymorphic and non-/polymorphic values do not overlap
                        if(ordering[compchar[ampersand.elements[k]]] == "UNORD") { # Case if character is unordered (max difference is 1)
                            firstrow[ampersand.elements[k]]<-0
                            secondrow[ampersand.elements[k]]<-1
                        }
                        if(ordering[compchar[ampersand.elements[k]]] == "ORD") { # Case if character is ordered (max difference is > 1)
                            firstrowvals<-as.numeric(strsplit(firstrow[ampersand.elements[k]],"&")[[1]]) # Get first row value(s)
                            secondrowvals<-as.numeric(strsplit(secondrow[ampersand.elements[k]],"&")[[1]]) # Get second row value(s)
                            poly.dist.mat<-matrix(0,nrow=length(firstrowvals),ncol=length(secondrowvals)) # Make mini distance matrix
                            for(l in 1:length(firstrowvals)) {
                                for(m in 1:length(secondrowvals)) poly.dist.mat[l,m]<-sqrt((firstrowvals[l]-secondrowvals[m])^2)
                            }
                            firstrow[ampersand.elements[k]]<-0
                            secondrow[ampersand.elements[k]]<-min(poly.dist.mat)
                        }
                    }
                }
            }
            # Get the difference between the two rows:
            diffs<-as.numeric(firstrow) - as.numeric(secondrow)
            # Get absolute differences:
            diffs<-sqrt(diffs*diffs)
            # If there are differences greater than 1 for unordered characters then rescore as 1:
            if(length(grep(TRUE,diffs > 1)) > 0) diffs[grep(TRUE,diffs > 1)[grep(TRUE,ordering[compchar[grep(TRUE,diffs > 1)]] == "UNORD")]]<-1
            # Get weighted differences:
            diffs<-as.numeric(weights[compchar])*diffs
            # Get actual distance:
            actual.dist<-dist(rbind(diffs,rep(0,length(diffs))),method=dist.method)
            # Work out maximum and minimum difference (again, checked against ordering) using compchar characters only:
            maxdiffs<-as.numeric(max.vals[compchar])-as.numeric(min.vals[compchar])
            if(length(grep(TRUE,maxdiffs > 1)) > 0) maxdiffs[grep(TRUE,maxdiffs > 1)[grep(TRUE,ordering[compchar[grep(TRUE,maxdiffs > 1)]] == "UNORD")]]<-1
            maxdiffs<-as.numeric(weights[compchar])*maxdiffs
            max.dist<-dist(rbind(maxdiffs,rep(0,length(maxdiffs))),method=dist.method)
            # Store data:
            dist.matrix[i,j]<-dist.matrix[j,i]<-actual.dist
            max.dist.matrix[i,j]<-max.dist.matrix[j,i]<-max.dist
            comp.char.matrix[i,j]<-comp.char.matrix[j,i]<-length(compchar)
        }
    }
    rownames(comp.char.matrix)<-colnames(comp.char.matrix)<-rownames(max.dist.matrix)<-colnames(max.dist.matrix)<-rownames(dist.matrix)<-colnames(dist.matrix)<-rownames(clad.matrix)
    result<-list(dist.matrix, max.dist.matrix, comp.char.matrix)
    names(result)<-c("dist.matrix","max.dist.matrix","comp.char.matrix")
    return(result)
}

get.anc.states<-function(nexus.matrix, tree) {
    require(ape)
    anc.lik.matrix<-matrix(nrow=Nnode(tree),ncol=length(nexus.matrix$matrix[1,])) # Create matrix to record ancestral state estimates
    rownames(anc.lik.matrix)<-c((Ntip(tree)+1):(Ntip(tree)+Nnode(tree))) # Label matrix to record ancestral state estimates
    for(i in 1:length(nexus.matrix$matrix[1,])) { # Cycle through characters:
        tipstogo<-rownames(nexus.matrix$matrix)[grep(TRUE,is.na(nexus.matrix$matrix[,i]))] # Find taxa with missing data:
        if(length(tipstogo) > 0) chartree<-drop.tip(tree,tipstogo) # Remove tips with missing data
        if(length(tipstogo) == 0) chartree<-tree # If no missing data use whole tree    
        tipvals<-nexus.matrix$matrix[chartree$tip.label,i]
        # Create state change matrix (important for ordered/unordered characters:
        minval<-as.numeric(nexus.matrix$min.vals[i])
        maxval<-as.numeric(nexus.matrix$max.vals[i])
        if(maxval-minval == 1) mymodel<-matrix(c(minval,maxval,maxval,minval),2) # Case for binary character (ordering irrelevant)
        if(maxval-minval > 1 && nexus.matrix$ordering[i] == "UNORD") { # Case for unordered multistate character
            mymodel<-matrix(1,ncol=(maxval-minval)+1,nrow=(maxval-minval)+1)
            for(j in 1:length(mymodel[1,])) mymodel[j,j]<-0
        }
        if(maxval-minval > 1 && nexus.matrix$ordering[i] == "ORD") { # Case for ordered multistate character
            mymodel<-matrix(0,ncol=(maxval-minval)+1,nrow=(maxval-minval)+1)
            for(j in 1:length(mymodel[1,])) {
                for(k in 1:length(mymodel[1,])) {
                    mymodel[j,k]<-sqrt((diff(c(j,k)))^2)
                }
            }
        }
        # Ascertain if polymorphisms are present and if so compute all possible prmutations
        if(length(grep("&",tipvals)) > 0) {
            permutation.counts<-vector(mode="numeric") # Vector to store permutation numbers for each taxon (i.e. number of possible states)
            for(j in 1:length(tipvals)) permutation.counts[j]<-length(strsplit(as.character(tipvals[j]),"&")[[1]]) # Get number of permutations for each taxon
            permutation.count<-prod(permutation.counts) # Get product of permutation counts (i.e. total number of permutations)
            permutations<-matrix(nrow=length(tipvals),ncol=permutation.count) # Make permutations matrix
            permutations[grep(TRUE,permutation.counts == 1),1:permutation.count]<-as.numeric(tipvals[grep(TRUE,permutation.counts == 1)]) # Fill in data for non-polymorphic characters
            phase.count<-1 # Permutation phase
            for(j in grep("&",tipvals)) {
                perm<-as.numeric(strsplit(as.character(tipvals[j]),"&")[[1]])
                if(phase.count == 1) permutations[j,]<-perm
                if(phase.count > 1) {
                    perm<-sort(rep(perm,phase.count))
                    permutations[j,]<-perm
                }
                phase.count<-length(perm)
            }
            rownames(permutations)<-names(tipvals) # Make sure tip value names are carried over
            # Ancestor state estimate for first permutation:
            if(length(unique(sort(permutations[,1]))) > 1) { # If there is more than one state
                if(length(mymodel[,1]) > length(sort(unique(as.numeric(permutations[,1]))))) { # Case if model has larger spread than actual tipvalues (e.g. spans 0-5 when 2 is not recorded)
                    anc.lik<-ace(permutations[,1],chartree,type="discrete",model=mymodel[sort(unique(as.numeric(permutations[,1])))+1,sort(unique(as.numeric(permutations[,1])))+1])$lik.anc
                }
                if(length(mymodel[,1]) <= length(sort(unique(as.numeric(permutations[,1]))))) { # Case if model has same spread as tip values (theoretically this should always be true, but in reality not so much)
                    anc.lik<-ace(permutations[,1],chartree,type="discrete",model=mymodel[])$lik.anc
                }
            }
            if(length(unique(sort(permutations[,1]))) == 1) {
                anc.lik<-matrix(0,nrow=length(permutations[,1])-1,ncol=max(c((unique(sort(as.numeric(permutations[,1])))+1),2)))
                colnames(anc.lik)<-c(0:max(c(unique(sort(as.numeric(permutations[,1]))),1)))
                anc.lik[,as.character(unique(sort(as.numeric(permutations[,1]))))]<-rep(1,length(permutations[,1])-1)
            }
            for(j in 2:length(permutations[1,])) { # Cycle through remaining permutations
                if(length(unique(sort(permutations[,j]))) > 1) {
                    if(length(mymodel[,1]) > length(sort(unique(as.numeric(permutations[,j]))))) { # Case if model has larger spread than actual tipvalues (e.g. spans 0-5 when 2 is not recorded)
                        anc.lik<-((anc.lik*(j-1))+ace(permutations[,j],chartree,type="discrete",model=mymodel[sort(unique(as.numeric(permutations[,j])))+1,sort(unique(as.numeric(permutations[,j])))+1])$lik.anc)/j # Get mean ancestral estimations
                    }
                    if(length(mymodel[,1]) <= length(sort(unique(as.numeric(permutations[,j]))))) { # Case if model has same spread as tip values (theoretically this should always be true, but in reality not so much)
                        anc.lik<-((anc.lik*(j-1))+ace(permutations[,j],chartree,type="discrete",model=mymodel)$lik.anc)/j # Get mean ancestral estimations
                    }
                }
                if(length(unique(sort(permutations[,j]))) == 1) {
                    anc.lik.perm<-matrix(0,nrow=length(permutations[,j])-1,ncol=max(c((unique(sort(as.numeric(permutations[,j])))+1),2)))
                    colnames(anc.lik.perm)<-c(0:max(c(unique(sort(as.numeric(permutations[,j]))),1)))
                    anc.lik.perm[,as.character(unique(sort(as.numeric(permutations[,j]))))]<-rep(1,length(permutations[,j])-1)
                    anc.lik<-((anc.lik*(j-1))+anc.lik.perm)/j # Get mean ancestral estimations
                }
            }
        }
        # If no polymorphisms present just get ancestral likelihoods:
        if(length(grep("&",tipvals)) == 0 && length(unique(sort(as.numeric(tipvals)))) > 1) {
            if(length(mymodel[,1]) > length(sort(unique(as.numeric(tipvals))))) { # Case if model has larger spread than actual tipvalues (e.g. spans 0-5 when 2 is not recorded)
                anc.lik<-ace(as.numeric(tipvals),chartree,type="discrete",model=mymodel[sort(unique(as.numeric(tipvals)))+1,sort(unique(as.numeric(tipvals)))+1])$lik.anc # Get ancestral state likelihoods for internal nodes
            }
            if(length(mymodel[,1]) <= length(sort(unique(as.numeric(tipvals))))) { # Case if model has same spread as tip values (theoretically this should always be true, but in reality not so much)
                anc.lik<-ace(as.numeric(tipvals),chartree,type="discrete",model=mymodel)$lik.anc # Get ancestral state likelihoods for internal nodes
            }
        }
        if(length(grep("&",tipvals)) == 0 && length(unique(sort(as.numeric(tipvals)))) == 1) { # Case if all states are the same
            anc.lik<-matrix(0,nrow=length(tipvals)-1,ncol=max(c((unique(sort(as.numeric(tipvals)))+1),2)))
            colnames(anc.lik)<-c(0:max(c(unique(sort(as.numeric(tipvals))),1)))
            anc.lik[,as.character(unique(sort(as.numeric(tipvals))))]<-rep(1,length(tipvals)-1)
        }
        # Convert to just most likely states (i.e. 0 OR 1 for binary):
        anc.lik.vector<-vector(mode="numeric")
        max.lik<-apply(anc.lik,1,max)
        for(j in 1:length(anc.lik[,1])) {
            anc.lik.vector[j]<-as.numeric(colnames(anc.lik)[match(max.lik[j],anc.lik[j,])]) # Case if single most likely state
            if(length(grep(max.lik[j],anc.lik[j,])) > 1) anc.lik.vector[j]<-NA # Case if multiple equally most likely states (i.e. replace with NA)
        }
        chartreenodes<-(Ntip(chartree)+1):(Ntip(chartree)+Nnode(chartree)) # Nodes for character tree
        names(anc.lik.vector)<-chartreenodes
        # Copy information across to ancestral state matrix for whole tree:
        for(j in chartreenodes) {
            descs<-sort(chartree$tip.label[find.descendants(j,chartree)]) # Find character tree descendants
            anc.lik.matrix[as.character(find.ancestor(descs,tree)),i]<-anc.lik.vector[as.character(j)]
        }
    }
    return(anc.lik.matrix)
}

get.bl<-function(state.matrix, tree, weights) {
    require(ape)
    tree.pd<-tree.nc<-tree.cc<-tree # Trees for patristic distance, number of changes and comparable characters
    for(i in 1:length(tree$edge[,1])) {
        node.1<-tree$edge[i,1] # Get nodes for branch
        node.2<-tree$edge[i,2] # Get nodes for branch
        compchar<-sort(intersect(grep(TRUE,!is.na(state.matrix[node.1,])),grep(TRUE,!is.na(state.matrix[node.2,])))) # Characters that can be compared
        if(length(compchar) > 0) {
            polychar<-sort(unique(c(grep("&",state.matrix[node.1,compchar]),grep("&",state.matrix[node.2,compchar])))) # Characters with polymorphisms
            if(length(polychar) == 0) {
                tree.nc$edge.length[i]<-sum(sqrt(((as.numeric(state.matrix[node.1,compchar])-as.numeric(state.matrix[node.2,compchar]))*weights[compchar])^2))
                tree.pd$edge.length[i]<-sum(sqrt(((as.numeric(state.matrix[node.1,compchar])-as.numeric(state.matrix[node.2,compchar]))*weights[compchar])^2))/(length(compchar)*weights[compchar])
            }
            if(length(polychar) > 0) {
                diffs<-sqrt(( as.numeric(state.matrix[node.1,compchar])-as.numeric(state.matrix[node.2,compchar])) ^2)
                for (j in 1:length(polychar)) {
                    # Get minimum difference between polymorphic characters:
                    mindiff<-min(sqrt((as.numeric(strsplit(as.character(state.matrix[node.1,compchar[polychar[j]]]),"&")[[1]])-as.numeric(strsplit(as.character(state.matrix[node.2,compchar[polychar[j]]]),"&")[[1]]))^2))
                    diffs[polychar[j]]<-mindiff # Add to differences vector
                }
                tree.nc$edge.length[i]<-sum(diffs*weights[compchar])
                tree.pd$edge.length[i]<-sum(diffs*weights[compchar])/(length(compchar)*weights[compchar])
            }
            tree.cc$edge.length[i]<-length(compchar)*weights[compchar]
        }
        if(length(compchar) == 0) {
            tree.nc$edge.length[i]<-0
            tree.pd$edge.length[i]<-0
            tree.cc$edge.length[i]<-0
        }
    }
    result<-list(tree.cc,tree.nc,tree.pd)
	names(result)<-c("tree.cc","tree.nc","tree.pd")
	return(result)
}

rbranch.phylo<-function(state.matrix, tree, ttree, cctree, permutations)
{
    charcorr<-(length(state.matrix[1,])-apply(is.na(state.matrix),1,sum))/length(state.matrix[1,]) # Get proportion of completeness of nodes/tips
    ttree$edge.length<-(cctree$edge.length/length(state.matrix[1,]))*ttree$edge.length # Normalise time tree branch lengths to reflect comparable characters number
    nchang<-sum(tree$edge.length) # Total number of character changes on tree
    permat<-matrix(0,nrow=length(tree$edge.length),ncol=permutations) # Matrix to store premutation results
    # Apportion branches to uniform distribution (between 0 and 1):
    brk<-vector(mode="numeric",length=length(ttree$edge.length))
	brk[1]<-ttree$edge.length[1]/length(ttree$edge.length)
    for (i in 2:length(brk)) brk[i]<-brk[i-1]+ttree$edge.length[i]/length(ttree$edge.length)
	brk<-brk/(sum(ttree$edge.length)/length(ttree$edge.length))
    # Assign changes to branches using a uniform distribution:
    for (i in 1:permutations) {
		randno<-runif(nchang,min=0,max=1)
		randno<-sort(randno)
		for (j in 1:length(tree$edge.length)) {
			breaker<-length(grep(TRUE,randno <= brk[j]))
			permat[j,i]<-breaker
			if (breaker > 0) randno<-randno[-(1:breaker)]
		}
	}
    sigs<-vector(mode="numeric") # Vector to store significant results
    for (i in 1:length(tree$edge.length)) ifelse(tree$edge.length[i] > sort(permat[i,])[ceiling(0.95*length(permat[i,]))],sigs[i]<-1,sigs[i]<-0) # Significance test for high rates
	return(sigs)
}

dec.latlong<-function(x) {
	degrees<-as.numeric(strsplit(x,"°")[[1]][1])
	minutes<-as.numeric(strsplit(strsplit(x,"°")[[1]][2],"′")[[1]][1])/60
	seconds<-as.numeric(strsplit(strsplit(strsplit(x,"°")[[1]][2],"′")[[1]][2],"\\.")[[1]][1])/6000
	decimals<-as.numeric(strsplit(strsplit(strsplit(strsplit(x,"°")[[1]][2],"′")[[1]][2],"\\.")[[1]][2],"″")[[1]][1])/1000000
	direction<-strsplit(strsplit(strsplit(strsplit(x,"°")[[1]][2],"′")[[1]][2],"\\.")[[1]][2],"″")[[1]][2]
	if(direction == "E" || direction == "N") out<-degrees+minutes+seconds+decimals
	if(direction == "W" || direction == "S") out<-(degrees+minutes+seconds+decimals)*-1
	return(out)
}

mrp.tree<-function(tree) {
    mrp<-matrix(0,ncol=(Nnode(tree)-1),nrow=Ntip(tree)) # MRP matrix
    rownames(mrp)<-tree$tip.label # Name rows as taxa
    colnames(mrp)<-as.character(c((Ntip(tree)+2):(Ntip(tree)+Nnode(tree)))) # Name columns as nodes
    for(j in c((Ntip(tree)+2):(Ntip(tree)+Nnode(tree)))) { # For each internal node...
        mrp[tree$tip.label[find.descendants(j, tree)],as.character(j)]<-1 # Fill MRP matrix
    }
    mrp<-mrp[sort(rownames(mrp)),]
    return(mrp)
}

find.descendant.edges<-function(n,tree)
{
    n.terminals<-length(find.descendants(n, tree)) # Find number of terminals (i.e. stopping point)
    nodes<-n # Store internal nodes
    edges<-grep(n,tree$edge[,1]) # Store edges
    while(length(grep(TRUE,tree$edge[edges,2] <= Ntip(tree))) < n.terminals) {# Keep going until all descendnat edges are found
        nodes<-tree$edge[edges,2][grep(TRUE,tree$edge[edges,2] > Ntip(tree))] # Get internal nodes found so far
        for(i in nodes) {
            edges<-sort(unique(c(edges,grep(TRUE,tree$edge[,1] == i))))
        }
    }
    return(edges)
}

mrp.trees<-function(trees) {
    ntrees<-length(summary(trees)[,1]) # How many trees are there?
    strict<-consensus(trees) # Make strict consensus tree
    strict<-root(strict,outgroup=strict$tip.label[1],resolve.root=TRUE) # Root it
    n.furcations<-rle(sort(strict$edge[,1]))$lengths # Number of branches emanating from each node
    names(n.furcations)<-rle(sort(strict$edge[,1]))$values # Add node names
    bi.nodes<-as.numeric(names(n.furcations)[grep(TRUE,n.furcations == 2)]) # Get bifurcating nodes
    poly.nodes<-as.numeric(names(n.furcations)[grep(TRUE,n.furcations > 2)]) # Get polytomous nodes
    if(length(grep(TRUE,(Ntip(strict)+1) == bi.nodes)) == 1) bi.nodes<-bi.nodes[-grep(TRUE,(Ntip(strict)+1) == bi.nodes)] # Remove root node if found
    if(length(grep(TRUE,(Ntip(strict)+1) == poly.nodes)) == 1) poly.nodes<-poly.nodes[-grep(TRUE,(Ntip(strict)+1) == bi.nodes)] # Remove root node if found
    all.nodes<-sort(c(bi.nodes,poly.nodes)) # Combination of above
    mrp.all<-matrix(0,nrow=Ntip(strict),ncol=length(all.nodes)) # Make global MRP matrix (initially just for bifurcating nodes)
    rownames(mrp.all)<-sort(strict$tip.label) # Name rows according to taxa
    mrp.weights<-rep(1, length(all.nodes)) # Make global MRP weights (initially just for bifurcating nodes)
    for(i in 1:length(all.nodes)) { # For each bifurcating node in strict consensus
        mrp.all[strict$tip.label[find.descendants(all.nodes[i], strict)],i]<-1 # Fill intial MRP matrix
    }
    tips.to.drop<-vector(mode="character") # Tips not relevant to polytomies
    for(i in 1:Ntip(strict)) { # For each tip
        links<-strict$edge[grep(TRUE,strict$edge[,2] == i),1] # Start to record nodes along path to tip
        if(links != (Ntip(strict)+1)) { # If not already at root
            while(links[length(links)] != (Ntip(strict)+1)) { # keep going until you are
                links<-c(links,strict$edge[grep(TRUE,strict$edge[,2] == links[length(links)]),1]) # Record nodes along path
            }
        }
        if(max(n.furcations[as.character(links)]) == 2) tips.to.drop<-c(tips.to.drop,strict$tip.label[i]) # If a droppable tip add to the list
    }
    if(length(tips.to.drop) > 0) { # Only if there are tips to drop
        pruned.trees<-rmtree(ntrees,Ntip(drop.tip(trees[[1]], tips.to.drop))) # Create new MPT list
        for(i in 1:ntrees) pruned.trees[[i]]<-drop.tip(trees[[i]], tips.to.drop) # without taxa irrelevant to polytomies
        strict<-drop.tip(strict, tips.to.drop) # Remove from consensus too
    }
    nodestogo<-(Ntip(strict)+1):(Ntip(strict)+Nnode(strict)) # Nodes in strict consensus to be ignored in MPTs
    for(i in nodestogo) assign(paste("desc.",i,sep=""),strict$tip.label[find.descendants(i,strict)])
    desc.groups<-vector(mode="character")
    for(i in 1:ntrees) { # For each MPT
        MPT.nodes<-(Ntip(pruned.trees[[i]])+1):(Ntip(pruned.trees[[i]])+Nnode(pruned.trees[[i]])) # Get node list
        for(j in nodestogo) { # For each node from consensus to be removed from MPT
            MPT.nodes<-MPT.nodes[-match(find.ancestor(get(paste("desc.",j,sep="")), pruned.trees[[i]]),MPT.nodes)] # Remove it
        }
        for(j in MPT.nodes) { # For each novel node in MPT
            desc.groups<-c(desc.groups,paste(sort(pruned.trees[[i]]$tip.label[find.descendants(j,pruned.trees[[i]])]),collapse="XQ")) # Add result to descendant groups
        }
    }
    new.nodes<-rle(sort(desc.groups))$lengths # Get clade frequencies across all MPTs
    names(new.nodes)<-rle(sort(desc.groups))$values # Label with clades, i.e. taxon clusters
    mrp.all<-cbind(mrp.all,matrix(rep(0,length(new.nodes)*length(mrp.all[,1])),ncol=length(new.nodes))) # Make space for results in mrp.all
    mrp.weights<-c(mrp.weights,rep(0,length(new.nodes))) # Make space for results in mrp.weights
    for(i in 1:length(new.nodes)) { # For each of the new nodes
        taxa<-strsplit(names(new.nodes[i]),"XQ")[[1]] # Get taxa that make up clade
        col<-(length(mrp.all[1,])+1)-i # Find appropraite column in MRP matrix and weights
        mrp.all[taxa,col]<-1 # Insert MRP data
        mrp.weights[col]<-new.nodes[i]/ntrees # Insert MRP results
    }
    out<-list(mrp.all,mrp.weights) # Create results list
    names(out)<-c("mrp","weights") # Add names
    return(out) # Return results
}

str.matrix<-function(clad.matrix) {
    full.matrix<-clad.matrix # Store extra copy of matrix
    pars.unif<-vector(mode="numeric") # Vector for storing constant characters
    for(i in 1:length(clad.matrix$matrix[1,])) { # For each character:
        if(length(unique(sort(clad.matrix$matrix[,i]))) <= 1) pars.unif[(length(pars.unif)+1)]<-i # Record constant, i.e. parsimony uninformative characters
    }
    zero.wts<-grep(TRUE,clad.matrix$weights == 0) # Record zero weight characters
    deletes<-sort(unique(c(pars.unif,zero.wts)))
    if(length(deletes) > 0) {
        clad.matrix$matrix<-clad.matrix$matrix[,-deletes]
        clad.matrix$ordering<-clad.matrix$ordering[-deletes]
        clad.matrix$weights<-clad.matrix$weights[-deletes]
        clad.matrix$max.vals<-clad.matrix$max.vals[-deletes]
        clad.matrix$min.vals<-clad.matrix$min.vals[-deletes]
    }
    dist.matrix<-dist.clad.matrix(clad.matrix, "euclidean")$dist.matrix # Get distance matrix
    pairs<-c(NA,NA) # Vector for stroing zero value pairs
    for(j in 1:length(dist.matrix[,1])) { # For each row
        for(k in 1:length(dist.matrix[,1])) { # For each column
            if(j > k) { # Make sure we are only looking at one diagonal
                if(is.na(dist.matrix[j,k]) == FALSE) { # For only distances that can be calculated
                    if(dist.matrix[j,k] == 0) { # If the distance is zero
                        pairs<-rbind(pairs,c(rownames(dist.matrix)[j],colnames(dist.matrix)[k])) # Record the taxon pairing
                    }
                }
            }
        }
    }
    rule.1a<-rule.1b<-rule.2a<-rule.2b<-c(NA,NA)
    if(length(pairs) > 2) { # As long as there are zerovalue pairs
        pairs<-pairs[-1,] # Delete first line which is empty
        if(length(pairs) == 2) pairs<-t(as.matrix(pairs)) # Case if only one pair
        for(j in 1:length(pairs[,1])) { # For each zero distance pair
            missing.1<-grep(TRUE,is.na(clad.matrix$matrix[match(pairs[j,1],rownames(clad.matrix$matrix)),])) # Get scored characters for first taxon
            missing.2<-grep(TRUE,is.na(clad.matrix$matrix[match(pairs[j,2],rownames(clad.matrix$matrix)),])) # Get scored characters for second taxon
            if(length(setdiff(missing.1,missing.2)) == 0 && length(setdiff(missing.2,missing.1)) == 0) { # STR Rule 1 test (equivalence; retain either taxon)
                if(length(missing.1) == length(clad.matrix$matrix[1,])) { # Meets Rule 1A (both taxa are known for all states)
                    rule.1a<-rbind(rule.1a,c(pairs[j,1],pairs[j,2]))
                } else { # Meets rule 1B (both taxa incompletely known)
                    rule.1b<-rbind(rule.1b,c(pairs[j,1],pairs[j,2]))
                }
            } else { # STR Rule 2 tests (asymmetric equivalence; retain more complete taxon)
                if(length(setdiff(missing.1,missing.2)) == 0) { # Taxon 2 redundant with respect to taxon 1
                    if(length(missing.1) == length(clad.matrix$matrix[1,])) { # If taxon 1 is completely known (Rule 2A)
                        rule.2a<-rbind(rule.2a,c(pairs[j,2],pairs[j,1])) # Add to vector, junior first
                    } else { # If taxon 1 is not completely known (Rule 2B)
                        rule.2b<-rbind(rule.2b,c(pairs[j,2],pairs[j,1])) # Add to vector, junior first
                    }
                }
                if(length(setdiff(missing.2,missing.1)) == 0) { # Taxon 1 redundant with respect to taxon 2
                    if(length(missing.2) == length(clad.matrix$matrix[1,])) { # If taxon 2 is completely known (Rule 2A)
                        rule.2a<-rbind(rule.2a,c(pairs[j,1],pairs[j,2])) # Add to vector, junior first
                    } else { # If taxon 2 is not completely known (Rule 2B)
                        rule.2b<-rbind(rule.2b,c(pairs[j,1],pairs[j,2])) # Add to vector, junior first
                    }
                }
            }
        }
    }
    # Remove empty first rows
    if(length(rule.1a) == 2) rule.1a<-rule.1a[-(1:2)] else rule.1a<-rule.1a[-1,]
	if(length(rule.1b) == 2) rule.1b<-rule.1b[-(1:2)] else rule.1b<-rule.1b[-1,]
    if(length(rule.2a) == 2) rule.2a<-rule.2a[-(1:2)] else rule.2a<-rule.2a[-1,]
    if(length(rule.2b) == 2) rule.2b<-rule.2b[-(1:2)] else rule.2b<-rule.2b[-1,]
    # Combine STR rules into single list:
    pairs<-rbind(rule.1a,rule.1b,rule.2a,rule.2b) # List taxon pairs
    rule<-c(rep("Rule 1A",length(rule.1a)/2),rep("Rule 1B",length(rule.1b)/2),rep("Rule 2A",length(rule.2a)/2),rep("Rule 2B",length(rule.2b)/2)) # List rules
    str.list<-cbind(pairs,rule) # Combine into single table
    if(length(str.list) > 0) {
        colnames(str.list)<-c("Junior","Senior","Rule") # Name columns
        removes<-sort(unique(str.list[,"Junior"])) # List junior taxa (i.e. those to remove)
        reduced.matrix<-removed.matrix<-full.matrix # New matrices for reduced and removed
        reduced.matrix<-reduced.matrix$matrix[-match(removes,rownames(reduced.matrix$matrix)),] # Remove STR taxa to create reduced matrix
        removed.matrix<-removed.matrix$matrix[match(removes,rownames(removed.matrix$matrix)),] # Isolate removed taxa for removed matrix
        result<-list(str.list,reduced.matrix,removed.matrix)
        names(result)<-c("str.list","reduced.matrix","removed.matrix")
    } else {
        result<-"No taxa can be safely removed"
    }
    return(result)
}

write.tnt<-function(clad.matrix, file.name) {
    head<-clad.matrix$header # Separate out header
    mat<-clad.matrix$matrix # Separate out matrix
    ord<-clad.matrix$ordering # Separate out ordering
    wts<-clad.matrix$weights # Separate out weights
    ntaxa<-length(mat[,1]) # Get number of taxa
    nchars<-length(mat[1,]) # Get number of characters
	head<-gsub("'","",head) # Remove any single quotes from header
    states<-sort(unique(as.vector(mat))) # Get states listed in matrix
    if(length(grep("&",states)) > 0) { # If polymorphisms are present
        while(length(grep("&",states)) > 0) { # For each polymorphism type
            temp<-strsplit(states[grep("&",states)[1]],"&")[[1]] # Split into separate components
            states<-states[-grep("&",states)[1]] # Remove type from list
            states<-unique(c(states,temp)) # Add components (if novel)
        }
    }
    if(length(states) > 8) { # If there are more than 8 states
        if(length(states) > 16) nstates<-32 # Set nstates at 32
		if(length(states) < 17) nstates<-16 # Set nstates at 16
    } else { # otherwise
        nstates<-8 # Set nstates at 8 (default)
    }
	if(max(nchar(rownames(clad.matrix$matrix))) > 32) { # Ifthe longest taxon name exceeds 32 characters
		taxname.length<-max(nchar(rownames(clad.matrix$matrix)))+1 # Set maximum taxon name length as that of maximum taxon
	} else { # Otherwise
		taxname.length<-32 # Just set it at 32 (TNT Ddefault)
	}
    if(head == "") { # If header text is absent:
        if(nstates == 8) toplines<-c("xread",paste(nchars,ntaxa,sep=" ")) # Make topline with no header and no states (8 or less)
        if(nstates == 16) toplines<-c("nstates num 16;","xread",paste(nchars,ntaxa,sep=" ")) # Make topline with no header and states (greater than 8, but less than 17)
        if(nstates == 32) toplines<-c("nstates num 32;","xread",paste(nchars,ntaxa,sep=" ")) # Make topline with no header and states (greater than 16, but less than 33)
    } else {
        if(nstates == 8) toplines<-c("xread",paste("'",paste(head,collapse="\n"),"'",sep=""),paste(nchars,ntaxa,sep=" ")) # Make topline with header and no states (8 or less)
        if(nstates == 16) toplines<-c("nstates num 16;","xread",paste("'",paste(head,collapse="\n"),"'",sep=""),paste(nchars,ntaxa,sep=" ")) # Make topline with header and states (greater than 8, but less than 17)
        if(nstates == 32) toplines<-c("nstates num 32;","xread",paste("'",paste(head,collapse="\n"),"'",sep=""),paste(nchars,ntaxa,sep=" ")) # Make topline with header and states (greater than 16, but less than 33)
    }
	if(taxname.length > 32) toplines<-c(paste("taxname +",taxname.length,";",sep=""),toplines) # Add extra line for extra-long taxon names
    taxa<-rownames(mat) # Get taxon names
    spaces.toadd<-(max(nchar(taxa))-nchar(taxa))+2 # Get number of spaces to add to character names
    for(i in 1:length(taxa)) taxa[i]<-paste(taxa[i],paste(rep(" ",spaces.toadd[i]),collapse=""),sep="") # Add spaces to taxon names
    # Convert alphabetic characters to numerical values:
    mat<-gsub("10","A",mat)
    mat<-gsub("11","B",mat)
    mat<-gsub("12","C",mat)
    mat<-gsub("13","D",mat)
    mat<-gsub("14","E",mat)
    mat<-gsub("15","F",mat)
    mat<-gsub("16","G",mat)
    mat<-gsub("17","H",mat)
    mat<-gsub("18","I",mat)
    mat<-gsub("19","J",mat)
    mat<-gsub("20","K",mat)
    mat<-gsub("21","L",mat)
    mat<-gsub("22","M",mat)
    mat<-gsub("23","N",mat)
    mat<-gsub("24","O",mat)
    mat<-gsub("25","P",mat)
    mat<-gsub("26","Q",mat)
    mat<-gsub("27","R",mat)
    mat<-gsub("28","S",mat)
    mat<-gsub("29","T",mat)
    mat<-gsub("30","U",mat)
    mat<-gsub("31","V",mat)
    mat<-gsub("32","W",mat)
    mat<-gsub("33","X",mat)
    mat<-gsub("34","Y",mat)
    mat<-gsub("35","Z",mat)
    if(length(grep("&",mat)) > 0) { # If there are polymorphic characters
        for(i in grep("&",mat)) { # For each polymorphic character
            mat[i]<-paste("[",gsub("&","",mat[i]),"]",sep="") # Convert to tnt format
        }
    }
    if(length(grep(TRUE,is.na(mat))) > 0) { # If missing characters are present
        missing<-grep(TRUE,is.na(mat)) # List missing values
        mat[missing]<-"?" # Replace with question mark
    }
    mtrx<-vector(mode="character") # Character matrix storage vector
    for(i in 1:ntaxa) { # For each taxon
        mtrx[i]<-paste(taxa[i],paste(mat[i,],collapse=""),sep="") # Store taxon lines in tnt format
    }
    ord<-gsub("UNORD","-",ord) # Replace UNORD with - for non-additive characters
    ord<-gsub("ORD","+",ord) # Replace ORD with + for additive characters
    ord.wts<-cbind(ord,wts) # Combine into single character block
    ord.wts<-paste(ord.wts[,1],ord.wts[,2],sep="[/") # Separate a la tnt format
    if(length(grep(TRUE,wts == 0)) > 0) { # Case if zero weight (inactive characters) are present
        for(i in grep(TRUE,wts == 0)) ord.wts[i]<-gsub("\\[/0","]/1",ord.wts[i]) # Make zero weights into inactive characters
    }
    charnos<-c(1:length(ord.wts))-1 # Get character numbers (starting from 0)
    ord.wts<-cbind(ord.wts,charnos) # Add character numbers to character block
    ord.wts<-paste(ord.wts[,1],ord.wts[,2],sep=" ") # Convert into vector
    ord.wts<-paste(ord.wts,collapse=" ") # Convert into single text string
    ccode<-paste("ccode ",ord.wts," ;",sep="") # Make character block
    tnt.file<-c(toplines,mtrx,";",ccode,"proc/;")
    write(tnt.file, file.name)
}

write.tnt.fa<-function(clad.matrix, file.name) {
    head<-clad.matrix$header # Separate out header
    mat<-clad.matrix$matrix # Separate out matrix
    ord<-clad.matrix$ordering # Separate out ordering
    wts<-clad.matrix$weights # Separate out weights
    ntaxa<-length(mat[,1]) # Get number of taxa
    nchars<-length(mat[1,]) # Get number of characters
	head<-gsub("'","",head) # Remove any single quotes from header
    states<-sort(unique(as.vector(mat))) # Get states listed in matrix
    if(length(grep("&",states)) > 0) { # If polymorphisms are present
        while(length(grep("&",states)) > 0) { # For each polymorphism type
            temp<-strsplit(states[grep("&",states)[1]],"&")[[1]] # Split into separate components
            states<-states[-grep("&",states)[1]] # Remove type from list
            states<-unique(c(states,temp)) # Add components (if novel)
        }
    }
    if(length(states) > 8) { # If there are more than 8 states
        if(length(states) > 16) nstates<-32 # Set nstates at 32
		if(length(states) < 17) nstates<-16 # Set nstates at 16
    } else { # otherwise
        nstates<-8 # Set nstates at 8 (default)
    }
	if(max(nchar(rownames(clad.matrix$matrix))) > 32) { # Ifthe longest taxon name exceeds 32 characters
		taxname.length<-max(nchar(rownames(clad.matrix$matrix)))+1 # Set maximum taxon name length as that of maximum taxon
	} else { # Otherwise
		taxname.length<-32 # Just set it at 32 (TNT Ddefault)
	}
    if(head == "") { # If header text is absent:
        if(nstates == 8) toplines<-c("xread",paste(nchars,ntaxa,sep=" ")) # Make topline with no header and no states (8 or less)
        if(nstates == 16) toplines<-c("nstates num 16;","xread",paste(nchars,ntaxa,sep=" ")) # Make topline with no header and states (greater than 8, but less than 17)
        if(nstates == 32) toplines<-c("nstates num 32;","xread",paste(nchars,ntaxa,sep=" ")) # Make topline with no header and states (greater than 16, but less than 33)
    } else {
        if(nstates == 8) toplines<-c("xread",paste("'",paste(head,collapse="\n"),"'",sep=""),paste(nchars,ntaxa,sep=" ")) # Make topline with header and no states (8 or less)
        if(nstates == 16) toplines<-c("nstates num 16;","xread",paste("'",paste(head,collapse="\n"),"'",sep=""),paste(nchars,ntaxa,sep=" ")) # Make topline with header and states (greater than 8, but less than 17)
        if(nstates == 32) toplines<-c("nstates num 32;","xread",paste("'",paste(head,collapse="\n"),"'",sep=""),paste(nchars,ntaxa,sep=" ")) # Make topline with header and states (greater than 16, but less than 33)
    }
	if(taxname.length > 32) toplines<-c(paste("taxname +",taxname.length,";",sep=""),toplines) # Add extra line for extra-long taxon names
    taxa<-rownames(mat) # Get taxon names
    spaces.toadd<-(max(nchar(taxa))-nchar(taxa))+2 # Get number of spaces to add to character names
    for(i in 1:length(taxa)) taxa[i]<-paste(taxa[i],paste(rep(" ",spaces.toadd[i]),collapse=""),sep="") # Add spaces to taxon names
	# Convert alphabetic characters to numerical values:
    mat<-gsub("10","A",mat)
    mat<-gsub("11","B",mat)
    mat<-gsub("12","C",mat)
    mat<-gsub("13","D",mat)
    mat<-gsub("14","E",mat)
    mat<-gsub("15","F",mat)
    mat<-gsub("16","G",mat)
    mat<-gsub("17","H",mat)
    mat<-gsub("18","I",mat)
    mat<-gsub("19","J",mat)
    mat<-gsub("20","K",mat)
    mat<-gsub("21","L",mat)
    mat<-gsub("22","M",mat)
    mat<-gsub("23","N",mat)
    mat<-gsub("24","O",mat)
    mat<-gsub("25","P",mat)
    mat<-gsub("26","Q",mat)
    mat<-gsub("27","R",mat)
    mat<-gsub("28","S",mat)
    mat<-gsub("29","T",mat)
    mat<-gsub("30","U",mat)
    mat<-gsub("31","V",mat)
    mat<-gsub("32","W",mat)
    mat<-gsub("33","X",mat)
    mat<-gsub("34","Y",mat)
    mat<-gsub("35","Z",mat)
    if(length(grep("&",mat)) > 0) { # If there are polymorphic characters
        for(i in grep("&",mat)) { # For each polymorphic character
            mat[i]<-paste("[",gsub("&","",mat[i]),"]",sep="") # Convert to tnt format
        }
    }
    if(length(grep(TRUE,is.na(mat))) > 0) { # If missing characters are present
        missing<-grep(TRUE,is.na(mat)) # List missing values
        mat[missing]<-"?" # Replace with question mark
    }
    mtrx<-vector(mode="character") # Character matrix storage vector
    for(i in 1:ntaxa) { # For each taxon
        mtrx[i]<-paste(taxa[i],paste(mat[i,],collapse=""),sep="") # Store taxon lines in tnt format
    }
    ord<-gsub("UNORD","-",ord) # Replace UNORD with - for non-additive characters
    ord<-gsub("ORD","+",ord) # Replace ORD with + for additive characters
    ord.wts<-cbind(ord,wts) # Combine into single character block
    ord.wts<-paste(ord.wts[,1],ord.wts[,2],sep="[/") # Separate a la tnt format
    if(length(grep(TRUE,wts == 0)) > 0) { # Case if zero weight (inactive characters) are present
        for(i in grep(TRUE,wts == 0)) ord.wts[i]<-gsub("\\[/0","]/1",ord.wts[i]) # Make zero weights into inactive characters
    }
    charnos<-c(1:length(ord.wts))-1 # Get character numbers (starting from 0)
    ord.wts<-cbind(ord.wts,charnos) # Add character numbers to character block
    ord.wts<-paste(ord.wts[,1],ord.wts[,2],sep=" ") # Convert into vector
    ord.wts<-paste(ord.wts,collapse=" ") # Convert into single text string
    ccode<-paste("ccode ",ord.wts," ;",sep="") # Make character block
	# Everything above is as in write.tnt(); below is the novel part:
    if(ntaxa <= 24) { # If there are few enough taxa for an exact solution
        anal<-c("collapse [;","ienum;") # Use implicit enumeration
    } else { # Otherwise
        anal<-c("rseed*;","hold 999;","xmult=rss fuse 50 drift 50 ratchet 50;","mult 50 =tbr drift;","tsave scratch.tre;","save;","tsave /;") # First iteration with new tech
		anal<-c(anal,rep(c("rseed*;","hold 999;","xmult=rss fuse 50 drift 50 ratchet 50;","mult 50 =tbr drift;","tsave scratch.tre +;","save;","tsave /;"),19)) # Iterations 2-20
		anal<-c(anal,c("hold 1000000;","shortread scratch.tre;","bbreak=tbr;")) # Read in new technology trees and finish with a heuristic (tbr) search for all mpts
    }
    out.file<-strsplit(strsplit(file.name,"/")[[1]][length(strsplit(file.name,"/")[[1]])],"\\.")[[1]][1] # Get stripped file name for use in export lines
    strict.name<-paste("export -",out.file,"tntmpts_plus_strict.nex;",sep="")
    mrp.name<-c("mrp;",paste("export ",out.file,"mrp.nex;",sep=""))
    tnt.file<-c("mxram 1024;",toplines,mtrx,";",ccode,anal,"nelsen*;",strict.name,mrp.name,"proc/;")
    write(tnt.file, file.name)
}

write.nexus.data<-function(clad.matrix, filename)
{
    topline<-"#NEXUS" # Top line
    if(!clad.matrix$header == "") { # If header text is present
        headlines<-vector(mode="character")
        for(i in length(clad.matrix$header)) { # For each line of header text
            headlines[(length(headlines)+1)]<-paste("[!",clad.matrix$header[1],"]",sep="") # Store as NEXUS text (i.e. non-readable)
        }
        headlines<-c("",headlines,"") # Fill header lines
    } else { # If no text
        headlines<-""  # Fill header line with a blank
    }
    states<-sort(unique(as.vector(clad.matrix$matrix))) # Get states in matrix
    if(length(grep("&",states)) > 0) { # If polymorphisms are present
        for(i in grep("&",states)) states<-c(states,strsplit(states[i],"&")[[1]]) # Find and split them
        states<-states[-grep("&",states)] # Remove them from states list
    }
    n.states<-max(sort(as.numeric(unique(states))))+1 # Get total number of states
    if(n.states > 10) { # If there are more states than can be represented by numbers (0-9)
        states<-c(0:9,LETTERS[1:(n.states-10)]) # List states as numbers and letters
        for(i in 10:(n.states-1)) { # For each number that needs to be replaced by a letter
            clad.matrix$matrix<-gsub(as.character(i),LETTERS[(i-9)],clad.matrix$matrix) # Replace it with a letter
        }
    } else {
        states<-0:(n.states-1) # List states as numbers
    }
    datablock<-c("BEGIN DATA;",paste("\t","DIMENSIONS  NTAX=",length(clad.matrix$matrix[,1])," NCHAR=",length(clad.matrix$matrix[1,]),";",sep=""),paste("\t","FORMAT SYMBOLS= \" ",paste(states,collapse=" "),"\" MISSING=? GAP=- ;",sep="")) # Fill out data block
    matrixblock.top<-c("MATRIX","") # Make top of matrix block
    name.spaces<-(max(nchar(rownames(clad.matrix$matrix)))+2)-nchar(rownames(clad.matrix$matrix)) # Get spaces required after each taxon name
    if(length(grep("&",clad.matrix$matrix)) > 0) { # If polymorphisms are present in the matrix
        for(i in grep("&",clad.matrix$matrix)) {
            clad.matrix$matrix[i]<-paste("(",paste(strsplit(clad.matrix$matrix[i],"&")[[1]],collapse=""),")",sep="")
        }
    }
    clad.matrix$matrix[grep(TRUE,is.na(clad.matrix$matrix))]<-"?" # Replace NAs with question marks
    matrixblock<-vector(mode="character") # Matrix block proper
    for(i in 1:length(clad.matrix$matrix[,1])) { # For each taxon
        matrixblock[i]<-paste(rownames(clad.matrix$matrix)[i],paste(rep(" ",name.spaces[i]),collapse=""),paste(clad.matrix$matrix[i,],collapse=""),sep="") # Get name, spaces and states
    }
    matrixblock.bottom<-c(";","END;","")
    if(length(unique(clad.matrix$ordering)) == 1 && unique(clad.matrix$ordering) == "UNORD") ord.block<-"\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;" # Case if all characters are unordered
    if(length(unique(clad.matrix$ordering)) == 1 && unique(clad.matrix$ordering) == "ORD") ord.block<-"\tOPTIONS  DEFTYPE=ord PolyTcount=MINSTEPS ;" # Case if all characters are ordered
    if(length(unique(clad.matrix$ordering)) == 2) { # Case if characters are a mix of ordered and unordered
        ord<-grep(TRUE,clad.matrix$ordering == "ORD") # List ordered characters
        unord<-grep(TRUE,clad.matrix$ordering == "UNORD") # List unordered characters
        if(min(diff(ord)) == 1) { # If there are contiguous ordered characters
            sort.ord<-vector(mode="numeric")
            ord.groups<-paste(ord[grep(TRUE,diff(ord) == 1)],ord[grep(TRUE,diff(ord) == 1)+1],sep="-") # Find contiguous characters
            ord<-ord[-sort(unique(c((grep(TRUE,diff(ord) == 1)+1),grep(TRUE,diff(ord) == 1))))] # Remove from series
            if(length(ord.groups) > 1) { # Nothing to collapse if only 1 group!
			   for(i in length(ord.groups):2) {
					if(diff(c(as.numeric(strsplit(ord.groups[(i-1)],"-")[[1]][2]),as.numeric(strsplit(ord.groups[i],"-")[[1]][1]))) == 0) {
						ord.groups[(i-1)]<-paste(strsplit(ord.groups[(i-1)],"-")[[1]][1],strsplit(ord.groups[i],"-")[[1]][2],sep="-")
						ord.groups<-ord.groups[-i]
					}
			   }
			}
            ord<-c(ord,ord.groups)
            for(i in 1:length(ord)) sort.ord[i]<-as.numeric(strsplit(ord[i],"-")[[1]][1])
            ord<-ord[order(sort.ord)]
        }
        if(min(diff(unord)) == 1) { # If there are contiguous ordered characters
            sort.unord<-vector(mode="numeric")
            unord.groups<-paste(unord[grep(TRUE,diff(unord) == 1)],unord[grep(TRUE,diff(unord) == 1)+1],sep="-") # Find contiguous characters
            unord<-unord[-sort(unique(c((grep(TRUE,diff(unord) == 1)+1),grep(TRUE,diff(unord) == 1))))] # Remove from series
			if(length(unord.groups) > 1) { # Nothing to collapse if only 1 group!
				for(i in length(unord.groups):2) {
					if(diff(c(as.numeric(strsplit(unord.groups[(i-1)],"-")[[1]][2]),as.numeric(strsplit(unord.groups[i],"-")[[1]][1]))) == 0) {
						unord.groups[(i-1)]<-paste(strsplit(unord.groups[(i-1)],"-")[[1]][1],strsplit(unord.groups[i],"-")[[1]][2],sep="-")
						unord.groups<-unord.groups[-i]
					}
                }
            }
            unord<-c(unord,unord.groups)
            for(i in 1:length(unord)) sort.unord[i]<-as.numeric(strsplit(unord[i],"-")[[1]][1])
            unord<-unord[order(sort.unord)]
        }
        ord.block<-c("\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;",paste("\tTYPESET * UNTITLED  = ",paste("unord: ",paste(unord,collapse=" "),sep=""),", ",paste("ord: ",paste(ord,collapse=" "),sep=""),";",sep=""))
    }
    if(length(grep(TRUE,unique(clad.matrix$weights) != 1)) > 0 || length(unique(clad.matrix$weights)) > 1) { # Case if not all weights are equal to 1
        wts<-sort(unique(clad.matrix$weights)) # Record weight numbers found
        wts.line<-vector(mode="character") # Make vector for storage
        for(i in wts) { # For each weight
            weight.nos<-grep(TRUE,clad.matrix$weights == i) # Get character numbers which have that weight
            if(min(diff(weight.nos)) == 1) { # If there are contiguous characters with the same weight
                sort.wt<-vector(mode="numeric") # Vector for sorting (used at end)
                wt.groups<-paste(weight.nos[grep(TRUE,diff(weight.nos) == 1)],weight.nos[grep(TRUE,diff(weight.nos) == 1)+1],sep="-") # Find contiguous characters
                weight.nos<-weight.nos[-sort(unique(c((grep(TRUE,diff(weight.nos) == 1)+1),grep(TRUE,diff(weight.nos) == 1))))] # Remove from series
                for(j in length(wt.groups):2) { # For each adjacent pair
                    if(diff(c(as.numeric(strsplit(wt.groups[(j-1)],"-")[[1]][2]),as.numeric(strsplit(wt.groups[j],"-")[[1]][1]))) == 0) { # If they overlap
                        wt.groups[(j-1)]<-paste(strsplit(wt.groups[(j-1)],"-")[[1]][1],strsplit(wt.groups[j],"-")[[1]][2],sep="-") # Collapse them into one
                        wt.groups<-wt.groups[-j] # And remove the now redundant one
                    }
                }
                weight.nos<-c(weight.nos,wt.groups) # Group contiguous with non-contiguous
                for(j in 1:length(weight.nos)) sort.wt[j]<-as.numeric(strsplit(weight.nos[j],"-")[[1]][1]) # Get correct order
                weight.nos<-weight.nos[order(sort.wt)] # And sort
            }
            wts.line<-c(wts.line,paste(i,": ",paste(weight.nos,collapse=" "),sep="")) # Store each weight line
        }
        ord.block<-c(ord.block,paste("\tWTSET * UNTITLED  = ",paste(wts.line,collapse=", "),";",sep="")) # Add to ordering block
    }
    ord.block<-c("BEGIN ASSUMPTIONS;",ord.block,"END;") # Top and tail ordering block
    out<-c(topline,headlines,datablock,matrixblock.top,matrixblock,matrixblock.bottom,ord.block)
    write(out, filename)
}

tntmrp.to.nexusmrp<-function(file)
{
	# Read in raw file:
	X<-scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
	start<-grep("ROOT",X)+1
	end<-grep(";",X[start:length(X)])[1]+start-2
	MATRIX<-X[start:end]
	while(length(grep("  ",MATRIX)) > 0) MATRIX<-gsub("  "," ",MATRIX)
	characters<-names<-vector(mode="character")
	for(i in 1:length(MATRIX)) {
		names[i]<-strsplit(MATRIX[i]," ")[[1]][1]
		characters[i]<-strsplit(MATRIX[i]," ")[[1]][2]
	}
	char.block<-matrix(nrow=length(names),ncol=nchar(characters[1]))
	for(i in 1:length(MATRIX)) char.block[i,]<-strsplit(characters[i],"")[[1]]
	char.block<-t(char.block)
	characters<-vector(mode="character")
	for(i in 1:length(char.block[,1])) characters[i]<-paste(char.block[i,],collapse="")
	characters<-unique(characters)
	char.block<-char.block[1:length(characters),]
	for(i in 1:length(characters)) char.block[i,]<-strsplit(characters[i],"")[[1]]
	char.block<-t(char.block)
	rownames(char.block)<-names
	# Make into clad.matrix format:
	header<-""
	ordering<-rep("UNORD",length(char.block[1,]))
	weights<-rep(1,length(char.block[1,]))
	max.vals<-rep(1,length(char.block[1,]))
	min.vals<-rep(0,length(char.block[1,]))
	result<-list(header,char.block,ordering,weights,max.vals,min.vals)
	names(result)<-c("header","matrix","ordering","weights","max.vals","min.vals")
	return(result)
}

read.tnt<-function(file)
{
    X<-scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in NEXUS file
    if(length(grep("mxram",X)) > 0) { # If there is an nstates line
        X<-X[-grep("mxram",X)] # Remove nstates line
    }
    if(length(grep("nstates",X)) > 0) { # If there is an nstates line
        X<-X[-grep("nstates",X)] # Remove nstates line
    }
    if(length(grep("taxname",X)) > 0) { # If there is a taxname line
        X<-X[-grep("taxname",X)] # Remove taxname line
    }
    X<-X[-grep(TRUE,X == "xread")] # Remove xread
    if(length(grep(",",X)) > 0) {
        header<-X[grep(",",X)] # Find header text
        header<-gsub("'","",header) # Remove quotes
        X<-X[-grep(",",X)] # Remove header line(s)
    } else {
        header<-""
    }
    nchar<-as.numeric(strsplit(X[1]," ")[[1]][1]) # Get number of characters
    ntax<-as.numeric(strsplit(X[1]," ")[[1]][2]) # Get number of taxa
    X<-X[-1] # Delete top line
    matrixblock<-X[1:ntax] # Get matrix block
    X<-X[-(1:ntax)] # Remove matrix block
    MATRIX<-matrix(nrow=ntax,ncol=nchar) # Make matrix
    rownames(MATRIX)<-c(1:ntax) # Dummy rownames
    for(i in 1:ntax) { # For each taxon
        rownames(MATRIX)[i]<-strsplit(matrixblock[i]," ")[[1]][1] # Extract taxon name
        matrixblock[i]<-strsplit(matrixblock[i]," ")[[1]][length(strsplit(matrixblock[i]," ")[[1]])] # And character block
        chars<-strsplit(matrixblock[i],"")[[1]] # Get character vector
        j<-1 # Start value for j
        while(j <= length(chars)) { # Go through character vector until all characters have been looked at
            if(chars[j] != "[") { # If the character is not a polymorphism
                MATRIX[i,grep(TRUE,is.na(MATRIX[i,]))[1]] <- chars[j] # Just store it in the matrix
                j<-j+1 # And move to the next one
			} else { # If it is a polymorphism
                MATRIX[i,grep(TRUE,is.na(MATRIX[i,]))[1]] <- paste(chars[(j+1):(j+grep(TRUE,chars[j:length(chars)] == "]")[1]-2)],collapse="&") # Paste all states into the matrix
                j<-j+grep(TRUE,chars[j:length(chars)] == "]")[1] # And move to the next character
			}
		}
	}
	MATRIX<-gsub("\\?",NA,MATRIX) # Replace question marks with NAs
	unq.states<-sort(unique(as.vector(MATRIX))) # Find unique states list
	if(length(grep("&",unq.states)) > 0) { # If polymorphisms are present
		polys<-unq.states[grep("&",unq.states)] # Isolate them
		unq.states<-unq.states[-grep("&",unq.states)] # Remove them
		for(i in 1:length(polys)) { # Go through each polymorphism
			unq.states<-c(unq.states,strsplit(polys[i],"&")[[1]]) # Split adn add to unqiue states list
		}
		unq.states<-sort(unique(unq.states)) # Get unique states again
	}
	if(length(sort(match(LETTERS,unq.states))) > 0) { # If there are letters (greater than ten states)
		lets.found<-unq.states[sort(match(LETTERS,unq.states))] # List leters found
		for(i in 1:length(lets.found)) { # For each leter found
			MATRIX<-gsub(lets.found[i],match(lets.found[i],LETTERS)+10,MATRIX) # Replace with appropriate number
		}
	}
	if(length(grep("ccode",X)) > 0) { # If weights and/or orderings are present
		ordering<-weights<-vector(mode="character",length=nchar)
		wt.ord<-gsub(";","",gsub("ccode ","",X[grep("ccode",X)])) # Get character weights/orderings
		wt.ord<-strsplit(wt.ord," ")[[1]]
		char.nos<-as.numeric(wt.ord[c(1:length(wt.ord))[grep(TRUE,c(1:length(wt.ord)) %% 2 == 0)]])+1 # Get character numbers (even cells)
		char.vals<-wt.ord[c(1:length(wt.ord))[grep(TRUE,c(1:length(wt.ord)) %% 2 != 0)]] # Get character values (odd cells)
		ordering[grep("-",char.vals)]<-"UNORD" # Get unordered characters
		ordering[grep("\\+",char.vals)]<-"ORD" # Get ordered characters
		for(i in 1:nchar) { # For each character
			if(length(grep("\\[",char.vals[i])) > 0) { # If character is active
                weights[i]<-as.numeric(strsplit(char.vals[i],"/")[[1]][2]) # Set weight
			} else { # If character is inactive
                weights[i]<-0 # Set weight to zero
			}
		}
	} else { # If they are not
		weights<-rep(1,nchar) # Make all characters weighted one
		ordering<-rep("UNORD",nchar) # Make all characters unordered
	}
	max.vals<-min.vals<-vector(mode="numeric",length=nchar) # Min and max value vectors
	for(i in 1:nchar) { # For each character
		unq.states<-sort(unique(as.vector(MATRIX[,i])))
		if(length(grep("&",unq.states)) > 0) { # If polymorphisms are present
			polys<-unq.states[grep("&",unq.states)] # Isolate them
			unq.states<-unq.states[-grep("&",unq.states)] # Remove them
			for(j in 1:length(polys)) { # Go through each polymorphism
				unq.states<-c(unq.states,strsplit(polys[j],"&")[[1]]) # Split and add to unqiue states list
			}
			unq.states<-sort(unique(unq.states)) # Get unique states again
		}
		min.vals[i]<-min(as.numeric(unq.states))
		max.vals[i]<-max(as.numeric(unq.states))
	}
	result<-list(header,MATRIX,ordering,weights,max.vals,min.vals)
	names(result)<-c("header","matrix","ordering","weights","max.vals","min.vals")
	return(result)
}

trees.to.mpts.plus.strict<-function(file)
{
	X<-scan(file = file, what = "", sep = "\n", quiet = TRUE) # Read in tnt treefile
	X<-X[grep("\\(",X)] # Get just the trees
	X<-gsub(" ","",X) # Remove spaces
	strict<-X[length(X)] # Get strict
	mpts<-X[1:(length(X)-1)] # Get mpts
	result<-list(mpts,strict) # Make output variable
	names(result)<-c("mpts","strict") # Add names
	return(result) # Return
}

pbdb.occs.in<-function(file,sep="\t",header=TRUE)
{
	Sys.setlocale('LC_ALL','C')
	X<-scan(file = file, what = "", sep = "\n", quiet = TRUE)
	for(i in 1:length(X)) {
		new.block<-strsplit(X[i],"\"")[[1]] # Break by quotation marks
		evens<-grep(TRUE,c(1:length(new.block)) %% 2 == 0) # Find within quotes parts (i.e. even numbers)
		new.block[evens]<-gsub(sep,"%sEpArAtOr%",new.block[evens]) # Replace seps within quotes
		X[i]<-paste(new.block,collapse="\"") # Reinsert quotes and store in X
	}
	seps<-vector(mode="numeric")
	for(i in 1:length(X)) {
		seps[i]<-length(grep(TRUE,strsplit(X[i],"")[[1]] == sep))
	}
	if(header == TRUE) {
		out<-matrix(ncol=min(seps)+1,nrow=length(X)-1)
		colnames(out)<-strsplit(X[1],sep)[[1]]
		X<-X[-1]
	}
	if(header == FALSE) out<-matrix(ncol=min(seps)+1,nrow=length(X))
	for(i in 1:length(X)) out[i,]<-strsplit(X[i],sep)[[1]]
	out<-gsub("%sEpArAtOr%",sep,out) # Reinsert non-separating separator characters
	out<-gsub("\"","",out) # Remove now redundant quotes
	return(out)
}