pgs.cols <- setNames(brewer.pal(length(pgs[[rgn]][1,,1]), "Set3"), names(pgs[[rgn]][1,,1]))

if (clade == "tetrapoda" & use.stages == TRUE) {
	focal_intervals <- rev(which(intervals$bin %in% c("Ypresian","Danian","Maastrichtian")))
} else if (clade == "tetrapoda" & use.stages == FALSE) {
	focal_intervals <- rev(which(intervals$bin %in% c("K8","Pg0","Pg2"))) #"K8","Pg0","Pg2"
} else if (clade == "dinosauromorpha" & use.stages == TRUE) {
	focal_intervals <- rev(which(intervals$bin %in% c("Tithonian","Campanian","Maastrichtian")))
} else if (clade == "dinosauromorpha" & use.stages == FALSE) {
	focal_intervals <- rev(which(intervals$bin %in% c("Tr4","J6","K8")))
}
space <- 1.2
legend_inset <- c(-0.65,0)
legend.cex <- 1.1
label.cex <- 0.8

#### Figures 1/3 (tetrapods and dinosaurs, respectively): SQS subsampled diversity adding occurrences ####
#run the next two lines of code to get a MEGAPLOT!!!
# focal_intervals <- which(sapply(numoccs[[rgn]], max) > 40)
# focal_intervals <- focal_intervals[which(names(focal_intervals) %in% intervals[pub_ints,'bin'])]
plot_SQS1 <- T
plot_SQS2 <- T
plot_iNEXT <- F
plot_TRiPS <- F
plot_samplingrate <- T
plot_completeness <- T
plot_ntons <- T
plot_pgs <- T
plot_extrapolators <- T
plot_CR <- F
fig1rows <- sum(plot_SQS1,plot_SQS2,plot_iNEXT,plot_CR,plot_TRiPS,plot_extrapolators,plot_samplingrate,plot_completeness,plot_ntons,plot_pgs)
fig1columns <- length(focal_intervals)
plot.nulls <- F

pdf(paste(folder.name, "/plots/Figure <", clade, "-diversity-plots> use.stages=", use.stages, ".pdf", sep = ""), width = fig1columns*3.3, height = fig1rows*1.5)
par(mfrow = c(fig1rows, fig1columns))
par(cex = 0.7)
par(mar = c(0, 3, 2, 0))
par(oma = c(5, 1, 0.5, 11))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
cex.lab <- 0.7
vlaboffset <- 0.09
hlaboffset <- 2.5
legend_inset <- c(-0.5,0)
legend.cex <- 1
colour.lwd <- 1.5; black.lwd <- colour.lwd*1.2
border <- 1 #1 = yes, 2 = no
border.lwd <- 0.05
null.lwd <- 10

myletters <- myLetters(length(focal_intervals) * fig1rows)
# myletters <- NA

for (i in focal_intervals) {

	xlim <- c(1,max(numoccs[[rgn]][[i]]))

	p <- 0

	if (run_SQS_Perl_1 == T && plot_SQS1 == T) {
		#plot subsampled diversity (SQS with all the trimmings) ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), function(x) {
			tmp <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[x], flavour == "SQS1")
			tmp$Subsampled.diversity/mean(tmp$Subsampled.diversity, na.rm = T)})))

		if (nulldist == TRUE && plot.nulls == TRUE) {
			sqs.null.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), function(q) {
				null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,] / mean(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,])
			})))
			sqs.ylim <- ylim.max(c(sqs.ylim, sqs.null.ylim))
		}
		pgs.ylim <- ylim.max(pgs[[rgn]][i,"mst_dist",])

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(0,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab)
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(SQS Perl Script V1)", side = 2, line = 2, cex = cex.lab)}
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (q in 1:length(quorum.levels)) {

				for (r in 1:nshuff) {
					lines((null.sqs.1[[rgn]][i,"Subsampled.diversity",q,r,] / mean(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,r,], na.rm = T)) ~ numoccs[[rgn]][[i]], col = adjustcolor(quorum.cols[q], alpha.f = 0.05), lwd = null.lwd)
				}
			}
		}
		for (j in 1:length(quorum.levels)) {
			X <- filter(df, region == names(regions)[rgn], bin == intervals$bin[i])
			Y <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[j], flavour == "SQS1")

			for (k in border:2) {points(X$numoccs, Y$Subsampled.diversity/mean(Y$Subsampled.diversity, na.rm = T), type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white", title = "Quorum Level")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex, title = "Quorum Level")
			par(xpd = FALSE)
		}
	}

	if (run_SQS_Perl_2 == T && plot_SQS2 == T) {
		#plot subsampled diversity (plain-vanilla SQS) ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), FUN = function(x) real.sqs.2[[rgn]][i,"Subsampled.diversity",x,]/mean(real.sqs.2[[rgn]][i,"Subsampled.diversity",x,], na.rm = T))))
		pgs.ylim <- ylim.max(pgs[[rgn]][i,"mst_dist",])

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(0,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, main = "")
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(SQS Perl Script V2)", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (q in 1:length(quorum.levels)) {
				lower_ci <- sapply(1:length(times), function(j) min(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,j]))
				upper_ci <- sapply(1:length(times), function(j) max(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,j]))
				noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
				lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
				polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[q], alpha.f = 0.2), border = FALSE, lty = 2, lwd = 0.5)

				for (r in 1:nshuff) {
					lines(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,r,] ~ numoccs[[rgn]][[i]], col = adjustcolor(quorum.cols[q], alpha.f = 0.05), lwd = null.lwd)
				}
			}
		}
		for (j in 1:length(quorum.levels)) {
			X <- filter(df, region == names(regions)[rgn], bin == intervals$bin[i])
			Y <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[j], flavour == "SQS2")

			for (k in border:2) {points(X$numoccs, Y$Subsampled.diversity/mean(Y$Subsampled.diversity, na.rm = T), type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (run_iNEXT == T && plot_iNEXT == T) {
		#plot subsampled diversity (iNEXT) ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), FUN = function(x) filter(iNEXT_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[x], order == 0, method != "extrapolated")$qD.UCL/mean(filter(iNEXT_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[x], order == 0, method != "extrapolated")$qD.UCL, na.rm = T))))
		pgs.ylim <- ylim.max(pgs[[rgn]][i,"mst_dist",])

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(0,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, main = "")
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(iNEXT CBR)", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (q in 1:length(quorum.levels)) {
				lower_ci <- sapply(1:length(times), function(j) min(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,j]))
				upper_ci <- sapply(1:length(times), function(j) max(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,j]))
				noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
				lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
				polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[q], alpha.f = 0.2), border = FALSE, lty = 2, lwd = 0.5)

				for (r in 1:nshuff) {
					lines(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,r,] ~ numoccs[[rgn]][[i]], col = adjustcolor(quorum.cols[q], alpha.f = 0.05), lwd = null.lwd)
				}
			}
		}
		for (j in 1:length(quorum.levels)) {
			X <- filter(df, region == names(regions)[rgn], bin == intervals$bin[i])
			Y2 <- filter(iNEXT_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[j], order == 0, method != "extrapolated")
			X <- left_join(Y2, df)

			lower_ci <- Y2$qD.LCL/mean(Y2$qD.LCL, na.rm = T)
			upper_ci <- Y2$qD.UCL/mean(Y2$qD.UCL, na.rm = T)
			lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
			polygon(c(X$numoccs[which(!is.na(lower_ci))], rev(X$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[j], alpha.f = 0.1), border = TRUE, lty = 2, lwd = border.lwd)
			for (k in border:2) {points(X$numoccs, Y2$qD/mean(Y2$qD, na.rm = T), type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (run_CR == T && plot_CR == T) {
		#plot iNEXT CR ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quota.levels), function(x) {
			tmp <- filter(CR_df, region == names(regions)[rgn], bin == intervals$bin[i], quota == names(quota.levels)[x])
			tmp$qD/mean(tmp$qD, na.rm = T)}
		)))

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(0,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, main = "")
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(iNEXT CR))", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (q in 1:length(quorum.levels)) {
				lower_ci <- sapply(1:length(times), function(j) min(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,j]))
				upper_ci <- sapply(1:length(times), function(j) max(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,j]))
				noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
				lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
				polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[q], alpha.f = 0.2), border = FALSE, lty = 2, lwd = 0.5)

				for (r in 1:nshuff) {
					lines(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,r,] ~ numoccs[[rgn]][[i]], col = adjustcolor(quorum.cols[q], alpha.f = 0.05), lwd = null.lwd)
				}
			}
		}
		X <- filter(df, region == names(regions)[rgn], bin == intervals$bin[i])

		for (j in 1:length(quota.levels)) {
			Y2 <- filter(CR_df, region == names(regions)[rgn], bin == intervals$bin[i], quota == names(quota.levels)[j], order == 0)
			lower_ci <- Y2$qD.LCL/mean(Y2$qD.LCL, na.rm = T)
			upper_ci <- Y2$qD.UCL/mean(Y2$qD.UCL, na.rm = T)
			lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
			polygon(c(X$numoccs[which(!is.na(lower_ci))], rev(X$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[j], alpha.f = 0.1), border = TRUE, lty = 2, lwd = border.lwd)
			for (k in border:2) {points(X$numoccs, Y2$qD/mean(Y2$qD, na.rm = T), type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (plot_extrapolators == T) {
		Y1 <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[1], flavour == "SQS2")
		Y2 <- filter(lambda5_df, region == names(regions)[rgn], bin == intervals$bin[i])
		Y3 <- filter(chao_df, region == names(regions)[rgn], bin == intervals$bin[i])
		X1 <- left_join(Y1, df)
		X2 <- left_join(Y2, df)
		X3 <- left_join(Y3, df)
		trips.upper.ylim <- ylim.max(real.trips[[rgn]]["Estimated richness","upper CI",i,])
		chao2.ylim <- ylim.max(c(Y1$Chao.2, Y2$lambda5.UCL))
		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(0, chao2.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == focal_intervals[1]) {mtext("Species Richness", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#lambda5
		lower_ci <- Y2$lambda5.LCL
		upper_ci <- Y2$lambda5.UCL
		lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
		polygon(c(X2$numoccs[which(!is.na(lower_ci))], rev(X2$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor("red", alpha.f = 0.1), border = TRUE, lty = 2, lwd = border.lwd)
		points(X2$numoccs, Y2$lambda5, type = "l", lty = 2, pch = 21, cex = 0.5, col = "red", lwd = colour.lwd)

		#chao
		lower_ci <- Y3$chao.LCL
		upper_ci <- Y3$chao.UCL
		lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
		polygon(c(X3$numoccs[which(!is.na(lower_ci))], rev(X3$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor("grey", alpha.f = 0.1), border = TRUE, lty = 2, lwd = border.lwd)
		points(X3$numoccs, Y3$chao, type = "l", lty = 2, pch = 21, cex = 0.5, col = "black", lwd = colour.lwd)

		#TRiPS
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (j in 1:nshuff) {
				lines(null.trips[[1]]["Estimated richness","MLE",i,j,] ~ numoccs[[rgn]][[i]], col = adjustcolor("blue", alpha.f = 0.1), lwd = null.lwd)
			}
		}
		lower_ci <- real.trips[[rgn]]["Estimated richness","lower CI",i,]
		upper_ci <- real.trips[[rgn]]["Estimated richness","upper CI",i,]
		noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
		lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
		polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor("blue", alpha.f = 0.1), border = TRUE, lty = 2, lwd = 0.5)
		points(numoccs[[rgn]][[i]], real.trips[[rgn]]["Estimated richness","MLE",i,], type = "l", lty = 2, pch = 21, cex = 0.5, col = "blue", lwd = colour.lwd)

		#legend
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = c("red","black","red","blue"), legend = c("Raw","Chao2","lambda-5","TRiPS"), lty = c(1,2,2,2), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (plot_TRiPS == T) {
		#### plot raw diversity ####
		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		raw.ylim <- ylim.max(real.sqs[[rgn]][i,"Raw.species",1,])
		trips.upper.ylim <- ylim.max(real.trips[[rgn]]["Estimated richness","upper CI",i,])
		trips.ylim <- ylim.max(real.trips[[rgn]]["Estimated richness","MLE",i,])

		if (nulldist == TRUE && plot.nulls == TRUE) {trips.null.ylim <- ylim.max(c(null.trips[[1]]["Estimated richness","MLE",i,,], trips.ylim))}
		if (trips.upper.ylim > raw.ylim * 3) {raw.ylim <- median(c(trips.ylim,trips.upper.ylim), na.rm = T)} else {raw.ylim <- trips.upper.ylim}
		plot(1, type = 'n', xlim = xlim, ylim = c(0,raw.ylim), bty = "l", ylab = "", xlab = "", main = "", log = logged, cex.lab = cex.lab)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		if (i == focal_intervals[1]) {mtext("Species Richness", side = 2, line = 2, cex = cex.lab)}
			points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#TRiPS
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (j in 1:nshuff) {
				lines(null.trips[[1]]["Estimated richness","MLE",i,j,] ~ numoccs[[rgn]][[i]], col = adjustcolor("blue", alpha.f = 0.1), lwd = null.lwd)
			}
		}
		lower_ci <- real.trips[[rgn]]["Estimated richness","lower CI",i,]
		upper_ci <- real.trips[[rgn]]["Estimated richness","upper CI",i,]
		noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
		lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
		polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor("blue", alpha.f = 0.1), border = TRUE, lty = 2, lwd = 0.5)
		points(numoccs[[rgn]][[i]], real.trips[[rgn]]["Estimated richness","MLE",i,], type = "l", lty = 2, pch = 21, cex = 0.5, col = "blue", lwd = colour.lwd)

		#raw legend ####
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = c("red","blue"), legend = c("Raw","TRiPS"), lty = c(1,2), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (plot_samplingrate == T) {
		#TRiPS sampling rate ####
		ylim <- ylim.max(real.trips[[rgn]]["Sampling rate","upper CI",i,])
		xlim <- c(1,max(numoccs[[rgn]][[i]]))
		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		if (is.infinite(ylim) == FALSE) {
			plot(1, type = 'n', xlim = xlim, ylim = c(0,ylim), bty = "l", ylab = "", xlab = "", main = "", cex.lab = cex.lab)
			abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
			lower_ci <- real.trips[[rgn]]["Sampling rate","lower CI",i,]
			upper_ci <- real.trips[[rgn]]["Sampling rate","upper CI",i,]
			noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
			lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
			polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor("red", alpha.f = 0.1), border = TRUE, lty = 2, lwd = 0.5)
			points(numoccs[[rgn]][[i]], real.trips[[rgn]]["Sampling rate","MLE",i,], type = "l", lty = 1, cex = 0.5, col = "red", lwd = black.lwd)
		} else {plot(1, type = 'n', xlim = xlim, ylim = c(0,1), yaxt = "n", xaxt = "n", bty = "l", ylab = "", xlab = "Occurrences", main = "", cex.lab = cex.lab); abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)}
		if (i == focal_intervals[1]) {mtext("Sampling Rate", side = 2, line = 2, cex = cex.lab)}
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == rev(focal_intervals)[1]) {par(xpd = NA); legend("right", inset = legend_inset, col = "red", legend = "TRiPS\nSampling\nRate\n", lty = 1, bty = "n", lwd = colour.lwd, cex = legend.cex); par(xpd = FALSE)}
	}

	if (plot_completeness == T) {
		#  TRiPS sampling probability ####
		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		if (all(is.na(real.trips[[rgn]]["Sampling probability","MLE",i,])) == FALSE) {
			plot(1, type = 'n', xlim = xlim, ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "l")
			abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
			lower_ci <- real.trips[[rgn]]["Sampling probability","lower CI",i,]
			upper_ci <- real.trips[[rgn]]["Sampling probability","upper CI",i,]
			noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
			lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
			polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor("black", alpha.f = 0.1), border = TRUE, lty = 2, lwd = 0.5)
			points(numoccs[[rgn]][[i]], real.trips[[rgn]]["Sampling probability","MLE",i,], type = "l", lty = 1, cex = 0.5, col = "black", lwd = colour.lwd)
		} else {plot(1, type = 'n', xlim = xlim, ylim = c(0,1), yaxt = "n", xaxt = "n", bty = "l", ylab = "", xlab = "Occurrences", main = ""); abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)}
		par(new = TRUE)
		plot(1, type = 'n', xlim = xlim, ylim = c(0,1), xlab = "", ylab = "", bty = "l", cex.lab = cex.lab)
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		if (i == focal_intervals[1]) {mtext("Sample Completeness", side = 2, line = 2, cex = cex.lab)}
			points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Good.s.u",1,], type = "l", lty = 1, cex = 0.5, col = "blue", lwd = black.lwd)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == rev(focal_intervals)[1]) {par(xpd = NA); legend("right", inset = legend_inset, col = c("black","blue"), legend = c("Sampling\nProbability\n","Good's u"), lty = 1, bty = "n", lwd = colour.lwd, cex = legend.cex); par(xpd = FALSE)}
	}

	if (plot_ntons == T) {
		# total singleton and multiton counts ####
		xlim <- c(1,max(numoccs[[rgn]][[i]]))
		ylim <- c(1,max(c(singleton_count[[rgn]][[i]], multiton_count[[rgn]][[i]], na.rm = T)))

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(singleton_count[[rgn]][[i]] ~ numoccs[[rgn]][[i]], type = 'l', lwd = colour.lwd, col = "red", ylim = ylim, xlim = xlim, xlab = "Species", ylab = "", bty = "l", main = "")
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		lines(multiton_count[[rgn]][[i]] ~ numoccs[[rgn]][[i]], type = 'l', lwd = colour.lwd, col = "blue")
		if (i == focal_intervals[1]) {mtext("Species", side = 2, line = 2, cex = cex.lab)}
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = c("red","blue"), legend = c("singletons","multitons"), bty = "n", lty = 1, lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (run_PGS_calcs == T && plot_pgs == T) {
		#palaeogeographic spread ####
		length.ylim <- c(0,ylim.max(pgs[[rgn]][i,c("mst_dist","max_gcd","mean_pairwise_GCD","median_pairwise_GCD"),]))
		pgs.xlim <- c(0, max(numoccs[[rgn]][[i]], na.rm = TRUE))
		occupancy.ylim <- c(0,ylim.max(pgs[[rgn]][i,"grid_cell_occupancy",]))
		STD.ylim <- c(0,100)

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		if (i == focal_intervals[1]) {}
		plot(1, type = "n", ylab = "", xlim = pgs.xlim, ylim = length.ylim, xlab = "Occurrences", main = "", bty = "l")
		if (i == focal_intervals[1]) {mtext("Distance (km)", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		points(pgs[[rgn]][i,"max_gcd",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = "black", lwd = black.lwd); points(pgs[[rgn]][i,"max_gcd",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = pgs.cols["max_gcd"], lwd = colour.lwd)
		points(pgs[[rgn]][i,"mst_dist",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = "black", lwd = black.lwd); points(pgs[[rgn]][i,"mst_dist",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = pgs.cols["mst_dist"], lwd = colour.lwd)
		par(new = T)
		plot(1, type = "n", yaxt = "n", ylab = "", xaxt = "n", xlim = pgs.xlim, ylim = occupancy.ylim, xlab = "", main = "", bty = "l", axes = FALSE)
		points(pgs[[rgn]][i,"grid_cell_occupancy",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = "black", lwd = black.lwd); points(pgs[[rgn]][i,"grid_cell_occupancy",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = pgs.cols["grid_cell_occupancy"], lwd = colour.lwd)
		if (i == rev(focal_intervals)[1]) {
			axis(4, line = 0)
			mtext("Occupied Grid Cells", side = 4, line = 1.5, cex = 0.7)
		}
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			pgs_legend_offset <- 1.3
			legend("right", inset = legend_inset * pgs_legend_offset, col = "black", legend = c("MST Length","GCO","Max GCD"), text.col = "white", bty = "n", lty = 1, lwd = black.lwd, cex = legend.cex)
			legend("right", inset = legend_inset * pgs_legend_offset, col = pgs.cols[-which(names(pgs.cols) %in% c("convex_hull","STD"))], legend = c("MST Length","GCO","Max GCD"), bty = "n", lty = 1, lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
		# par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * 3), par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3]), paste("(",letters[3 + which(focal_intervals == i)],")", sep = ""), cex = 1.5, font = 2, family = 'sans'); par(xpd = F)
	}

	par(mfg = c(1,which(focal_intervals == i)))
	title(intervals$bin[i], line = 1.1)
	axis(side = 3, at = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], labels =  rev(decades), col = "white", padj = 1.5, las = 1, cex.axis = 0.75)

}

graphics.off()


#### Figure X - raw and standardised species-accumulation curves ####
quorum <- "Quorum 0.4"
if (clade == "tetrapoda" & use.stages == T) {
	pub_ints <- 19:52
} else if (clade == "tetrapoda" & use.stages == FALSE) {
	pub_ints <- 8:29
	# pub_ints <- 1:22
} else {
	pub_ints <- 1:nrow(intervals)
}
log_xy <- c("","xy")
xlim_occs <- c(1, max(sapply(numoccs[[rgn]][pub_ints], max)))
pdf(paste(folder.name, "/plots/Figure <", clade, "-", "collector-curves-all-intervals> for ", names(regions)[rgn], ", stages = ", use.stages, ", ", quorum, ".pdf", sep = ""), width = 16, height = 8)
par(mfrow = c(2,4), oma = c(0,0,0,8.5), pty = "s")
for (l in 1:2) {
	label.height <- par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3])
	#raw diversity
	plot(1, type = 'n', xlim = xlim_occs, ylim = c(1,max(real.sqs[[rgn]][pub_ints,"Raw.species",1,], na.rm = T)), bty = "l", ylab = "Raw diversity", xlab = "Occurrences", log = log_xy[l])
	if (nulldist == TRUE) {
	}
	for (i in pub_ints) {
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",quorum,], type = "l", pch = 16, cex = 0.5, col = "black", lwd = black.lwd)
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",quorum,], type = "l", pch = 16, cex = 0.5, col = interval.cols[i], lwd = colour.lwd)
	}
	if (l == 1) {par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * 2.5), par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3]), "(a)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}
	if (l == 2) {par(xpd = NA); text(par("usr")[2] * 0.2 + (par("usr")[1] * 2.5), label.height, "(e)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}

	#SQS
	#Quorum 0.4
	quorum <- "Quorum 0.4"
	ylim <- ylim.max(real.sqs[[rgn]][pub_ints,"Subsampled.diversity",quorum,])
	plot(1, type = 'n', xlim = xlim_occs, ylim = c(1,ylim), bty = "l", ylab = paste("SQS Subsampled Richness for", quorum, sep = " "), xlab = "Occurrences", log = log_xy[l])
	if (nulldist == TRUE) {
	}
	for (i in pub_ints) {
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Subsampled.diversity",quorum,], type = "l", pch = 16, cex = 0.5, col = "black", lwd = black.lwd)
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Subsampled.diversity",quorum,], type = "l", pch = 16, cex = 0.5, col = interval.cols[i], lwd = colour.lwd)
	}
	if (l == 1) {par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * 2.5), par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3]), "(b)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}
	if (l == 2) {par(xpd = NA); text(par("usr")[2] * 0.2 + (par("usr")[1] * 2.5), label.height, "(f)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}

	#Quorum 0.5
	quorum <- "Quorum 0.5"
	ylim <- ylim.max(real.sqs[[rgn]][pub_ints,"Subsampled.diversity",quorum,])
	plot(1, type = 'n', xlim = xlim_occs, ylim = c(1,ylim), bty = "l", ylab = paste("SQS Subsampled Richness for", quorum, sep = " "), xlab = "Occurrences", log = log_xy[l])
	if (nulldist == TRUE) {
	}
	for (i in pub_ints) {
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Subsampled.diversity",quorum,], type = "l", pch = 16, cex = 0.5, col = "black", lwd = black.lwd)
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Subsampled.diversity",quorum,], type = "l", pch = 16, cex = 0.5, col = interval.cols[i], lwd = colour.lwd)
	}
	if (l == 1) {par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * 2.5), par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3]), "(c)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}
	if (l == 2) {par(xpd = NA); text(par("usr")[2] * 0.2 + (par("usr")[1] * 2.5), label.height, "(g)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}

	#TRiPS
	par(mar = c(5.1, 4.1, 4.1, 2.1))
	plot(1, type = 'n', xlim = xlim_occs, ylim = c(1,max(real.trips[[rgn]]["Estimated richness","MLE",pub_ints,], na.rm = T)), bty = "l", ylab = "TRiPS Estimated Richness", xlab = "Occurrences", log = log_xy[l])
	for (i in pub_ints) {
		points(numoccs[[rgn]][[i]], real.trips[[rgn]]["Estimated richness","MLE",i,], type = "l", pch = 16, cex = 0.5, col = "black", lwd = black.lwd)
		points(numoccs[[rgn]][[i]], real.trips[[rgn]]["Estimated richness","MLE",i,], type = "l", pch = 16, cex = 0.5, col = interval.cols[i], lwd = colour.lwd)
	}
	if (l == 1) {par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * 2.5), par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3]), "(d)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}
	if (l == 2) {par(xpd = NA); text(par("usr")[2] * 0.2 + (par("usr")[1] * 2.5), label.height, "(h)", cex = 1.5, font = 2, family = 'sans'); par(xpd = F)}

	if (l == 1)
	{
		par(xpd = NA)
		legend("topright", inset = c(-0.4,0.4), text.col = "white", col = "black", legend = intervals$bin[pub_ints], lty = 1, bty = "n", lwd = black.lwd, cex = 0.9)
		legend("topright", inset = c(-0.4,0.4), col = interval.cols[pub_ints], legend = intervals$bin[pub_ints], lty = 1, bty = "n", lwd = colour.lwd, cex = 0.9)
		par(xpd = FALSE)
	}
}
graphics.off()




#### Figure X (singletons/multitons and palaeogeographic spread for focal intervals) ####
legend.cex <- 0.55
pdf(paste(folder.name, "/plots/Figure <", clade, "-singletons-multitons-and-spread> for ", names(regions)[rgn], " ", clade, ", stages = ", use.stages, ".pdf", sep = ""), width = 14, height = 8)
par(mfrow = c(2,3), oma = c(0,0,0,0), mar = c(4.1, 5, 4.1, 6), pty = "s")

for (i in focal_intervals) {
	xlim <- c(1,max(numoccs[[rgn]][[i]]))
	ylim <- c(1,max(c(singleton_count[[rgn]][[i]], multiton_count[[rgn]][[i]], na.rm = T)))

	plot(singleton_count[[rgn]][[i]] ~ numoccs[[rgn]][[i]], type = 'l', lwd = 2, col = "red", ylim = ylim, xlab = "Total Occurrences", ylab = "Singletons/Multitons", bty = "l", main = intervals$bin[i])
	abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
	axis(side = 3, at = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], labels =  rev(decades), col = "white", padj = 1, las = 1)
	lines(multiton_count[[rgn]][[i]] ~ numoccs[[rgn]][[i]], type = 'l', lwd = 2, col = "blue")
	par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * 3), par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3]), paste("(",letters[which(focal_intervals == i)],")", sep = ""), cex = 1.5, font = 2, family = 'sans'); par(xpd = F)
	if (i == focal_intervals[1]) {legend("topleft", legend = c("total singletons","total multitons"), lty = 1, lwd = 2, col = c("red","blue"), bty = "n", cex = 0.7)}
}
for (i in focal_intervals) {
	length.ylim <- c(0,ylim.max(pgs[[rgn]][i,c("mst_dist","max_gcd","mean_pairwise_GCD","median_pairwise_GCD"),]))
	pgs.xlim <- c(0, max(numoccs[[rgn]][[i]], na.rm = TRUE))
	occupancy.ylim <- c(0,ylim.max(pgs[[rgn]][i,"grid_cell_occupancy",]))
	STD.ylim <- c(0,100)

	plot(1, type = "n", ylab = "Distance (km)", xlim = pgs.xlim, ylim = length.ylim, xlab = "Occurrences", main = intervals$bin[i], bty = "l")
	axis(1, line = 0)
	abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
	axis(side = 3, at = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], labels =  rev(decades), col = "white", padj = 1, las = 1)
	points(pgs[[rgn]][i,"max_gcd",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = "black", lwd = 4); points(pgs[[rgn]][i,"max_gcd",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = pgs.cols["max_gcd"], lwd = 3.5)
	points(pgs[[rgn]][i,"mst_dist",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = "black", lwd = 4); points(pgs[[rgn]][i,"mst_dist",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = pgs.cols["mst_dist"], lwd = 3.5)
	points(pgs[[rgn]][i,"median_pairwise_GCD",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = "black", lwd = 3.5); points(pgs[[rgn]][i,"median_pairwise_GCD",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = pgs.cols["median_pairwise_GCD"], lwd = 3)
	par(new = T)
	plot(1, type = "n", yaxt = "n", ylab = "", xaxt = "n", xlim = pgs.xlim, ylim = occupancy.ylim, xlab = "", main = intervals$bin[i], bty = "l", axes = FALSE)
	points(pgs[[rgn]][i,"grid_cell_occupancy",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = "black", lwd = 3.5); points(pgs[[rgn]][i,"grid_cell_occupancy",] ~ numoccs[[rgn]][[i]], type = "l", pch = 16, col = pgs.cols["grid_cell_occupancy"], lwd = 3)
	axis(4, line = 0.5)
	mtext("Occupied cells (n)", side = 4, line = 2.5, cex = 0.7)
	par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * 3), par("usr")[4] + 0.05 * (par("usr")[4] - par("usr")[3]), paste("(",letters[3 + which(focal_intervals == i)],")", sep = ""), cex = 1.5, font = 2, family = 'sans'); par(xpd = F)
	if (i == focal_intervals[1]) {legend("topleft", legend = gsub(x = names(pgs.cols[-which(names(pgs.cols) %in% c("convex_hull","STD","mean_pairwise_GCD"))]), "_", replacement = " "), lty = 1, lwd = 3, col = pgs.cols[-which(names(pgs.cols) %in% c("convex_hull","STD","mean_pairwise_GCD"))], bty = "n", cex = 0.7)}
}
graphics.off()

#### Figure with null distributions ####
plot_SQS1 <- T
plot_SQS2 <- F
plot_iNEXT <- F
plot_chao <- T
plot_lambda5 <- T
plot_TRiPS <- T
plot_multitonsubsampling <- F
plot_CR <- F
fig1rows <- sum(plot_SQS1,plot_SQS2,plot_iNEXT,plot_CR,plot_multitonsubsampling,plot_chao,plot_lambda5,plot_TRiPS)
fig1columns <- length(focal_intervals)

ci.alpha <- 0.1
plot.nulls <- T
plot.CI.border <- F
null.lwd <- 3
null.alpha <- 0.1

graphics.off()

pdf(paste(folder.name, "/plots/Figure <", clade, "-diversity-plots> use.stages=", use.stages, " with null distributions.pdf", sep = ""), width = fig1columns*3.5, height = fig1rows*3)
par(mfrow = c(fig1rows, fig1columns))
par(cex = 0.7)
par(mar = c(0, 3, 2, 0))
par(oma = c(5, 1, 0.5, 11))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
cex.lab <- 0.7
vlaboffset <- 0.09
hlaboffset <- 2.5
legend_inset <- c(-0.5,0)
legend.cex <- 1
colour.lwd <- 1.5; black.lwd <- colour.lwd * 1.2
border <- 1 #1 = yes, 2 = no
logged <- "y"

myletters <- myLetters(length(focal_intervals) * fig1rows)
# myletters <- NA

for (i in focal_intervals) {

	xlim <- c(1,max(numoccs[[rgn]][[i]]))
	null.x <- rev(round(seq(from = 1, to = max(numoccs[[rgn]][[i]]), length.out = length(times))))

	p <- 0

	if (run_SQS_Perl_1 == T && plot_SQS1 == T) {
		#plot subsampled diversity (SQS with all the trimmings) ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), function(x) {
			tmp <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[x], flavour == "SQS1")
			tmp$Subsampled.diversity})))

		raw.ylim <- ylim.max(real.sqs[[rgn]][i,"Raw.species",1,])
		sqs.ylim <- ylim.max(c(sqs.ylim, raw.ylim))

		if (nulldist == TRUE) {
			sqs.null.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), function(q) {
				null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,]
			})))
			sqs.ylim <- ylim.max(c(sqs.ylim, sqs.null.ylim))
		}

		p <- p + 1

		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(1,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, log = logged)
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(SQS Perl Script V1)", side = 2, line = 2, cex = cex.lab)}
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		if (nulldist == TRUE) {
			for (q in 1:length(quorum.levels)) {
				lower_ci <- sapply(1:length(times), function(j) min(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,j,]))
				upper_ci <- sapply(1:length(times), function(j) max(null.sqs.1[[rgn]][i,"Subsampled.diversity",q,j,]))
				noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
				lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
				polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[q], alpha.f = null.alpha), border = FALSE, lty = 2, lwd = 0.5)

				for (r in 1:nshuff) {
					lines((null.sqs.1[[rgn]][i,"Subsampled.diversity",q,,r]) ~ null.x, col = adjustcolor(quorum.cols[q], alpha.f = null.alpha), lwd = null.lwd)
				}
			}
		}

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#sqs richness
		for (j in 1:length(quorum.levels)) {
			X <- filter(df, region == names(regions)[rgn], bin == intervals$bin[i])
			Y <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[j], flavour == "SQS1")

			for (k in border:2) {points(X$numoccs, Y$Subsampled.diversity, type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		#legend
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white", title = "Quorum Level")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex, title = "Quorum Level")
			par(xpd = FALSE)
		}
	}

	if (run_SQS_Perl_2 == T && plot_SQS2 == T) {
		#plot subsampled diversity (plain-vanilla SQS) ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), FUN = function(x) real.sqs.2[[rgn]][i,"Subsampled.diversity",x,])))
		raw.ylim <- ylim.max(real.sqs[[rgn]][i,"Raw.species",1,])
		sqs.ylim <- ylim.max(c(sqs.ylim, raw.ylim))

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(1,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, main = "", log = logged)
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(SQS Perl Script V2)", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (nulldist == TRUE) {
			for (q in 1:length(quorum.levels)) {
				lower_ci <- sapply(1:length(times), function(j) min(null.sqs.2[[rgn]][i,"Subsampled.diversity",q,j,]))
				upper_ci <- sapply(1:length(times), function(j) max(null.sqs.2[[rgn]][i,"Subsampled.diversity",q,j,]))
				noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
				lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
				polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[q], alpha.f = null.alpha), border = FALSE, lty = 2, lwd = 0.5)

				for (r in 1:nshuff) {
					lines((null.sqs.2[[rgn]][i,"Subsampled.diversity",q,,r]) ~ null.x, col = adjustcolor(quorum.cols[q], alpha.f = null.alpha), lwd = null.lwd)
				}
			}
		}

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#sqs richness
		for (j in 1:length(quorum.levels)) {
			X <- filter(df, region == names(regions)[rgn], bin == intervals$bin[i])
			Y <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[j], flavour == "SQS2")

			for (k in border:2) {points(X$numoccs, Y$Subsampled.diversity, type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (run_iNEXT == T && plot_iNEXT == T) {
		#plot subsampled diversity (iNEXT) ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quorum.levels), FUN = function(x) filter(iNEXT_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[x], order == 0, method != "extrapolated")$qD.UCL)))
		raw.ylim <- ylim.max(real.sqs[[rgn]][i,"Raw.species",1,])
		sqs.ylim <- ylim.max(c(sqs.ylim, raw.ylim))

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(1,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, main = "", log = logged)
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(iNEXT CBR)", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#sqs richness
		for (j in 1:length(quorum.levels)) {
			Y2 <- filter(iNEXT_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[j], order == 0, method != "extrapolated")
			X <- left_join(Y2, df)

			lower_ci <- Y2$qD.LCL
			upper_ci <- Y2$qD.UCL
			lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
			polygon(c(X$numoccs[which(!is.na(lower_ci))], rev(X$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[j], alpha.f = null.alpha), border = plot.CI.border, lty = 2, lwd = 0.5)
			for (k in border:2) {points(X$numoccs, Y2$qD, type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (run_CR == T && plot_CR == T) {
		#plot iNEXT CR ####
		#ylims for per-interval maxima
		sqs.ylim <- ylim.max(unlist(lapply(1:length(quota.levels), function(x) {
			tmp <- filter(CR_df, region == names(regions)[rgn], bin == intervals$bin[i], quota == names(quota.levels)[x])
			tmp$qD.UCL}
		)))
		raw.ylim <- ylim.max(real.sqs[[rgn]][i,"Raw.species",1,])
		sqs.ylim <- ylim.max(c(sqs.ylim, raw.ylim))

		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(1,sqs.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, main = "", log = logged)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == focal_intervals[1]) {mtext("Rescaled Richness\n(iNEXT CR))", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)
		X <- filter(df, region == names(regions)[rgn], bin == intervals$bin[i])

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#sqs richness

		for (j in 1:length(quota.levels)) {
			Y2 <- filter(CR_df, region == names(regions)[rgn], bin == intervals$bin[i], quota == names(quota.levels)[j], order == 0)
			lower_ci <- Y2$qD.LCL
			upper_ci <- Y2$qD.UCL
			lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
			polygon(c(X$numoccs[which(!is.na(lower_ci))], rev(X$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor(quorum.cols[j], alpha.f = null.alpha), border = plot.CI.border, lty = 2, lwd = 0.5)
			for (k in border:2) {points(X$numoccs, Y2$qD, type = "l", lty = 1, cex = 0.5, col = c("black",quorum.cols[j])[k], lwd = c(black.lwd,colour.lwd)[k])}
		}

		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = rep("black", length(quorum.levels)), lty = 1, legend = rev(paste("q =",quorum.levels)), bty = "n", lwd = rep(black.lwd*1.1, length(quorum.levels)), cex = legend.cex, text.col = "white")
			legend("right", inset = legend_inset, col = rev(quorum.cols), legend = rev(paste("q =",quorum.levels)), lty = rep(1, length(quorum.levels)), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (plot_chao == T) {
		Y1 <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[1], flavour == "SQS2")
		Y3 <- filter(chao_df, region == names(regions)[rgn], bin == intervals$bin[i])
		X1 <- left_join(Y1, df)
		X3 <- left_join(Y3, df)
		ylim <- ylim.max(c(Y1$Chao.2, Y3$chao.UCL))
		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(1, ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, log = logged)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == focal_intervals[1]) {mtext("Species Richness", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)

		#null
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (j in 1:nshuff) {
				lines(null.chao[[rgn]][,i,,j] ~ null.x, col = adjustcolor("blue", alpha.f = null.alpha), lwd = null.lwd)
			}
		}

		#Chao2 ####
		points(X1$numoccs, Y1$Chao.2, type = "l", lty = 2, pch = 21, cex = 0.5, col = "blue", lwd = 1.5)

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#chao
		lower_ci <- Y3$chao.LCL
		upper_ci <- Y3$chao.UCL
		lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
		polygon(c(X3$numoccs[which(!is.na(lower_ci))], rev(X3$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor("blue", alpha.f = ci.alpha), border = plot.CI.border, lty = 2, lwd = border.lwd)
		points(X3$numoccs, Y3$chao, type = "l", lty = 2, pch = 21, cex = 0.5, col = "blue", lwd = colour.lwd)

		#legend
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = c("red","black"), legend = c("Raw","Chao1"), lty = c(1,2), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (plot_lambda5 == T) {
		Y1 <- filter(SQS_df, region == names(regions)[rgn], bin == intervals$bin[i], quorum == names(quorum.levels)[1], flavour == "SQS2")
		Y2 <- filter(lambda5_df, region == names(regions)[rgn], bin == intervals$bin[i])
		X1 <- left_join(Y1, df)
		X2 <- left_join(Y2, df)
		ylim <- ylim.max(Y2$lambda5.UCL)
		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(1, ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, log = logged)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == focal_intervals[1]) {mtext("Species Richness", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)

		#null
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (j in 1:nshuff) {
				lines(null.lambda5[[rgn]][,i,,j] ~ null.x, col = adjustcolor("blue", alpha.f = null.alpha), lwd = null.lwd)
			}
		}

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#lambda5
		lower_ci <- Y2$lambda5.LCL
		upper_ci <- Y2$lambda5.UCL
		lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
		polygon(c(X2$numoccs[which(!is.na(lower_ci))], rev(X2$numoccs[which(!is.na(upper_ci))])), c(lower_ci,rev(upper_ci)), col = adjustcolor("blue", alpha.f = ci.alpha), border = plot.CI.border, lty = 2, lwd = border.lwd)
		points(X2$numoccs, Y2$lambda5, type = "l", lty = 2, pch = 21, cex = 0.5, col = "blue", lwd = colour.lwd)

		#legend
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = c("red","black"), legend = c("Raw","lambda-5"), lty = c(1,2,2,2), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}

	if (plot_TRiPS == T) {
		trips.upper.ylim <- ylim.max(real.trips[[rgn]]["Estimated richness","upper CI",i,])
		p <- p + 1
		par(mfg = c(p,which(focal_intervals == i)))
		plot(1, type = 'n', xlim = xlim, ylim = c(1, chao2.ylim), bty = "l", xlab = "", ylab = "", cex.lab = cex.lab, log = logged)
		par(xpd = NA); text(par("usr")[1] + (par("usr")[1] * hlaboffset), par("usr")[4] + vlaboffset * (par("usr")[4] - par("usr")[3]), paste("(", myletters[(match(row(intervals)[i], focal_intervals)) + ((p - 1) * fig1columns)], ")", sep = ""), cex = 1.05, font = 2, family = 'sans'); par(xpd = FALSE)
		if (i == focal_intervals[1]) {mtext("Species Richness", side = 2, line = 2, cex = cex.lab)}
		abline(v = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], lty = 2, col = "lightgrey", lwd = 0.5)

		#TRiPS
		if (nulldist == TRUE && plot.nulls == TRUE) {
			for (j in 1:nshuff) {
				lines(null.trips[[rgn]]["Estimated richness","MLE",i,,j] ~ null.x, col = adjustcolor("blue", alpha.f = null.alpha), lwd = null.lwd)
			}
		}

		lower_ci <- real.trips[[rgn]]["Estimated richness","lower CI",i,]
		upper_ci <- real.trips[[rgn]]["Estimated richness","upper CI",i,]
		noccs <- numoccs[[rgn]][[i]][!is.na(upper_ci)]
		lower_ci <- lower_ci[!is.na(lower_ci)]; upper_ci <- upper_ci[!is.na(upper_ci)]
		polygon(c(noccs, rev(noccs)), c(lower_ci,rev(upper_ci)), col = adjustcolor("blue", alpha.f = ci.alpha), border = plot.CI.border, lty = 2, lwd = border.lwd)
		points(numoccs[[rgn]][[i]], real.trips[[rgn]]["Estimated richness","MLE",i,], type = "l", lty = 2, pch = 21, cex = 0.5, col = "blue", lwd = colour.lwd)

		#raw richness
		points(numoccs[[rgn]][[i]], real.sqs[[rgn]][i,"Raw.species",1,], type = "l", lty = 1, col = "red", lwd = colour.lwd, xlim = xlim)

		#legend
		if (i == rev(focal_intervals)[1]) {
			par(xpd = NA)
			legend("right", inset = legend_inset, col = c("red","black"), legend = c("Raw","TRiPS"), lty = c(1,2), bty = "n", lwd = colour.lwd, cex = legend.cex)
			par(xpd = FALSE)
		}
	}



	par(mfg = c(1,which(focal_intervals == i)))
	title(intervals$bin[i], line = 1.1)
	axis(side = 3, at = numoccs[[rgn]][[i]][names(numoccs[[rgn]][[i]]) %in% names(decades)], labels =  rev(decades), col = "white", padj = 1.5, las = 1, cex.axis = 0.75)

}

graphics.off()


#### Figure S<region-clade-diversity-curves>
pdf(paste(folder.name, "/plots/Fig. S<",names(regions)[rgn], "-", clade, "-diversity-curves>, stages = ", use.stages, ".pdf", sep = ""), width = 6, height = 12)
#merge data objects
q <- "Quorum 0.4"
tmp_SQS_1_df <- SQS_1_df; tmp_SQS_1_df[which(tmp_SQS_1_df$Raw.species < 2), "Raw.species"] <- NA
tmp_SQS1 <- tmp_SQS_1_df %>% filter(quorum == q & year == 2016) %>% select(-method, -midpoint) %>% rename(`Face-value Richness` = Raw.species, `SQS (Perl script v1)` = Subsampled.diversity, midpoint = Midpoint) %>% gather("method", "qD", c(`Face-value Richness`,`SQS (Perl script v1)`)); tmp_SQS1$qD.LCL <- tmp_SQS1$qD.UCL <- NA
tmp_SQS2 <- SQS_2_df %>% filter(quorum == q & year == 2016) %>% rename(`SQS (Perl script v2)` = Subsampled.diversity) %>% select(-method) %>% gather("method", "qD", `SQS (Perl script v2)`); tmp_SQS2$qD.LCL <- tmp_SQS2$qD.UCL <- NA
tmp_iNEXT <- iNEXT_df %>% filter(quorum == q & year == 2016 & method != "extrapolated" & order == 0); tmp_iNEXT$method <- "SQS (iNEXT)"
tmp_CR <- CR_df %>% filter(quota == "Quota 50" & year == 2016 & method != "extrapolated" & order == 0); tmp_CR$method <- "Classical Rarefaction"
tmp_TRiPS <- TRiPS_df %>% filter(year == 2016); tmp_TRiPS$method <- "TRiPS"
tmp_chao <- chao_df %>% filter(year == 2016) %>% rename(qD = chao, qD.LCL = chao.LCL, qD.UCL = chao.UCL); tmp_chao$method <- "Chao2"; tmp_chao[which(tmp_chao$qD < 2), "qD"] <- tmp_chao[which(tmp_chao$qD < 2), "qD.LCL"] <- tmp_chao[which(tmp_chao$qD < 2), "qD.UCL"] <- NA
tmp_lambda5 <- lambda5_df %>% filter(year == 2016) %>% rename(qD = lambda5, qD.LCL = lambda5.LCL, qD.UCL = lambda5.UCL); tmp_lambda5$method <- "lambda-5"; tmp_lambda5[which(tmp_lambda5$qD < 2), "qD.LCL"] <- tmp_lambda5[which(tmp_lambda5$qD < 2), "qD.UCL"] <- NA; tmp_lambda5[which(tmp_lambda5$qD < 2), "qD"] <- NA

all_df <- full_join(tmp_lambda5,
		    full_join(tmp_chao,
		    	  full_join(tmp_TRiPS,
		    	  	  full_join(tmp_SQS1,
		    	  	  	  full_join(tmp_SQS2,
		    	  	  	  	  full_join(tmp_iNEXT, tmp_CR))))))

all_df$method = factor(all_df$method, levels = c("Face-value Richness","TRiPS","Chao2","lambda-5","SQS (Perl script v1)","SQS (Perl script v2)","SQS (iNEXT)","Classical Rarefaction")) #reorder method factor

all_df %>%
	filter(!is.na(qD) & qD > 2) %>%
	filter(midpoint >= 51.9 & midpoint <= 251.685) %>% #for NAm tetrapods
	ggplot(aes(x = midpoint, y = qD, ymin = qD.LCL, ymax = qD.UCL)) +
	geom_line(lty = 2, col = "darkgrey") +
	geom_point() +
	theme_bw() +
	geom_vline(xintercept = period.boundaries[which(period.boundaries < 250)], lty = 2, col = "grey") +
	geom_vline(xintercept = 66) +
	scale_x_reverse() +
	labs(colour = "Richness Standardisation Method", pch = "Richness Standardisation Method", x = "Age (myr)", y = "Estimated Species Richness") +
	scale_colour_brewer(palette = "Set1") +
	geom_errorbar() +
	facet_grid(method~., scales = "free_y") +
	scale_y_log10()
graphics.off()


#### COMPARISON OF SQS RICHNESS ESTIMATES FROM iNEXT AND THE SQS PERL SCRIPT VERSION 2
pdf("./simulation-figures/Fig. S<see-they-look-the-same>.pdf", width = 5, height = 5)
tmp1 <- iNEXT_df %>% filter(method == "interpolated", quorum == "Quorum 0.5", region == "global", order == 0)
tmp2 <- SQS_2_df %>% filter(quorum == "Quorum 0.5", region == "global")
tmp2 <- tmp2 %>% rename(qD = Subsampled.diversity)
tmp3 <- full_join(select(tmp1, bin, year, flavour, qD), select(tmp2, bin, year, flavour, qD))
tmp3 %>% spread(flavour, qD) %>% ggplot(aes(x = iNEXT, y = SQS2)) + geom_point() + theme_bw() + geom_abline(slope = 1) + coord_equal() + labs(x = "Coverage-standardised Richness (iNEXT)", y = "Coverage-standardised Richness SQS Perl Script V2")
graphics.off()

#### COMPARISON OF SQS RICHNESS ESTIMATES FROM SQS PERL SCRIPT VERSIONS 1 & 2
tmp1 <- SQS_1_df %>% filter(quorum == "Quorum 0.5", region == "global")
tmp2 <- SQS_2_df %>% filter(quorum == "Quorum 0.5", region == "global")
tmp1 <- tmp1 %>% rename(qD = Subsampled.diversity)
tmp2 <- tmp2 %>% rename(qD = Subsampled.diversity)
tmp3 <- full_join(select(tmp1, bin, year, flavour, qD), select(tmp2, bin, year, flavour, qD))
pdf("./simulation-figures/Fig. S<plain-vanilla-sqs-vs-collsperref-sqs>.pdf", width = 5, height = 5)
tmp3 %>% spread(flavour, qD) %>% ggplot(aes(x = SQS1, y = SQS2)) + geom_point() + theme_bw() + geom_abline(slope = 1) + coord_equal() + labs(x = "Coverage-standardised Richness SQS Perl Script V1", y = "Coverage-standardised Richness SQS Perl Script V2")
graphics.off()
