##############################################################################
##### flatten results contained in recursive lists into tidy data frames #####
##############################################################################

df <- list(lapply(multiton_count, function(x) do.call("rbind", x)) %>% melt() %>% rename(bin = Var1, year = Var2, multiton_count = value, region = L1),
	   lapply(singleton_count, function(x) do.call("rbind", x)) %>% melt() %>% rename(bin = Var1, year = Var2, singleton_count = value, region = L1),
	   lapply(newoccs, function(x) do.call("rbind", x)) %>% melt() %>% rename(bin = Var1, year = Var2, newoccs = value, region = L1),
	   lapply(newtax, function(x) do.call("rbind", x)) %>% melt() %>% rename(bin = Var1, year = Var2, newtax = value, region = L1),
	   lapply(numoccs, function(x) do.call("rbind", x)) %>% melt() %>% rename(bin = Var1, year = Var2, numoccs = value, region = L1),
	   goodsu %>% melt() %>% rename(bin = Var1, year = Var2, goodsu = value, region = L1), #lapply(goodsu, function(x) do.call("rbind", x))
	   melt(pgs) %>% spread(Var2, value) %>% rename(bin = Var1, year = Var3, region = L1)
) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2), .)

df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = df$year))
df <- tbl_df(df)

#### flatten SQS Perl script with collsperref into tidy data frame ####
if (run_SQS_Perl_1 == TRUE) {
	SQS_1_df <- melt(real.sqs.1)
	SQS_1_df <- spread(SQS_1_df, Var2, value)
	SQS_1_df <- tbl_df(rename(SQS_1_df, bin = Var1, quorum = Var3, year = Var4, region = L1))
	SQS_1_df <- tbl_df(cbind(SQS_1_df, select(intervals[match(SQS_1_df$bin, intervals$bin), ], LAD, FAD, midpoint, duration)))
	SQS_1_df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = SQS_1_df$year))
	SQS_1_df$order <- rep(0, nrow(SQS_1_df))
	SQS_1_df$method <- rep("interpolated", nrow(SQS_1_df))
	SQS_1_df$flavour <- rep("SQS1", nrow(SQS_1_df))
	SQS_1_df$estimator <- rep("SQS", nrow(SQS_1_df))
	# SQS_1_df %>% filter(year == 2016 & quorum == "Quorum 0.4") %>% ggplot(., aes(x = Midpoint, y = Subsampled.diversity)) + geom_line(color = "darkgrey") + geom_point(cex = 2) + scale_x_reverse() + geom_vline(xintercept = 66, lty = 2)
	# filter(SQS_1_df, bin == "Maastrichtian" & quorum == "Quorum 0.4") %>% ggplot(., aes(x = year, y = Subsampled.diversity)) + geom_line(color = "darkgrey") + geom_point(cex = 2)
}

#### flatten SQS Perl script without collsperref into tidy data frame ####
if (run_SQS_Perl_2 == TRUE) {
	SQS_2_df <- melt(real.sqs.2)
	SQS_2_df <- spread(SQS_2_df, Var2, value)
	SQS_2_df <- tbl_df(rename(SQS_2_df, bin = Var1, quorum = Var3, year = Var4, region = L1))
	SQS_2_df <- tbl_df(cbind(SQS_2_df, select(intervals[match(SQS_2_df$bin, intervals$bin), ], LAD, FAD, midpoint, duration)))
	SQS_2_df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = SQS_2_df$year))
	SQS_2_df$order <- rep(0, nrow(SQS_2_df))
	SQS_2_df$method <- rep("interpolated", nrow(SQS_2_df))
	SQS_2_df$flavour <- rep("SQS2", nrow(SQS_2_df))
	SQS_2_df$estimator <- rep("SQS", nrow(SQS_2_df))
}

#### flatten lambda-5 recursive list into tidy data frame
if (run_lambda_5 == TRUE) {

	lambda5_df <- melt(real.lambda5) %>% spread(., Var2, value) %>% rename(bin = Var1, year = L2, region = L1) %>% tbl_df()
	lambda5_df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = lambda5_df$year))
	lambda5_df <- tbl_df(cbind(lambda5_df, select(intervals[match(lambda5_df$bin, intervals$bin), ], LAD, FAD, midpoint, duration)))

}

#### flatten Chao recursive list into tidy data frame
if (run_chao == TRUE) {

	chao_df <- melt(real.chao) %>% spread(., Var2, value) %>% rename(bin = Var1, year = L2, region = L1) %>% tbl_df()
	chao_df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = chao_df$year))
	chao_df <- tbl_df(cbind(chao_df, select(intervals[match(chao_df$bin, intervals$bin), ], LAD, FAD, midpoint, duration)))

}

#### flatten iNEXT recursive list into tidy data frame
if (run_iNEXT == TRUE) {

	tmp <- real.inext
	for (rgn in 1:length(regions)) {
		for (yr in 1:length(times)) {
			for (int in 1:nrow(intervals)) {
				tmp[[rgn]][[yr]][[int]] <- as.data.frame(real.inext[[rgn]][[yr]][[int]])
				tmp[[rgn]][[yr]][[int]]$region <- rep(names(regions)[rgn], nrow(real.inext[[rgn]][[yr]][[int]]))
				tmp[[rgn]][[yr]][[int]]$year <- rep(names(times)[yr], nrow(real.inext[[rgn]][[yr]][[int]]))
				tmp[[rgn]][[yr]][[int]]$bin <- rep(intervals$bin[int], nrow(real.inext[[rgn]][[yr]][[int]]))
			}
		}
	}
	tmp <- do.call('rbind', unlist(unlist(tmp, recursive = F), recursive = FALSE))
	rownames(tmp) <- NULL
	tmp <- tbl_df(tmp)
	iNEXT_df <- tbl_df(cbind(tmp, select(intervals[match(tmp$bin, intervals$bin), ], LAD, FAD, midpoint, duration)))
	rm(list = "tmp")
	iNEXT_df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = iNEXT_df$year))
	iNEXT_df$flavour <- rep("iNEXT", nrow(iNEXT_df))
	iNEXT_df$estimator <- rep("SQS", nrow(iNEXT_df))
}

if (run_CR == TRUE) {

	tmp <- real.CR
	for (rgn in 1:length(regions)) {
		for (yr in 1:length(times)) {
			for (int in 1:nrow(intervals)) {
				tmp[[rgn]][[yr]][[int]] <- as.data.frame(real.CR[[rgn]][[yr]][[int]])
				tmp[[rgn]][[yr]][[int]]$region <- rep(names(regions)[rgn], nrow(real.CR[[rgn]][[yr]][[int]]))
				tmp[[rgn]][[yr]][[int]]$year <- rep(names(times)[yr], nrow(real.CR[[rgn]][[yr]][[int]]))
				tmp[[rgn]][[yr]][[int]]$bin <- rep(intervals$bin[int], nrow(real.CR[[rgn]][[yr]][[int]]))
			}
		}
	}
	tmp <- do.call('rbind', unlist(unlist(tmp, recursive = F), recursive = FALSE))
	rownames(tmp) <- NULL
	tmp <- tbl_df(tmp)
	CR_df <- tbl_df(cbind(tmp, select(intervals[match(tmp$bin, intervals$bin), ], LAD, FAD, midpoint, duration)))
	rm(list = "tmp")
	CR_df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = CR_df$year))
	CR_df$flavour <- rep("iNEXT", nrow(CR_df))
	CR_df$estimator <- rep("CR", nrow(CR_df))
}

#### flatten TRiPS recursive list into tidy data frame
if (run_TRiPS == TRUE) {
	TRiPS_df <- melt(real.trips); colnames(TRiPS_df) <- c("TRiPS_variable","TRiPS_statistic","bin","year","value","region")
	TRiPS_df$Var <- paste(TRiPS_df$TRiPS_variable, TRiPS_df$TRiPS_statistic); TRiPS_df <- select(TRiPS_df, -TRiPS_variable, -TRiPS_statistic)
	TRiPS_df <- spread(TRiPS_df, Var, value)
	TRiPS_df <- tbl_df(cbind(TRiPS_df, select(intervals[match(TRiPS_df$bin, intervals$bin), ], LAD, FAD, midpoint, duration)))
	TRiPS_df$year <- as.integer(gsub(pattern = "year_", replacement = "", ignore.case = T, x = TRiPS_df$year))
	TRiPS_df$order <- rep(0, nrow(TRiPS_df))
	TRiPS_df$estimator <- rep("TRiPS", nrow(TRiPS_df))
	TRiPS_df$flavour <- rep("TRiPS", nrow(TRiPS_df))
	TRiPS_df <- rename(TRiPS_df, qD = `Estimated richness MLE`, qD.LCL = `Estimated richness lower CI`, qD.UCL = `Estimated richness upper CI`)
}

SQS_df <- full_join(SQS_1_df, SQS_2_df)
