#### set parameters for simulations ####
minss <- 2 #minimum sample size
maxss <- 10000 #maximum sample size
knots <- 50 #number of logarithmically-spaced intermediate steps between min and max
richnesses <- c(50,100,200,400) #true richnesses
sd_values <- c(0,1,1.5,2) #true evenness

# function to generate a logarithmically-spaced sequence of sample sizes
lseq <- function(from = 10, to = 10000, length.out = 5) {
	exp(seq(log(from), log(to), length.out = length.out))
}

#make object with information about communities
pools <- tibble(
	assemblage = LETTERS[1:(length(richnesses)*length(sd_values))],
	spp = do.call("c", lapply(richnesses, function(x) rep(x, length(sd_values)))),
	sd = c(rep(sd_values, length(richnesses))))

spp <- pools$spp; sd <- pools$sd; names(sd) <- names(spp) <- pools$assemblage

#sample sizes to use
sizes <- lseq(from = minss, to = maxss, length.out = knots) %>% round() %>% unique()

#generate underlying frequency distributions for simulated communities
underlying <- lapply(1:length(spp), function(k) {
	if (sd[k] > 0) {
		lapply(1:10000, function(x) sort(exp(rnorm(spp[k], sd = sd[k])))) %>% do.call("rbind", .) %>% colMeans()
	} else if (sd[k] == 0) {
		rep(1/spp[k], spp[k]) #perfectly uniform SAD
	}
}); names(underlying) <- pools$assemblage

#choose quorums and quotas
quora <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,0.999)
quotas <- c(50,100,200,400,800,1600)#,3200,6400,12800)
nreps <- 200 #number of replicate simulations
plotting_quorum <- 0.9 #quorum for plotting SQS results when only one quorum level shown
plotting_quota <- 800 #quota for plotting CR results when only one quota shown

#### do main simulations #####
rem <- list()
for (k in 1:length(spp)) {
	rem[[k]] <- list()
	for (i in 1:length(sizes)) {
		cat(paste("Assemblage",names(spp)[k], "richness", spp[k],"sample size",sizes[i],"\n"))

		#generate replicate samples of abundance data
		obs <- lapply(1:nreps, function(z) as.vector(unname(table(sample(1:spp[k], size = sizes[i], replace = T, prob = underlying[[k]])))))

		#calculate SQS (iNEXT) q = 0
		abbot <- tbl_df(do.call("rbind", lapply(quora, function(q) tryCatch(quantile(do.call("c", mclapply(1:nreps, function(x) {
			sqs_inext <- estimateD(x = obs[[x]], datatype = "abundance", base = "coverage", level = q, conf = NULL)
			if (sqs_inext$method == "extrapolated") {return(NA)} else {return(as.numeric(sqs_inext[1,4]))}
		}, mc.cores = detectCores(), mc.preschedule = TRUE)), probs = c(0.025,0.5,0.975), na.rm = F), error = function(err) c(NA,NA,NA)))))
		colnames(abbot) <- c("qD.LCL","qD","qD.UCL")
		abbot$method <- "sqs_inext"
		abbot$quorum <- quora
		abbot$sample_size <- sizes[i]
		abbot$true_richness <- spp[k]
		abbot$sd <- sd[k]
		abbot$assemblage <- names(spp)[k]

		#calculate CR
		farrage <- tbl_df(do.call("rbind", lapply(quotas, function(q) tryCatch(quantile(do.call("c", mclapply(1:nreps, function(x) {
			cr_inext <- estimateD(x = obs[[x]], datatype = "abundance", base = "size", level = q, conf = NULL)
			if (cr_inext$method == "extrapolated") {return(NA)} else {return(as.numeric(cr_inext[1,4]))}
		}, mc.cores = detectCores(), mc.preschedule = TRUE)), probs = c(0.025,0.5,0.975)), error = function(err) c(NA,NA,NA)))))
		colnames(farrage) <- c("qD.LCL","qD","qD.UCL")
		farrage$method <- "cr_inext"
		farrage$quota <- quotas
		farrage$sample_size <- sizes[i]
		farrage$true_richness <- spp[k]
		farrage$sd <- sd[k]
		farrage$assemblage <- names(spp)[k]

		#calculate other richness estimators
		may <- do.call("rbind", mclapply(1:nreps, function(x) {
			otherEstimators(s1 = obs[[x]])
		}, mc.cores = detectCores(), mc.preschedule = TRUE))

		may1 <- apply(may, MARGIN = 2, FUN = function(l) {
			if (all(is.na(l))) {
				c(NA,NA,NA)
			} else {
				quantile(l,probs=c(0.025,0.5,0.975), na.rm = T)
			}
		})
		may2 <- melt(may1) %>% spread(Var1, value); colnames(may2) <- c("method","qD.LCL","qD","qD.UCL")
		may2 <- tbl_df(may2)
		may2$sample_size <- sizes[i]
		may2$true_richness <- spp[k]
		may2$sd <- sd[k]
		may2$assemblage <- names(spp)[k]

		#calculate doubleton coverage
		farron <- t(do.call("c", mclapply(1:nreps, function(x) sampleCoverage(obs[[x]]))) %>% quantile(., probs = c(0.025,0.5,0.975), na.rm = T))
		farron2 <- tbl_df(farron)
		colnames(farron2) <- c("SC.LCL","SC","SC.UCL")
		farron2$sample_size <- sizes[i]
		farron2$true_richness <- spp[k]
		farron$sd <- sd[k]
		farron2$assemblage <- names(spp)[k]

		#merge tibbles
		coalition <- full_join(abbot, farrage) %>% full_join(., may2) %>% full_join(., farron2)

		f1 <- tbl_df(t(quantile(sapply(obs, function(x) sum(x == 1)), probs = c(0.025,0.5,0.975), na.rm = T))); colnames(f1) <- c("f1.LCL","f1","f1.UCL"); rownames(f1) <- NULL
		f2 <- tbl_df(t(quantile(sapply(obs, function(x) sum(x == 2)), probs = c(0.025,0.5,0.975), na.rm = T))); colnames(f2) <- c("f2.LCL","f2","f2.UCL"); rownames(f2) <- NULL
		fm <- tbl_df(t(quantile(sapply(obs, function(x) sum(x > 1)), probs = c(0.025,0.5,0.975), na.rm = T))); colnames(fm) <- c("fm.LCL","fm","fm.UCL"); rownames(fm) <- NULL
		coalition <- tbl_df(cbind(coalition, f1, f2, fm))
		rem[[k]][[i]] <- filter(coalition, method != "u")

	}; rem[[k]] <- do.call("rbind", rem[[k]])
}; rem <- do.call("rbind", rem)

system("say Analysis finished, fool!")

simdf <- rem
simdf$assemblage <- paste(simdf$assemblage, simdf$true_richness, simdf$sd)

simdf$sd_labels <- factor(simdf$sd, labels = c("flat", "lognorm (sd = 1)", "lognorm (sd = 1.5)", "lognorm (sd = 2)"))

simdf$method <- factor(simdf$method, levels = c(
	"raw",
	"trips",
	"chao",
	"lambda",
	"sqs_inext",
	"cr_inext"
))

simdf$method_labels <- factor(simdf$method, labels = c(
	"Face-value Counts",
	"TRiPS",
	"Chao1",
	"Lambda-5",
	"SQS",
	"CR"
))

simdf$true_richness_labels <- factor(paste(simdf$true_richness, "species"))

simdf$true_richness_labels <- factor(simdf$true_richness_labels, labels = c("50 species","100 species","200 species","400 species"))

coverage_vars <- str_subset(pattern = "sqs", string = unique(simdf$method))

#calculate relative richness ratios (relative to least diverse assemblage)
ass <- unique(simdf$assemblage)
for (j in 1:length(sd_values)) {
	for (i in 1:length(richnesses)) {
		simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "true_richness.ratio"] <- filter(simdf, true_richness == richnesses[i] & sd == sd_values[j])$true_richness/filter(simdf, true_richness == richnesses[1] & sd == sd_values[j])$true_richness
		simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD.ratio"]            <- filter(simdf, true_richness == richnesses[i] & sd == sd_values[j])$qD/filter(simdf, true_richness == richnesses[1] & sd == sd_values[j])$qD
		simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD.LCL.ratio"]        <- simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD.ratio"]/simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD"] * simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD.LCL"]
		simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD.UCL.ratio"]        <- simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD.ratio"]/simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD"] * simdf[which(simdf$true_richness == richnesses[i] & simdf$sd == sd_values[j]), "qD.UCL"]
	}
}

