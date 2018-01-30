#### see how well Good's u performs in simulations ####
moar_sizes <- lseq(from = 1, to = 10000, length.out = 40) %>% round() %>% unique()
moar_nreps <- 500
tim <- list()
for (k in 1:nrow(pools)) {
	distr <- setNames(underlying[[k]], paste("sp", 1:length(underlying[[k]]), sep = "_"))
	distr <- distr/sum(distr)
	tim[[k]] <- lapply(1:length(moar_sizes), function(i) {
		s1 <- lapply(1:moar_nreps, function(x) table(sample(names(distr), size = moar_sizes[i], replace = T, prob = distr)))
		u_reps <- do.call("c", mclapply(1:moar_nreps, function(x) sampleCoverage(s1[[x]]), mc.cores = detectCores(), mc.preschedule = TRUE))
		mr_reps <- do.call("c", mclapply(1:moar_nreps, function(x) multitonRatio(s1[[x]]), mc.cores = detectCores(), mc.preschedule = TRUE))
		true_u_reps <- do.call("c", mclapply(1:moar_nreps, function(x) sum(distr[which(names(distr) %in% names(s1[[x]]))]), mc.cores = detectCores(), mc.preschedule = TRUE))
		error_reps <- u_reps - true_u_reps
		u <- as.numeric(quantile(u_reps, c(0.025,0.5,0.975), na.rm = T))
		mr <- as.numeric(quantile(mr_reps, c(0.025,0.5,0.975), na.rm = T))
		true_u <- as.numeric(quantile(true_u_reps, c(0.025,0.5,0.975), na.rm = T))
		error <- as.numeric(quantile(error_reps, c(0.025,0.5,0.975), na.rm = T))
		u.LCL <- u[1]; u.UCL <- u[3]; u <- u[2]
		mr.LCL <- mr[1]; mr.UCL <- mr[3]; mr <- mr[2]
		true_u.LCL <- true_u[1]; true_u.UCL <- true_u[3]; true_u <- true_u[2]
		error.LCL <- error[1]; error.UCL <- error[3]; error <- error[2]
		tibble(sample_size = moar_sizes[i], assemblage = pools$assemblage[k], true_richness = pools$spp[k], sd = pools$sd[k], mr.LCL = mr.LCL, mr = mr, mr.UCL = mr.UCL, u.LCL = u.LCL, u = u, u.UCL = u.UCL, true_u.LCL = true_u.LCL, true_u = true_u, true_u.UCL = true_u.UCL, error.LCL = error.LCL, error = error, error.UCL = error.UCL)
	}); tim[[k]] <- do.call("rbind", tim[[k]])
}; tim <- do.call("rbind", tim)

#### plot Good's u simulation results ####
pdf("./simulation-figures/Fig. S<coverage-simulations-bivariate>.pdf", width = simfigheight, height = simfigwidth * 0.7)
tim %>%
	filter(true_richness %in% c(50,400)) %>%
	filter(sd != 1.5) %>%
	mutate(sd = as.factor(sd)) %>% mutate(sd = recode(sd, "0" = "flat", "1" = "sigma = 1", "2" = "sigma = 2")) %>%
	mutate(true_richness = paste("True Richness =", true_richness)) %>% mutate(true_richness = factor(true_richness, levels = sort(unique(true_richness), decreasing = T))) %>%
	ggplot(aes(x = true_u, y = u)) +
	geom_point(alpha = 0.5) +
	coord_equal() +
	geom_abline(slope = 1, intercept = 0, colour = "blue") +
	# theme_bw() +
	xlim(0,1) +
	ylim(0,1) +
	xlab("True Coverage") +
	ylab("Good's u") +
	facet_grid(true_richness~sd)
graphics.off()

pdf("./simulation-figures/Fig. S<coverage-simulations-rarefaction-curves>.pdf", width = simfigheight, height = simfigwidth * 0.7)
tim %>%
	filter(true_richness %in% c(50,400)) %>%
	filter(sd != 1.5) %>%
	mutate(sd = as.factor(sd)) %>% mutate(sd = recode(sd, "0" = "flat", "1" = "sigma = 1", "2" = "sigma = 2")) %>%
	mutate(true_richness = paste("True Richness =", true_richness)) %>% mutate(true_richness = factor(true_richness, levels = sort(unique(true_richness), decreasing = T))) %>%
	ggplot(aes(x = sample_size, y = true_u, ymin = true_u.LCL, ymax = true_u.UCL, fill = "True coverage")) +
	geom_line(aes(colour = "True coverage")) +
	geom_ribbon(aes(x = sample_size, ymin = true_u.LCL, ymax = true_u.UCL), alpha = 0.2, size = 0.1, show.legend = FALSE) +
	geom_line(aes(x = sample_size, y = u, colour = "Good's u")) +
	geom_ribbon(aes(x = sample_size, ymin = u.LCL, ymax = u.UCL, fill = "Good's u"), alpha = 0.2, size = 0.1, show.legend = FALSE) +
	scale_x_continuous(trans = "log10") + #, breaks = c(1,10,100,1000)
	# annotation_logticks(side="tb") +
	labs(colour = 'Coverage') +
	ylim(c(0,1)) +
	# theme_bw() +
	theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal") +
	xlab("Sample Size") + ylab("Coverage") +
	facet_grid(true_richness~sd)
graphics.off()

pdf("./simulation-figures/Fig. S<coverage-simulations-deviation>.pdf", width = simfigheight, height = simfigwidth)
tim %>% ggplot(aes(x = sample_size, y = error, ymin = error.LCL, ymax = error.UCL)) +
	geom_line() +
	geom_ribbon(aes(x = sample_size, ymin = error.LCL, ymax = error.UCL), alpha = 0.2, size = 0.1, show.legend = FALSE) +
	scale_x_log10() +
	# theme_bw() +
	# annotation_logticks(side="tb") +
	ylab("Deviation from true coverage") +
	xlab("Sample Size") +
	facet_grid(sd~true_richness)
graphics.off()

