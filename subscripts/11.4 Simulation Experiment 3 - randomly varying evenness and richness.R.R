#### simulation tests looking at effect of evenness ####
spp1 <- 100
reps <- 500
n <- c(25,50,100,500,1000,5000,10000)
plotting_n <- c(50,500,5000)
quora <- c(0.1,0.4,0.7,0.9,0.99,0.999)
bob1 <- jim1 <- list()
for (i in 1:length(n)) {
	sds <- runif(min = 1, max = 2, n = reps)
	obs <- lapply(1:reps, function(k) as.vector(unname(table(sample(1:spp1, size = n[i], replace = T, prob = exp(rnorm(spp1, sd = sds[[k]])))))))

	bob <- mclapply(1:reps, function(x) {
		tibble(method = "sqs", qD = sapply(quora, function(q) simpleSQS(n = obs[[x]], quorum = q)), quorum = quora, sd = sds[x], n = n[i], truerichness = spp1) #[-which.max(obs[[x]])]
	}, mc.cores = detectCores(), mc.preschedule = TRUE)
	bob1[[i]] <- bob %>% do.call("rbind", .) %>% tbl_df()

	jim <- mclapply(1:reps, function(x) {
		tmp <- otherEstimators(s1 = obs[[x]])
		tibble(method = names(tmp), qD = as.numeric(tmp), sd = sds[x], n = n[i], truerichness = spp1)
	}, mc.cores = detectCores(), mc.preschedule = TRUE)
	jim1[[i]] <- jim %>% do.call("rbind", .) %>% tbl_df()
}

jim2 <- do.call("rbind", jim1); bob2 <- do.call("rbind", bob1)

u <- jim2 %>% filter(method == "u") %>% dplyr::rename(u = qD) %>% select(-method)
jim2 <- jim2 %>% filter(method != "u")
jim3 <- full_join(jim2, bob2)
jim3 <- full_join(jim3, u)
jim3$method <- as.factor(jim3$method)

jim3$method = factor(jim3$method, levels = c("raw","sqs","multi","trips","chao","lambda")) #reorder method factor

tmp <- jim3 %>% filter(method == "sqs") %>% mutate(method = as.character(method))
tmp$method <- paste(tmp$method, tmp$quorum)

tmp <- tmp %>% spread(method, qD) %>% select(-quorum) %>% gather("method", "qD", 5:10) %>% mutate(method = as.factor(method))

jim4 <- full_join(filter(jim3, method != "sqs"), tmp) %>% select(-quorum)
jim4$method = factor(jim4$method, levels = c("raw", str_subset(pattern = "sqs", string = unique(jim4$method)), "multi","trips","chao","lambda")) #reorder method factor

#### plot results of simulations fixing richness and varying evenness ####
pdf("./simulation-figures/Fig. S<coverage-sample-size-evenness-relationship>.pdf", width = 8, height = 2.8)
p <- list()
for (i in 1:(length(plotting_n))) {
	x <- plotting_n[i] #because ggplot gets confused by n in tibble and in vector
	p[[i]] <- jim3 %>% filter((method %in% c("raw")) & n == x) %>%
		mutate(n = paste("Sample Size =", n)) %>% mutate(n = factor(n, levels = sort(unique(n), decreasing = T))) %>%
		ggplot(aes(x = sd, y = u)) +
		geom_point(alpha = 0.8, size = 0.5) +
		# theme_bw() +
		facet_grid(~n) +
		geom_smooth(method = "lm") +
		labs(x = "Evenness (SD)", y = "Coverage (Good's u)") +
		theme(legend.position = "none")
}
p1 <- do.call(grid.arrange,c(p,list(nrow = 1)))
graphics.off()


#### varying evenness and richness, too ####
reps <- 500
n <- c(25,50,100,500,1000,5000,10000)
quora <- c(0.1,0.4,0.7,0.9,0.99,0.999)
bill1 <- tom1 <- list()
for (i in 1:length(n)) {
	sds <- runif(min = 1, max = 2, n = reps)
	spps <- round(exp(rnorm(reps, mean = 5, sd = 1)))
	obs <- lapply(1:reps, function(k) as.vector(unname(table(sample(1:spps[k], size = n[i], replace = T, prob = exp(rnorm(spps[k], sd = sds[k])))))))

	bill <- mclapply(1:reps, function(x) {
		tibble(method = "sqs", qD = sapply(quora, function(q) simpleSQS(n = obs[[x]][-which.max(obs[[x]])], quorum = q)), quorum = quora, sd = sds[x], n = n[i], truerichness = spps[x]) #
	}, mc.cores = detectCores(), mc.preschedule = TRUE)
	bill1[[i]] <- bill %>% do.call("rbind", .) %>% tbl_df()

	tom <- mclapply(1:reps, function(x) {
		tmp <- otherEstimators(s1 = obs[[x]])
		tibble(method = names(tmp), qD = as.numeric(tmp), sd = sds[x], n = n[i], truerichness = spps[x])
	}, mc.cores = detectCores(), mc.preschedule = TRUE)
	tom1[[i]] <- tom %>% do.call("rbind", .) %>% tbl_df()
}

tom2 <- do.call("rbind", tom1); bill2 <- do.call("rbind", bill1)
bill2$qD <- as.numeric(bill2$qD); bill2$sd <- as.numeric(bill2$sd); bill2$quorum <- as.numeric(bill2$quorum); bill2$quorum <- as.numeric(bill2$quorum); bill2$n <- as.numeric(bill2$n); bill2$truerichness <- as.numeric(bill2$truerichness); bill2$method <- as.character(bill2$method)

u <- tom2 %>% filter(method == "u") %>% dplyr::rename(u = qD) %>% select(-method)
tom2 <- tom2 %>% filter(method != "u")
tom3 <- full_join(tom2, bill2)
tom3 <- full_join(tom3, u)
tom3$method <- as.factor(tom3$method)

tom3$method = factor(tom3$method, levels = c("raw","sqs","multi","trips","chao","lambda")) #reorder method factor

tmp <- tom3 %>% filter(method == "sqs") %>% mutate(method = as.character(method))
tmp$method <- paste(tmp$method, tmp$quorum)

tmp <- tmp %>% spread(method, qD) %>% select(-quorum) %>% gather("method", "qD", 5:10) %>% mutate(method = as.factor(method))

tom4 <- full_join(filter(tom3, method != "sqs"), tmp) %>% select(-quorum)
tom4$method = factor(tom4$method, levels = c("raw","trips","chao","lambda",str_subset(pattern = "sqs", string = unique(tom4$method)), "multi")) #reorder method factor

pdf("./simulation-figures/Fig. S<rf-sd-vs-qD-rnorm>.pdf", width = simfigwidth * 1.2, height = simfigheight * 1.3)
tom4 %>%
	filter(!grepl("sqs", method), n %in% plotting_n) %>%
	mutate(method = recode(method, "raw" = "Face-value\nCounts", "chao" = "Chao1", "lambda" = "Lambda-5", "trips" = "TRiPS")) %>%
	ggplot(aes(x = sd, y = qD)) +
	geom_point(aes(colour = truerichness), alpha = 0.9, size = 0.5) +
	scale_color_viridis(trans = "log", breaks = c(1, 10, 50, 400)) +
	scale_y_log10() +
	stat_smooth(method = "lm") +
	theme(legend.position = "bottom") +
	labs(x = "Evenness (SD)", y = "Size-standardised Richness", colour = "True Richness") +
	facet_wrap(method~n, scales = "free", labeller = labeller(n = label_both))
graphics.off()

pdf("./simulation-figures/Fig. <rf-trueD-vs-qD-rnorm>.pdf", width = simfigwidth * 0.9, height = simfigheight)
plotting_n.bak <- plotting_n
plotting_n <- c(50,500,5000)
tom4 %>%
	filter(!grepl("sqs", method), n %in% plotting_n) %>%
	mutate(method = recode(method, "raw" = "Face-value\nCounts", "chao" = "Chao1", "lambda" = "Lambda-5", "trips" = "TRiPS")) %>%
	mutate(n = paste("Quota =", n)) %>% mutate(n = factor(n, levels = c("Quota = 50", "Quota = 500", "Quota = 5000"))) %>%
	ggplot(aes(x = truerichness, y = qD)) +
	geom_point(aes(colour = sd), alpha = 0.9, size = 0.5) +
	scale_color_viridis() +
	scale_y_log10() +
	scale_x_log10() +
	# geom_smooth(se = F) +
	geom_abline(slope = 1) +
	theme(legend.position = "bottom") +
	labs(colour = "SD", x = "True Richness", y = "Size-standardised Richness") +
	facet_grid(method~n, scales = "free")
plotting_n <- plotting_n.bak
graphics.off()
