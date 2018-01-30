#### settings ####
quora <- c(seq(from = 0, to = 0.9, by = 0.1), 0.95, 0.99, 0.999, 0.9999, 1)
plotting_quora <- quora#seq(from = 0, to = 1, by = 0.1)
ss <- sapply(1:length(spp), function(k) {simdf %>% filter(sd == sd[k], true_richness == spp[k] & SC.LCL == 1) %>% select(sample_size) %>% min})
ss <- ss * 3
trials <- 2000

#### analyses ####
qf <- list()
for (k in 1:length(spp)) {
	cat(paste("spp =", spp[k], "sd =", sd[k]), "\n")
	occurrences <- sample(1:spp[k], size = ss[k], replace = T, prob = underlying[[k]])
	tmp <- tbl_df(cbr(occurrences = occurrences, quora = quora, trials = trials))
	tmp$quorum <- quora
	tmp <- gather(tmp, "method", "qD", 1:3)
	tmp$true_richness <- spp[k]
	tmp$sd <- sd[k]
	qf[[k]] <- tmp
}; qf <- qf %>% do.call("rbind", .)
qf$method <- factor(qf$method, levels = c(
	"raw",
	# "trips",
	"chao",
	"lambda"
))


#### coverage-standardised richness ####
reps = 500
# reps = 4
sds <- runif(min = 1, max = 2, n = reps)
spps <- round(exp(rnorm(reps, mean = 5, sd = 1)))
obs <- lapply(1:reps, function(k) sample(1:spps[k], size = 5000, replace = T, prob = exp(rnorm(spps[k], sd = sds[k]))))
moarquora <- c(seq(from = 0.1, to = 0.9, by = 0.1), 0.99)

neil <- lapply(1:reps, function(k) {
	tmp <- tbl_df(cbr(occurrences = obs[[k]], quora = moarquora, trials = 500))
	tmp$quorum <- moarquora
	tmp <- gather(tmp, "method","qD", 1:3)
	tmp$sd <- sds[k]
	tmp$true_richness <- spps[k]
	tmp
})
neil1 <- neil %>% do.call("rbind", .) %>% tbl_df()
neil1$method = factor(neil1$method, levels = c("raw","chao","lambda")) #reorder method factor


simfigwidth <- 8; simfigheight <- 8


pdf("./simulation-figures/Fig. S<qf-sd-vs-qD-rnorm>.pdf", width = simfigwidth, height = simfigheight)
neil1 %>% filter(quorum %in% c(0.1,0.5,0.9,0.99) & !is.na(qD)) %>%
	ggplot(aes(x = sd, y = qD)) +
	geom_point(aes(colour = true_richness), alpha = 0.9, size = 0.5) +
	scale_color_viridis(trans = "log", breaks = c(50,500)) +
	scale_y_log10() +
	stat_smooth(method = "lm") +
	labs(colour = "True Richness", x = "Evenness (SD)", y = "Coverage-standardised Richness") +
	facet_wrap(method~quorum, scales = "free") +
	# theme_minimal() +
	theme(legend.position = "bottom")
graphics.off()

pdf("./simulation-figures/Fig. <qf-trueD-vs-qD-rnorm>.pdf", width = simfigwidth, height = simfigheight * 1.1)
neil1 %>% filter(quorum %in% c(0.1,0.5,0.9,0.99) & !is.na(qD)) %>%
	mutate(method = recode(method, "raw" = "Face-value\nCounts", "chao" = "Chao1", "lambda" = "Lambda-5")) %>%
	mutate(quorum = paste("Quorum = ", quorum, sep = "")) %>%
	ggplot(aes(x = true_richness, y = qD)) +
	geom_point(aes(colour = sd), alpha = 0.9, size = 0.5) +
	scale_color_viridis() +
	scale_y_log10() +
	scale_x_log10() +
	geom_abline(slope = 1) +
	labs(colour = "Evenness (SD)", x = "True Richness", y = "Coverage-standardised Richness") +
	facet_grid(method~quorum, scales = "free") +
	theme(legend.position = "bottom")
graphics.off()

#### multiton-ratio-standardised richness ####
moar_mr <- c(seq(from = 0.05, to = 0.95, by = 0.05))
joschi <- lapply(1:reps, function(k) {
	tmp <- tbl_df(mbr(occurrences = obs[[k]], mr = moar_mr, trials = 10))
	tmp$mr <- moar_mr
	tmp <- gather(tmp, "method","qD", 1:3)
	tmp$sd <- sds[k]
	tmp$true_richness <- spps[k]
	tmp
})
joschi1 <- joschi %>% do.call("rbind", .) %>% tbl_df()
joschi1$method = factor(joschi1$method, levels = c("raw","chao","lambda")) #reorder method factor
joschi1$mr_factor <- as.factor(joschi1$mr)

pdf("./simulation-figures/Fig. <mf-sd-vs-qD-rnorm>.pdf")
joschi1 %>% filter(mr_factor %in% c(0.2,0.4,0.7,0.9) & !is.na(qD)) %>%
	ggplot(aes(x = sd, y = qD)) +
	geom_point(aes(colour = true_richness), alpha = 0.9, size = 0.5) +
	scale_color_viridis(trans = "log") +
	scale_y_log10() + stat_smooth(method = "lm") +
	labs(colour = "True Richness", x = "Evenness", y = "Estimated Richness") +
	facet_grid(method~mr_factor, scales = "free")
graphics.off()

pdf("./simulation-figures/Fig. <mf-trueD-vs-qD-rnorm>.pdf", width = width, height = height * 0.8)
joschi1 %>% filter(mr_factor %in% c(0.2,0.4,0.7,0.9) & !is.na(qD)) %>%
	ggplot(aes(x = true_richness, y = qD)) +
	geom_point(aes(colour = sd), alpha = 0.9, size = 0.5) +
	scale_color_viridis() + scale_y_log10(limits = c(1,1000)) +
	scale_x_log10(limits = c(1,1000)) +
	geom_abline(slope = 1) +
	labs(colour = "SD", x = "True Richness", y = "Estimated Richness") +
	facet_grid(method~mr_factor, scales = "free")
graphics.off()


#### important plots ####
width = 8; height = 8

pdf("./simulation-figures/Fig. <coverage-rarefied-estimators>.pdf", width = simfigwidth, height = simfigheight * 0.5)
qf %>% filter(quorum %in% quora) %>%
	filter(true_richness %in% c(400)) %>%
	mutate(method = recode(method, "raw" = "Face-value\nCounts", "chao" = "Chao1", "lambda" = "Lambda-5")) %>%
	mutate(sd = as.factor(sd)) %>% mutate(sd = recode(sd, "0" = "Flat", "1" = "1", "1.5" = "1.5", "2" = "2")) %>%
	mutate(true_richness = paste("True Richness =", true_richness)) %>% mutate(true_richness = factor(true_richness, levels = sort(unique(true_richness), decreasing = F))) %>%
	ggplot(aes(x = quorum, y = qD, colour = as.factor(sd))) +
	geom_line() +
	scale_y_continuous(trans = "log2", breaks = c(1,2,5,10,25,50,100,200,400,800)) +
	scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
	scale_color_brewer(palette = "Set1") +
	labs(x = "Sample Coverage", y = "Coverage-standardised Richness", colour = "Evenness (lognormal variance)") +
	theme(legend.position = "bottom") +
	facet_grid(true_richness~method, scales = "free")
graphics.off()

pdf("./simulation-figures/Fig. S<qf-q-vs-qD-via-sd>.pdf", width = simfigwidth, height = simfigheight)
qf %>% filter(quorum %in% plotting_quora) %>%
	ggplot(aes(x = quorum, y = qD, colour = as.factor(true_richness))) +
	geom_line() +
	scale_y_log10() +
	scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
	labs(x = "Quorum", y = "Coverage-standardised Richness", colour = "True Richness") +
	theme(legend.position = "bottom") +
	facet_grid(method~sd, scales = "free")
graphics.off()

pdf("./simulation-figures/Fig. S<coverage-rarefied-estimators-superimposed>.pdf", width = simfigwidth, height = simfigheight * 0.4)
qf %>% filter(quorum %in% plotting_quora) %>%
	filter(true_richness %in% c(400)) %>%
	mutate(method = recode(method, "raw" = "Face-value Counts", "chao" = "Chao1", "lambda" = "Lambda-5")) %>%
	mutate(sd = as.factor(sd)) %>% mutate(sd = recode(sd, "0" = "Flat", "1" = "sigma = 1", "1.5" = "sigma = 1.5", "2" = "sigma = 2")) %>%
	mutate(true_richness = paste("True Richness =", true_richness)) %>% mutate(true_richness = factor(true_richness, levels = sort(unique(true_richness), decreasing = F))) %>%
	ggplot(aes(x = quorum, y = qD, colour = as.factor(method))) +
	geom_line() +
	scale_y_log10() +
	scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
	theme(legend.position = "bottom") +
	labs(x = "Quorum", y = "Coverage-standardised Richness", colour = "Richness Estimator") +
	facet_grid(true_richness~sd, scales = "free")
graphics.off()

pdf("./simulation-figures/Fig. S<qf-sd-vs-qD>.pdf", width = simfigwidth * 1.2, height = simfigheight)
qf %>% filter(quorum %in% c(0.1,0.5,0.9,0.99)) %>%
	ggplot(aes(x = sd, y = qD, colour = as.factor(true_richness))) +
	geom_point() +
	geom_line() +
	scale_y_continuous(trans = "log2") +
	scale_x_continuous(breaks = c(0, 1, 1.5, 2)) +
	theme(legend.position = "bottom") +
	labs(x = "Evenness (SD)", y = "Coverage-standardised Richness", colour = "True Richness") +
	facet_wrap(method~quorum, scales = "free", nrow = 3)
graphics.off()

pdf("./simulation-figures/Fig. <qf-sd-vs-qD-superimposed>.pdf", width = 5, height = 4)
qf %>% filter(quorum %in% c(0.9,0.99)) %>%
	filter(true_richness %in% c(200,400)) %>%
	mutate(method = recode(method, "raw" = "Face-value Counts", "chao" = "Chao1", "lambda" = "Lambda-5")) %>%
	ggplot(aes(x = sd, y = qD, colour = as.factor(method), lty = as.factor(true_richness))) +
	geom_point() +
	geom_line() +
	scale_y_continuous(trans = "log2", breaks = c(1,2,5,10,25,50,100,200,400,800)) +
	scale_x_continuous(breaks = c(0, 1, 1.5, 2)) +
	theme(legend.position = "bottom", legend.box = "vertical") +
	labs(lty = "True Richness", colour = "Richness Estimator", x = "Evenness", y = "Coverage-standardised Richness") +
	facet_grid(.~quorum, scales = "free")
graphics.off()


pdf("./simulation-figures/Fig. S<qf-trueD-vs-qD>.pdf", width = simfigwidth * 1.2, height = simfigheight)
qf %>% filter(quorum %in% c(0.1,0.5,0.9,0.99)) %>%
	ggplot(aes(x = true_richness, y = qD, colour = as.factor(sd))) +
	geom_point() +
	geom_line() +
	scale_y_continuous(trans = "log2", breaks = c(3,6,12,25,50,100,200,400)) +
	scale_x_continuous(trans = "log2", breaks = c(50,100,200,400)) +
	scale_color_brewer(palette = "Set1") +
	labs(x = "True Richness", y = "Coverage-standardised Richness", colour = "Evenness (SD)") +
	theme(legend.position = "bottom") +
	facet_wrap(method~quorum, scales = "free")
graphics.off()

pdf("./simulation-figures/Fig. S<sd-vs-qD-via-quorum>.pdf", height = 3.5, width = 9)
simdf %>% filter(method %in% c("sqs_alroy") & quorum %in% c(0.100,0.500,0.900,0.990) & sample_size == 10000) %>%
	transform(method = factor(method, labels = c("SQS"))) %>%
	ggplot(aes(x = sd, y = qD, colour = as.factor(true_richness))) +
	geom_line() +
	scale_y_log10() +
	theme(legend.position = "bottom") +
	labs(x = "Evenness (SD)", y = "Coverage-standardised Richness", colour = "True Richness") +
	facet_grid(method~quorum, scales = "free", labeller = labeller(method = as_labeller(method_names)))
graphics.off()

pdf("./simulation-figures/Fig. S<rf-sd-vs-qD>.pdf", width = simfigwidth * 1.2, height = simfigheight)
plotting_sizes <- sizes[seq(from = 15, to = (length(sizes) - 6), length.out = 4)]
simdf %>%
	filter(sample_size %in% plotting_sizes, method %in% c("raw","chao","lambda"), quorum %in% c(NA,0.990), quota %in% c(NA,800)) %>%
	ggplot(aes(x = sd, y = qD, colour = as.factor(true_richness))) +
	geom_point() +
	geom_line() +
	scale_y_continuous(trans = "log2") +
	theme(legend.position = "bottom") +
	labs(x = "Evenness (SD)", y = "Size-standardised Richness", colour = "True Richness") +
	facet_wrap(method~sample_size, scales = "free")
graphics.off()

pdf("./simulation-figures/Fig. S<rf-sd-vs-qD-superimposed>.pdf", width = 5, height = 4)
plotting_sizes <- sizes[seq(from = length(sizes) - 15, to = (length(sizes) - 9), length.out = 2)]
simdf %>%
	filter(sample_size %in% plotting_sizes, method %in% c("raw","chao","lambda"), quorum %in% c(NA,0.990), quota %in% c(NA,800)) %>%
	filter(true_richness %in% c(200,400)) %>%
	mutate(method = recode(method, "raw" = "Face-value Counts", "chao" = "Chao1", "lambda" = "Lambda-5")) %>%
	mutate(sample_size = paste("Quota =", sample_size)) %>% mutate(sample_size = factor(sample_size, levels = sort(unique(sample_size), decreasing = T))) %>%
	ggplot(aes(x = sd, y = qD, colour = as.factor(method), lty = as.factor(true_richness))) +
	geom_point() +
	geom_line() +
	scale_y_continuous(trans = "log2", breaks = c(1,2,5,10,25,50,100,200,400,800)) +
	scale_x_continuous(breaks = c(0, 1, 1.5, 2)) +
	theme(legend.position = "bottom", legend.box = "vertical") +
	labs(lty = "True Richness", colour = "Richness Estimator", x = "Evenness", y = "Size-standardised Richness") +
	facet_grid(.~sample_size, scales = "free")
graphics.off()

pdf("./simulation-figures/Fig. S<rf-trueD-vs-qD>.pdf", width = simfigwidth * 1.2, height = simfigheight)
plotting_sizes <- sizes[seq(from = 15, to = (length(sizes) - 6), length.out = 4)]
simdf %>% filter(sample_size %in% plotting_sizes, method %in% c("raw","chao","lambda"), quorum %in% c(NA,0.990), quota %in% c(NA,800)) %>%
	ggplot(aes(x = true_richness, y = qD, colour = as.factor(sd_labels))) +
	geom_point() +
	geom_line() +
	scale_y_continuous(trans = "log2", breaks = c(3,6,12,25,50,100,200,400,800)) +
	scale_x_continuous(trans = "log2", breaks = c(50,100,200,400,800)) +
	scale_color_brewer(palette = "Set1") +
	labs(x = "True Richness", y = "Size-standardised Richness", colour = "Evenness") +
	theme(legend.position = "bottom") +
	facet_wrap(method_labels~sample_size, scales = "free")
graphics.off()
