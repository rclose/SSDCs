#### plot results of Simulations 1 ####

pdf("./simulation-figures/Fig. S<rarefied-coverage>.pdf", width = 6, height = 5)
simdf %>% filter(method == "raw") %>%
	filter(true_richness %in% c(400)) %>%
	ggplot(aes(x = sample_size, y = 1 - SC, ymax = 1 - SC.UCL, ymin = 1 - SC.LCL, colour = as.factor(sd_labels))) +
	geom_line() +
	geom_ribbon(aes(fill = as.factor(sd_labels)), alpha = 0.1, colour = NA, colour = NA) +
	# theme_bw() +
	theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal") +
	scale_x_log10() +
	scale_y_log10() +
	annotation_logticks(side = "lb") +
	labs(colour = "Evenness", fill = "Evenness", y = "1 - Good's u (=Coverage Deficit)", x = "Sample Size")
	# facet_grid(true_richness~.)
graphics.off()

#extrapolator methods
pdf("./simulation-figures/Fig. S<rarefied_extrapolators>.pdf", width = simfigwidth, height = simfigheight)
simdf %>% filter(method %in% c("chao","lambda","raw","trips")) %>%
	ggplot(aes(x = sample_size, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = as.factor(true_richness), fill = as.factor(true_richness))) +
	geom_line() +
	geom_ribbon(alpha = 0.2, colour = NA) +
	geom_line(aes(x = sample_size, y = true_richness, colour = as.factor(true_richness)), lty = 2) +
	# theme_bw() +
	theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal") +
	scale_x_log10() +
	scale_y_log10() +
	xlab("Sample Size") +
	ylab("qD") +
	# scale_y_continuous(trans='log2', breaks=c(1,2,4,8,16)) +
	labs(colour = "True Richness", fill = "True Richness", lty = "Evenness", y = "Size-standardised Richness", x = "Sample Size") +
	facet_grid(method_labels~sd_labels)
graphics.off()



#richness rarefaction curves
pdf("./simulation-figures/Fig. S<rarefied-estimators>.pdf", width = 9, height = 1.8*7) #uncomment and comment bits to make version with CIs
# plotting_richnesses <- c(50,100,200,400)
plotting_richnesses <- c(50,200)
plotting_sds <- c(0,1,2)
voffset <- 0.65
i <- 20; lab1 <- ddply(filter(simdf, sample_size == sizes[i], quorum %in% c(NA, plotting_quorum), quota %in% c(NA,plotting_quota), !grepl("sqs_alroy|q1|q2", method), true_richness %in% plotting_richnesses, sd %in% plotting_sds),.(method_labels,true_richness,sd_labels), summarise, qD.ratio = round(qD.ratio, 1), qD = mean(qD)); lab1$sample_size <- sizes[i]; lab1$qD.ratio <- as.character(round(lab1$qD.ratio, 1)); lab1$qD.ratio[is.na(lab1$qD.ratio)] = "NA"
i <- 30; lab2 <- ddply(filter(simdf, sample_size == sizes[i], quorum %in% c(NA, plotting_quorum), quota %in% c(NA,plotting_quota), !grepl("sqs_alroy|q1|q2", method), true_richness %in% plotting_richnesses, sd %in% plotting_sds),.(method_labels,true_richness,sd_labels), summarise, qD.ratio = round(qD.ratio, 1), qD = mean(qD)); lab2$sample_size <- sizes[i]; lab2$qD.ratio <- as.character(round(lab2$qD.ratio, 1)); lab2$qD.ratio[is.na(lab2$qD.ratio)] = "NA"
i <- 45; lab3 <- ddply(filter(simdf, sample_size == sizes[i], quorum %in% c(NA, plotting_quorum), quota %in% c(NA,plotting_quota), !grepl("sqs_alroy|q1|q2", method), true_richness %in% plotting_richnesses, sd %in% plotting_sds),.(method_labels,true_richness,sd_labels), summarise, qD.ratio = round(qD.ratio, 1), qD = mean(qD)); lab3$sample_size <- sizes[i]; lab3$qD.ratio <- as.character(round(lab3$qD.ratio, 1)); lab3$qD.ratio[is.na(lab3$qD.ratio)] = "NA"
p <- simdf %>%
	filter(quorum %in% c(NA, plotting_quorum), quota %in% c(NA,plotting_quota), !grepl("sqs_alroy|q1|q2", method), true_richness %in% plotting_richnesses, sd %in% plotting_sds) %>%
	# mutate(tmp_labels = recode(method_labels, 'SQS (-dominant)' = paste("SQS", plotting_quorum), 'SQS q=0' = paste("SQS", plotting_quorum), 'CR' = paste("CR", plotting_quota))) %>%
	ggplot(aes(x = sample_size, y = qD, colour = as.factor(true_richness), fill = as.factor(true_richness))) +
	geom_line() +
	# geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.2, colour = NA) +
	geom_point(data = lab1, aes(x = sample_size, y = qD, group = NULL)) +
	geom_point(data = lab2, aes(x = sample_size, y = qD, group = NULL)) +
	geom_point(data = lab3, aes(x = sample_size, y = qD, group = NULL)) +
	scale_x_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
	scale_y_continuous(trans = 'log2', breaks = c(1,2,5,10,25,50,100,200,400,800,1600)) +
	geom_line(aes(x = sample_size, y = true_richness, colour = as.factor(true_richness)), lty = 3, cex = 0.5) +
	xlab(label = "Sample Size") +
	ylab(label = "Richness") +
	theme_bw() +
	theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal") +
	labs(colour = "True Richness", fill = "True Richness", y = "Size-standardised Richness") +
	facet_grid(method_labels~sd_labels)
p +
	geom_text(data = lab1, aes(x = sample_size, y = qD * voffset, group = NULL, label = paste("(",qD.ratio,")", sep = "")), cex = 3) +
	geom_text(data = lab2, aes(x = sample_size, y = qD * voffset, group = NULL, label = paste("(",qD.ratio,")", sep = "")), cex = 3) +
	geom_text(data = lab3, aes(x = sample_size, y = qD * voffset, group = NULL, label = paste("(",qD.ratio,")", sep = "")), cex = 3)
# p
graphics.off()

pdf("./simulation-figures/Fig. <rarefied-estimators-superimposed>.pdf", width = 8, height = 1.6*2)
# plotting_richnesses <- c(50,100,200,400)
plotting_richnesses <- c(400)
p <- simdf %>%
	filter(quorum %in% c(NA, plotting_quorum), quota %in% c(NA,plotting_quota), !grepl("sqs|inext", method), true_richness %in% plotting_richnesses) %>%
	# filter(quorum %in% c(NA, plotting_quorum), quota %in% c(NA,plotting_quota), true_richness %in% c(100,400), !grepl("sqs|inext", method)) %>%
	mutate(tmp_labels = recode(method_labels, 'SQS (-dominant)' = paste("SQS", plotting_quorum), 'CR' = paste("CR", plotting_quota))) %>%
	ggplot(aes(x = sample_size, y = qD, colour = as.factor(method_labels), fill = as.factor(method_labels))) +
	geom_line() +
	# geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.2, colour = NA) +
	# geom_point(data = lab1, aes(x = sample_size, y = qD, group = NULL)) +
	# geom_point(data = lab2, aes(x = sample_size, y = qD, group = NULL)) +
	# geom_point(data = lab3, aes(x = sample_size, y = qD, group = NULL)) +
	scale_x_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
	scale_y_continuous(trans = 'log2', breaks = c(1,5,10,25,50,100,200,400,800,1600)) +
	geom_line(aes(x = sample_size, y = true_richness), lty = 3, cex = 0.5) +
	xlab(label = "Sample Size") +
	ylab(label = "Richness") +
	# theme_bw() +
	theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal") +
	labs(colour = "Richness Estimator", y = "Size-standardised Richness", x = "Sample Size") +
	facet_grid(true_richness~sd_labels)
p
graphics.off()



pdf("./simulation-figures/Fig. S<SQS-sd-qD>.pdf", width = 10, height = 4.5)
simdf %>%
	filter(method %in% coverage_vars & !grepl("sqs_inext", method)) %>%
	filter(true_richness %in% c(50,100,200)) %>%
	filter(quorum %in% c(0.1, 0.6, 0.95)) %>%
	mutate(quorum = paste("Quorum =", quorum)) %>% mutate(quorum = factor(quorum, levels = sort(unique(quorum), decreasing = F))) %>%
	ggplot(aes(x = sample_size, y = qD.ratio, ymin = qD.LCL.ratio, ymax = qD.UCL.ratio, lty = as.factor(sd_labels), colour = as.factor(true_richness), fill = as.factor(true_richness))) +
	geom_line(aes(x = sample_size, y = true_richness.ratio), lty = 1, colour = "darkgrey") +
	geom_line() +
	geom_ribbon(alpha = 0.2, colour = NA) +
	theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "horizontal") +
	scale_x_log10() +
	scale_y_continuous(trans = 'log2', breaks = c(1,2,4,8,16)) +
	# geom_line(aes(x = sample_size, y = true_richness.ratio, colour = as.factor(true_richness)), lty = 2) +
	# ggtitle(coverage_vars[i]) +
	# geom_text(data = filter(simdf, sample_size == maxss, method == coverage_vars[i]), aes(x = sample_size + 1000, y = qD.ratio, label = round(qD.ratio, 1)), inherit.aes = F, nudge_x = 0.3, cex = 2.5) + #+ geom_dl(aes(label = assemblage), method = list(dl.combine("last.points"), cex = 0.8))
	# theme_bw() +
	labs(colour = "True Richness", fill = "True Richness", lty = "Evenness", y = "Richness Ratio", x = "Sample Size") +
	facet_wrap(~quorum)
graphics.off()

pdf("./simulation-figures/Fig. S<n-ton-curves>.pdf", width = simfigwidth, height = simfigheight * 0.6)
simdf %>%
	filter(true_richness %in% c(50,400)) %>%
	filter(method == "raw") %>%
	rename(richness = qD) %>%
	mutate(sd = as.factor(sd)) %>% mutate(sd = recode(sd, "0" = "flat", "1" = "sigma = 1", "1.5" = "sigma = 1.5", "2" = "sigma = 2")) %>%
	mutate(true_richness = paste("True Richness =", true_richness)) %>% mutate(true_richness = factor(true_richness, levels = sort(unique(true_richness), decreasing = T))) %>%
	gather("ton_type","ton_counts", c(f1,f2,f3,fm,richness)) %>%
	ggplot(aes(x = sample_size, y = ton_counts, colour = ton_type)) +
	geom_line() +
	scale_y_log10(breaks = c(1,5,10,25,50,100,200,400,800,1600)) +
	scale_x_log10() +
	# theme_bw() +
	theme(legend.position = "bottom") +
	labs(colour = "Variable") +
	facet_grid(true_richness ~ as.factor(sd)) +
	ylab("Count") +
	xlab("Sample Size")
graphics.off()


