clades <- c(
	"dinosauromorpha"
	,"tetrapoda"
)
for (cld in 1:length(clades)) {
	clade <- clades[cld]
	save(clades, clade, cld, file = "./saved_settings.RData")
	source("./DTTTT.R", echo = TRUE)
	load("./saved_settings.RData")
}
