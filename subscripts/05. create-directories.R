#### create directory structure ####

folder.name <- "./output"
if (new.output.folder == TRUE) {
	folder.name <- paste(folder.name, "/", str_replace_all(str_replace_all(Sys.time(), fixed(":"), "-"), fixed(" "), "-"), sep = "")
}

if (clear_plots == TRUE) {
	do.call(file.remove,list(list.files(paste(folder.name, "/plots/", clade, "/", sep = ""), full.names = TRUE, recursive = T)))
}

for (rgn in 1:length(regions)) {
	dir.create(file.path(folder.name, "plots", clade, names(regions)[rgn]), recursive = T, showWarnings = F)
}

dir.create(file.path(folder.name, "plots", clade, "global"), recursive = T, showWarnings = F)
dir.create(file.path(folder.name, "plots", clade, "pan-regional"), recursive = T, showWarnings = F)

