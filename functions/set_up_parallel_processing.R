# set up parallel processing depending on architecture
library(foreach)
if ( Sys.info()[['sysname']] == 'Darwin' | Sys.info()[['sysname']] == 'Linux' ) {   # Mac OS X or Linux
	library(doMC)
	cores <- detectCores()
	registerDoMC(cores)
}

if ( Sys.info()[['sysname']] == 'Windows' ) {
	# Windows
	library(doParallel)
	library(parallel)
	cores <- detectCores()
	cl <- makeCluster(cores - 1, outfile = "")
	registerDoParallel(cl)
}
