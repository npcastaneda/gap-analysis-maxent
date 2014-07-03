#Take the list of species to model and create a new file with occurrences of selected species and climate data

createOccFiles <- function(occ, taxfield, outDir) {
	
	cat("\n")
	
	if (!file.exists(outDir)) {
		dir.create(outDir)
	}
	
	cat("Loading occurrence file \n")
	
	occ <- read.csv(occ)
	
	spL <- unique(occ[,taxfield])
	
	spcounter <- 1
	
	cat("Now printing \n")
	
	for (sp in spL) {
		
		nspp <- length(spL)
		cat("      ...", paste(round(spcounter/nspp*100,2),"% Completed", sep=""), "\n")
		
		spData <- occ[which(occ[,taxfield] == sp),]
		
		csvName <- paste(outDir, "/", sp, ".csv", sep="")
		write.csv(spData, csvName, row.names=F, quote=F)
		rm(spData)
		
		spcounter <- spcounter+1
	}
	return("Done")
}
