require(rgdal)
require(raster)

source(paste(src.dir,"/000.zipRead.R",sep=""))

# Script to calculate proportion of the dist. range with SD above 0.15 (ASD15)

# idir <- crop_dir
# spID <- "Avena_abyssinica"

calcASD15 <- function(idir, spID) {
	cat("Taxon", spID, "\n")
	spFolder <- paste(idir, "/maxent_modeling/models/", spID, sep="")
	projFolder <- paste(spFolder, "/projections", sep="")
	
	esdCpt <- paste(spID, "_worldclim2_5_ESD.asc.gz", sep="")
	esdThr <- paste(spID, "_worldclim2_5_ESD_PR.asc.gz", sep="")
	
	dumm <- paste(spID, "_worldclim2_5_EMN.asc.gz", sep="")
	
  if (file.exists(paste(projFolder,"/",dumm,sep=""))) {
  	cat("..Reading raster files \n")
  	dumm <- zipRead(projFolder, dumm)
  	esdCpt <- zipRead(projFolder, esdCpt)
  	esdThr <- zipRead(projFolder, esdThr)
  	
  	esdCpt[which(dumm[] < 0.001)] <- NA
  	
  	rm(dumm)
  	
  	esdThr[which(esdThr[] == 0)] <- NA
  	
  	cat("..Calculating \n")
  	szCpt <- length(which(esdCpt[] >= 0))
  	szCptUncertain <- length(which(esdCpt[] >= 0.15))
  	rateCpt <- szCptUncertain / szCpt * 100
  	
  	szThr <- length(which(esdThr[] >= 0))
  	szThrUncertain <- length(which(esdThr[] >= 0.15))
  	rateThr <- szThrUncertain / szThr * 100
  	
  	cat("..Writing results \n")
  	dfOut <- data.frame(taxon=spID, sizeComplete=szCpt, sizeCompleteUncertain=szCptUncertain, rateComplete=rateCpt, sizeThresholded=szThr, sizeThresholdedUncertain=szThrUncertain, rateThresholded=rateThr)
  } else {
    cat("..Writing results \n")
    dfOut <- data.frame(taxon=spID, sizeComplete=NA, sizeCompleteUncertain=NA, rateComplete=NA, sizeThresholded=NA, sizeThresholdedUncertain=NA, rateThresholded=NA)
  }
  	
	oFile <- paste(spFolder, "/metrics/ASD15.csv", sep="")
	write.csv(dfOut, oFile, quote=F, row.names=F)
	
	return(dfOut)
}

summarizeASD15 <- function(idir) {
	
	odir <- paste(idir, "/summary-files", sep="")
	if (!file.exists(odir)) {
		dir.create(odir)
	}
	
	spList <- list.files(paste(idir, "/occurrence_files", sep=""))
	
	sppC <- 1
	for (spp in spList) {
		spp <- unlist(strsplit(spp, ".", fixed=T))[1]
		fdName <- spp #paste("sp-", spp, sep="")
		spFolder <- paste(idir, "/maxent_modeling/models/", fdName, sep="")
		
		if (file.exists(spFolder)) {
			
			res <- calcASD15(idir, spp)
			
			metFile <- paste(spFolder, "/metrics/ASD15.csv", sep="")
			metrics <- read.csv(metFile)
			
			if (sppC == 1) {
				outSum <- metrics
			} else {
				outSum <- rbind(outSum, metrics)
			}
			sppC <- sppC + 1
		}
	}
	outFile <- paste(idir, "/maxent_modeling/summary-files/ASD15.csv", sep="")
	write.csv(outSum, outFile, quote=F, row.names=F)
	
}
