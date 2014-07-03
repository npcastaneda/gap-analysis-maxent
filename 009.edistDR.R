require(rgdal)
require(raster)

source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.bufferPoints.R",sep=""))

#Calculate the size of the DR, of the convexhull in km2, of the native area, and of the herbarium samples
#based on the area of the cells

edistDR <- function(bdir, spID) {
	
	idir <- paste(bdir, "/maxent_modeling", sep="")
	ddir <- paste(bdir, "/samples_calculations", sep="")
	pcdir <- paste(bdir, "/biomod_modeling/current-clim/wwf_eco", sep="")
	
	#Creating the directories
	if (!file.exists(ddir)) {
		dir.create(ddir)
	}
	
	spOutFolder <- paste(ddir, "/", spID, sep="")
	if (!file.exists(spOutFolder)) {
		dir.create(spOutFolder)
	}
	
	#Read the thresholded raster (PA), multiply it by the area raster and sum up those cells that are != 0
	cat("Taxon", spID, "\n")
	spFolder <- paste(idir, "/models/", spID, sep="")
	projFolder <- paste(spFolder, "/projections", sep="")
	
	cat("Loading wwf terrestrial ecosystems \n")
	pc1 <- raster(paste(pcdir,"/wwf_eco_terr.asc",sep=""))
# 	rs <- raster(xmn=-179.125, xmx=179.75, ymn=-56, ymx=83.62501, ncols=8613, nrows=3351)
# 	pc1 <- setExtent(pc1, extent(rs), keepres=FALSE, snap=FALSE)
#   rm(rs)
 
	#Edist of the DR
	cat("Reading presence/absence surface \n")
	grd <- paste(spID, "_worldclim2_5_EMN_PA.asc.gz", sep="")
	spList <- read.csv(paste(bdir, "/summary-files/taxaForRichness.csv", sep=""))
	isValid <- spList$IS_VALID[which(spList$TAXON == paste(spID))]
	
	if (isValid == 1) {
		grd <- zipRead(projFolder, grd)
		
		cat("Env. distribution of the DR \n")
		grda <- grd * pc1 #PC1
		edistDR1 <- unique(grda[])
		edistDR1 <- edistDR1[which(edistDR1 != 0 & !is.na(edistDR1))]
		rm(grda)
				
		rm(grd)
	} else {
		edistDR1 <- NULL
	}
	
	#Edist of the convex-hull
	if (file.exists(paste(ddir, "/", spID, "/convex-hull.asc.gz", sep=""))) {
		cat("Reading convex hull \n")
		grd <- zipRead(paste(ddir, "/", spID, sep=""), "convex-hull.asc.gz")
		
		cat("Env. distribution of the convex hull \n")
		grda <- grd * pc1
		edistCH1 <- unique(grda[])
		edistCH1 <- edistCH1[which(edistCH1 != 0 & !is.na(edistCH1))]
		rm(grda)
		
		rm(grd)
	} else {
		edistCH1 <- NULL
# 		edistCH2 <- NULL
	}
  
  
	#Edist of the native area
	naFolder <- paste(bdir, "/biomod_modeling/native-areas/asciigrids/", spID, sep="")
	if (file.exists(paste(naFolder, "/narea.asc.gz", sep=""))) {
		cat("Reading native area \n")
		grd <- zipRead(naFolder, "narea.asc.gz")
		
		cat("Env. distribution of the native area \n")
		grda <- grd * pc1
		edistNA1 <- unique(grda[])
		edistNA1 <- edistNA1[which(edistNA1 != 0 & !is.na(edistNA1))]
		rm(grda)
		
		rm(grd)
	} else {
		edistNA1 <- NULL
	}
	
	#Edist of the herbarium samples CA50
	cat("Reading h-samples buffer \n")
	if (file.exists(paste(ddir, "/", spID, "/hsamples-buffer.asc.gz", sep=""))) {
		grd <- zipRead(paste(ddir, "/", spID, sep=""), "hsamples-buffer.asc.gz")
		
		cat("Env. distribution of h-samples buffer \n")
		grda <- grd * pc1
		edistHB1 <- unique(grda[])
		edistHB1 <- edistHB1[which(edistHB1 != 0 & !is.na(edistHB1))]
		rm(grda)
		
		rm(grd)
	} else {
		edistHB1 <- NULL
	}
	
	#Size of the germplasm samples CA50
	cat("Reading g-samples buffer \n")
	if (file.exists(paste(ddir, "/", spID, "/gsamples-buffer.asc.gz", sep=""))) {
		grd <- zipRead(paste(ddir, "/", spID, sep=""), "gsamples-buffer.asc.gz")
		
		cat("Env. distribution of g-samples buffer \n")
		grda <- grd * pc1
		edistGB1 <- unique(grda[])
		edistGB1 <- edistGB1[which(edistGB1 != 0 & !is.na(edistGB1))]
		rm(grda)
		
		rm(grd)
	} else {
		edistGB1 <- NULL
	}
	#Writing results
	outDF <- data.frame(DRDist.PC1=length(edistDR1), CHDist.PC1=length(edistCH1), NADist.PC1=length(edistNA1), HBDist.PC1=length(edistHB1), GBDist.PC1=length(edistGB1))
  
	write.csv(outDF, paste(spOutFolder, "/edist_wwf.csv", sep=""), quote=F, row.names=F)
	return(outDF)
}


summarizeDR_env <- function(idir) {
	
	ddir <- paste(idir, "/samples_calculations", sep="")
	
	odir <- paste(idir, "/maxent_modeling/summary-files", sep="")
	if (!file.exists(odir)) {
		dir.create(odir)
	}
	
	spList <- list.files(paste(idir, "/occurrence_files", sep=""))
	
	sppC <- 1
	for (spp in spList) {
		spp <- unlist(strsplit(spp, ".", fixed=T))[1]
		fdName <- spp #paste("sp-", spp, sep="")
		spFolder <- paste(idir, "/maxent_modeling/models/", fdName, sep="")
		spOutFolder <- paste(ddir, "/", spp, sep="")
		
		if (file.exists(spFolder)) {
			
			res <- edistDR(idir, spp)
			
			metFile <- paste(spOutFolder, "/edist_wwf.csv", sep="")
			metrics <- read.csv(metFile)
			metrics <- cbind(taxon=spp, metrics)
			
			if (sppC == 1) {
				outSum <- metrics
			} else {
				outSum <- rbind(outSum, metrics)
			}
			sppC <- sppC + 1
		} else {
			cat("The taxon was never modeled \n")
		}
	}
	
	outFile <- paste(odir, "/edist_wwf.csv", sep="")
	write.csv(outSum, outFile, quote=F, row.names=F)
	return(outSum)
}
