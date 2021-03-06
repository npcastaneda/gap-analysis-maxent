require(rgdal)
require(raster)

source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.bufferPoints.R",sep=""))

#bdir <- "F:/gap_analysis_publications/gap_phaseolus"
#When presence surface not available, get all the populations and use the 50km buffer

speciesRichness <- function(bdir) {
	idir <- paste(bdir, "/biomod_modeling", sep="")
	ddir <- paste(bdir, "/samples_calculations", sep="")
	
	outFolder <- paste(bdir, "/species_richness", sep="")
	if (!file.exists(outFolder)) {
		dir.create(outFolder)
	}
	
	spList <- read.csv(paste(idir, "/summary-files/taxaForRichness.csv", sep=""))
	allOcc <- read.csv(paste(bdir, "/occurrences/",crop,".csv", sep=""))
	
	sppC <- 1
	rcounter <- 1
	scounter <- 1
	for (spp in spList$Taxon) {
		
		cat("Processing taxon", spp, "\n")
		
		isValid <- spList$ValidModel[spList$Taxon == paste(spp)]
		
		if (isValid == 1) {
			exOcc <- T
			cat("Load presence/absence raster and sd raster \n")
			sppFolder <- paste(idir, "/models/proj.", spp, sep="")
			projFolder <- paste(sppFolder, "/grdfiles", sep="")
			
			pagrid <- paste(spp, "_PA_NA.asc", sep="")
			#pagrid <- zipRead(projFolder, pagrid)
			pagrid <- raster(paste(projFolder,"/",pagrid,sep=""))
			
			assign(paste("sdgrid",sppC,sep=""), pagrid)
      #assign(paste("sdgrid",sppC,sep=""), paste(spp, "_PA_NA.asc", sep=""))
			#assign(paste("sdgrid",sppC,sep=""), zipRead(projFolder, get(paste("sdgrid",sppC,sep=""))))
      
		} else {
			cat("Samples buffer \n")
			tallOcc <- allOcc[which(allOcc$Taxon == paste(spp)),]
			if (nrow(tallOcc) != 0) {
				exOcc <- T
				spOutFolder <- paste(ddir, "/", spp, sep="")
				if (!file.exists(spOutFolder)) {
					dir.create(spOutFolder)
				}
				tallOcc <- as.data.frame(cbind(as.character(tallOcc$Taxon), tallOcc$lon, tallOcc$lat))
				names(tallOcc) <- c("taxon", "lon", "lat")
				
				write.csv(tallOcc, paste(spOutFolder, "/samples.csv", sep=""), quote=F, row.names=F)
				rm(tallOcc)
				pagrid <- createBuffers(paste(spOutFolder, "/samples.csv", sep=""), spOutFolder, "samples-buffer.asc", 50000, paste(bdir, "/masks/mask.asc", sep=""))
			} else {
				exOcc <- F
			}
		}
		
		#Species richness
		if (rcounter == 1) {
			if (isValid == 1 | exOcc) {
				sprich <- pagrid
				rm(pagrid)
				rcounter <- rcounter + 1
			}
		} else {
			if (isValid == 1 | exOcc) {
				sprich <- sprich + pagrid
				rm(pagrid)
				rcounter <- rcounter + 1
			}
		}
		
		#Species sd
		if (scounter == 1) {
			if (isValid == 1) {
				sdlist <- get(paste("sdgrid",sppC,sep=""))
				scounter <- scounter + 1
			}
		} else {
			if (isValid == 1) {
				sdlist <- c(sdlist, get(paste("sdgrid",sppC,sep="")))
				scounter <- scounter + 1
			}
		}
		
		sppC <- sppC + 1
	}

	cat("Writing richness raster \n")
	dumm <- zipWrite(sprich, outFolder, "species-richness.asc.gz")
	
	cat("Calculating mean sd raster \n")
	sdmean <- sum(stack(sdlist)) / sprich
	cat("Writing \n")
	dumm <- zipWrite(sdmean, outFolder, "species-richness-sdmean.asc.gz")

	cat("Calculating max sd raster \n")
	sdmax <- max(stack(sdlist))
	cat("Writing \n")
	dumm <- zipWrite(sdmax, outFolder, "species-richness-sdmax.asc.gz")

	cat("Done! \n")

	return(stack(sprich,sdmean,sdmax))
}
