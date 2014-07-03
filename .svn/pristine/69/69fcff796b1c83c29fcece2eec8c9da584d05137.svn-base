require(SDMTools)
require(rgdal)
require(maptools)
require(raster)

gpclibPermit()

source(paste(src.dir,"/000.zipWrite.R",sep=""))

#idir <- "F:/gap_analysis_publications/gap_phaseolus/modeling_data"
#idir <- "G:/ncastaneda/gap-analysis-tomato/gap_tomato/occurrences/test"
#outName <- "gsamples-buffer.asc"
#spID <- "Solanum_huaylasense"
#outFolder <- paste(idir, "/", sep="")
#outFolder <- paste("F:/gap_analysis_publications/gap_phaseolus/samples_calculations/", spID, sep="")
#buffDist <- 50000
#spFile <- paste(idir, "/occurrence_files/", spID, ".csv", sep="")
#spFile <- paste(idir, "/", spID, ".csv", sep="")
#msk <- "G:/ncastaneda/clim/bio_10m_esri/bio_1"

createBuffers <- function(spFile, outFolder, outName, buffDist, msk) {
	occ <- read.csv(spFile)
	
	cat("Buffering the points \n")
	
	pb <- pbCreate(nrow(occ), "text", style=3)
	for (pnt in 1:nrow(occ)) {
		#cat("Point", pnt, "out of", nrow(occ), "\n")
		
		lons = lats = NULL #objects to store points buffering hull points
		for (bear in seq(0,length=360,by=1)) { #cycle through 360 directions and get points at 50 km from hull points
			tt = destination(occ$lat[pnt],occ$lon[pnt],bearing = bear,distance = buffDist)
			lons = c(lons,tt$lon2); lats = c(lats,tt$lat2)
		}
		
		buff <- cbind(lons, lats)
		buff <- rbind(buff, buff[1,])
		
		assign(paste("pol",pnt, sep=""), Polygons(list(Polygon(buff)), pnt))
		
		if (pnt == 1) {
			polgrp <- c(get(paste("pol",pnt,sep="")))
		} else {
			polgrp <- c(polgrp, get(paste("pol",pnt,sep="")))
		}
		pbStep(pb, pnt)
	}
	pbClose(pb)
	
	cat("Creating the raster \n")
	msk <- raster(msk)
	pa <- rasterize(SpatialPolygons(polgrp), msk)

	pa[which(!is.na(pa[]))] <- 1
	pa[which(is.na(pa[]) & msk[] == 1)] <- 0
	pa[which(is.na(msk[]))] <- NA

	cat("Writing output \n")
	paName <- zipWrite(pa, outFolder, paste(outName, ".gz", sep=""))

	return(pa)
}
