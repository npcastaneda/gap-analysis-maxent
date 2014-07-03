require(rgdal)
require(raster)

source(paste(src.dir,"/000.zipWrite.R",sep=""))

setOptions(overwrite=T)

#bdir <- "F:/gap_analysis_publications/gap_phaseolus"
#spID <- "Phaseolus_acutifolius"

populationDistance <- function(bdir, spID) {
  idir <- paste(bdir, "/maxent_modeling", sep="")
  odir <- paste(bdir, "/samples_calculations", sep="")
  spOutFolder <- paste(odir, "/", spID, sep="")
  
  cat("Loading occurrences \n")
  occ <- read.csv(paste(bdir, "/occurrence_files/", spID, ".csv", sep=""))
  xy <- occ[,2:3]
  
  cat("Loading mask \n")
  msk <- raster(paste(bdir, "/masks/mask.asc", sep=""))
  
  cat("Distance from points \n")
  dgrid <- distanceFromPoints(msk, xy)
  dgrid[which(is.na(msk[]))] <- NA
  
  cat("Writing output \n")
  dumm <- zipWrite(dgrid, spOutFolder, "pop-dist.asc.gz")
  return(dgrid)
}

summarizeDistances <- function(bdir) {
  spList <- list.files(paste(bdir, "/occurrence_files", sep=""))
  sppC <- 1
  
  for (spp in spList) {
    
    spp <- unlist(strsplit(spp, ".", fixed=T))[1]
    pdir <- paste(bdir,"/samples_calculations/",spp,sep="")
    #pop=sum(pdir=="pop-dist.asc.gz")
    if(file.exists(paste(pdir,"/pop-dist.asc.gz",sep=""))){
      cat("The file already exists",spp,"\n")
    }else{
      cat("Processing taxon", spp, "\n")
      if(!file.exists(pdir)){dir.create(pdir)}
      dg <- populationDistance(bdir, spp)}
  }
  return(spList)
}


