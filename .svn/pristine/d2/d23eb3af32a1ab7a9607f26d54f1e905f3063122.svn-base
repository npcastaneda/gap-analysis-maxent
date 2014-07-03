# Distance to poins - parallelized
# H. Achicanoy
# Jul - 2013

source(paste(src.dir,"/000.zipWrite.R",sep=""))

setOptions(overwrite=T)

populationDistance <- function(bdir,spID) {
  idir <- paste(bdir, "/maxent_modeling", sep="")
  odir <- paste(bdir, "/samples_calculations", sep="")
  spOutFolder <- paste(odir,"/",spID, sep="")
  
  cat("Loading occurrences \n")
  occ <- read.csv(paste(bdir, "/occurrence_files/",spID, ".csv", sep=""))
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

ParProcess <- function(bdir, ncpu) {
	
	spList <- list.files(paste(bdir, "/occurrence_files", sep=""),pattern=".csv")
	
	pDist_wrapper <- function(i) {
    library(raster)
    library(rgdal)
    sp <- spList[i]
		sp <- unlist(strsplit(sp, ".", fixed=T))[1]
		cat("\n")
		cat("...Species", sp, "\n")
		out <- summarizeDistances(bdir)
	}
  
  library(snowfall)
  sfInit(parallel=T,cpus=ncpu)
  
  sfExport("populationDistance")
  sfExport("summarizeDistances")
  sfExport("zipWrite")
  sfExport("bdir")
  sfExport("src.dir")
  sfExport("spList")
  
  #run the control function
  system.time(sfSapply(as.vector(1:length(spList)), pDist_wrapper))
  
  #stop the cluster
  sfStop()
  
  return("Done!")
  
}
