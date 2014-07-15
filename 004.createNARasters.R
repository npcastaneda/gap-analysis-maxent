require(rgdal)
require(raster)
require(maptools)
gpclibPermit(); options(warn=-1)

source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.createChullBuffer.R",sep=""))

createNARaster <- function(spID,inDir) {
  
  cat("\n")
  cat("Taxon", spID ,"\n")
  
  inNADir <- paste(inDir, "/native-areas/polyshps", sep="")
  outNADir <- paste(inDir, "/native-areas/asciigrids", sep="")
  outFolder <- paste(outNADir, "/", spID, sep="")
  allOcc <- read.csv(paste(crop_dir, "/occurrences/",crop,".csv", sep=""))
  
  if (!file.exists(paste(outFolder, "/narea.asc.gz", sep=""))) {
    
    cat("Not processed, thus processing \n")
    
    if (!file.exists(outNADir)) {
      dir.create(outNADir)
    }
    
    if (!file.exists(outFolder)) {
      dir.create(outFolder)
    }
    
    shpName <- paste(inNADir, "/", spID, "/narea.shp", sep="")
    
    if(!file.exists(shpName)){
      #Creating convex hull instead
      cat("No shapefile available, thus creating convex hull \n")
      occFile <- paste(crop_dir, "/occurrence_files/", spID, ".csv", sep="")
      tallOcc <- allOcc[which(allOcc$Taxon == paste(spID)),]
      if(nrow(tallOcc)==0){
        cat("No coordinates for this species /n")
      }else{
        NAGrid <- chullBuffer(crop_dir, occFile, outFolder, 500000)
        zipWrite(NAGrid, outFolder, "narea.asc.gz")
      }
    }else{
      #Reading polygon shapefile and mask  
      cat("Reading and converting \n")
      pol <- readShapePoly(shpName)
      rs <- raster(paste(crop_dir, "/masks/mask.asc", sep=""))
      
      #pa <- polygonsToRaster(pol, rs)
      pa <- rasterize(pol,rs)
      
      pa[which(!is.na(pa[]))] <- 1
      pa[which(is.na(pa[]) & rs[] == 1)] <- 0
      
      cat("Writing output \n")
      pa <- zipWrite(pa, outFolder, "narea.asc.gz")
      return(pa)
    }
    
  } else {
    cat("Already processed \n")
    #pa <- zipRead(outFolder, "narea.asc.gz")
  }
}

inDir <- paste(crop_dir,"/biomod_modeling",sep="")
#spList <- list.files(paste(inDir, "/native-areas/polyshps", sep=""))
spList <- read.csv(paste(crop_dir,"/sample_counts/sample_count_table.csv",sep=""))
spList <- spList$TAXON

for (spp in spList) {
  cat("Processing", spp, "\n")
  ot <- createNARaster(spp, inDir)
}
