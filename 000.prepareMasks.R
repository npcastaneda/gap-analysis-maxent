require(rgdal)
require(raster)
require(maptools)
gpclibPermit()

source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.zipRead.R",sep=""))

createMasks <- function(inDir,env_dir){
  cat("Preparing the mask! \n")
  
  inPolyDir <- paste(inDir,"/polyshps",sep="")
  
  if (!file.exists(paste(inDir, "/mask.asc", sep=""))) {
    
    cat("Layer is not processed, thus processing \n")
    shpName <- paste(inPolyDir, "/mask.shp", sep="")
    
    #Reading polygon shapefile and mask
    cat("Reading and converting \n")
    pol <- readShapePoly(shpName)
    ls <- list.files(env_dir, pattern="bio")    
    rs <- raster(paste(env_dir, "/", ls[1],sep=""))
    
    pa <- rasterize(pol,rs)
    pa <- trim(pa)
    
#     pa[which(!is.na(pa[]))] <- 1
    pa <- reclassify(pa,c(minValue(pa), maxValue(pa), 1))
    #pa[which(is.na(pa[]) & rs[] == 1)] <- 0
    
    # Prepare cellArea.asc
    rs_a <- area(pa)
    rs_a <- mask(rs_a, pa)

    cat("Writing outputs \n")
    writeRaster(rs_a,paste(inDir,"/cellArea.asc",sep=""), overwrite=T)
    writeRaster(pa,paste(inDir,"/mask.asc",sep=""), overwrite=T)
    
  } else {
    cat("Already processed \n")
  }
  
}
