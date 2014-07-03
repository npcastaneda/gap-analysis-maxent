# Prepare collecting uncertainties maps
# Castaneda - 2014

# Initial settings
require(rgdal)
require(raster)

source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.zipWrite.R",sep=""))

# The function
# bdir <- crop_dir
maxpdist <- function(bdir){
  
  pointDir <- paste(bdir,"/samples_calculations", sep="")
  
  priorityLevel <- c("HPS","MPS","LPS","NFCR")
  
  for(p in priorityLevel){
    gdir <- paste(bdir, "/gap_spp/",p,sep="")
    odir <- paste(bdir, "/gap_richness/", p, sep="")
    spList <- read.csv(paste(bdir, "/priorities/priorities.csv", sep=""))
    
    names(spList)[1]="TAXON"
    spList <- spList[which(spList$FPCAT == p),]
    
    cat("\n")
    cat("Finding the maximum geographic distance to a known population among", nrow(spList), p, "\n")
    
    if(nrow(spList) != 0){
      spList <- list.files(gdir, pattern=".asc.gz")
      spList <- strsplit(spList,".asc.gz") # MASO!
      spList <- paste(pointDir,"/",spList,sep="")
      
      path = spList
      fname="pop-dist.asc.gz"
      dplist <- lapply(spList,zipRead,fname)
      
      cat("Calculating max pdist raster \n")
      dpmax <- max(stack(dplist))
      cat("Writing \n")
      dumm <- zipWrite(dpmax, odir, "gap-richness-dpmax.asc.gz")
      
    } else {
      cat("No taxa under this prioritization category \n")
    } 
  }
}

