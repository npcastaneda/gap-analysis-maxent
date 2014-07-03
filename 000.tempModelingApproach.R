source(paste(src.dir,"/000.getMetrics.R",sep=""))
source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.createChullBuffer.R",sep=""))

VerificationProcess <- function(spID, OSys, inputDir) {
  mxe_out <- paste(inputDir,"/maxent_modeling/models",sep="")
  if (!file.exists(mxe_out)) {dir.create(mxe_out)}
  spDir <- paste(mxe_out,"/",spID,sep="")
  
  outFolder <- paste(inputDir,"/maxent_modeling/models/",spID, sep="")
  
  distMeanPA <- paste(outFolder, "/projections/", spID, "_worldclim2_5_EMN_PA.asc.gz", sep="")
  distMeanPR <- paste(outFolder, "/projections/", spID, "_worldclim2_5_EMN_PR.asc.gz", sep="")
  distStdvPR <- paste(outFolder, "/projections/", spID, "_worldclim2_5_ESD_PR.asc.gz", sep="")
  
  if(!file.exists(distMeanPA) | !file.exists(distMeanPR) | !file.exists(distStdvPR)){
    cat("No existen los archivos de modelacion!")
    #source(paste(src.dir,"/005.modelingApproach.R",sep=""))
    #GapProcess(inputDir=crop_dir, OSys="linux", ncpu=3)
    
  }else{
    cat("Removing and recalculating the native-aread distributions \n")
    
    unlink(distMeanPA)
    unlink(distMeanPR)
    unlink(distStdvPR)
    
    NADir <- paste(inputDir, "/biomod_modeling/native-areas/asciigrids", sep="")
    
    cat("Taxon ", spID, "\n")
    
    #1. Load species data
    occFile <- paste(inputDir, "/occurrence_files/", spID, ".csv", sep="")
    
    if (file.exists(occFile)) {
      #1.1 Load the data
      inData <- read.csv(occFile)
      nOcc <- nrow(inData)
      
      if (nOcc > 0) {
        
        suffix <- "worldclim2_5"
        
        distMean <- zipRead(paste(outFolder, "/projections", sep=""), paste(spID, "_", suffix, "_EMN.asc.gz", sep=""))
        distStdv <- zipRead(paste(outFolder, "/projections", sep=""), paste(spID, "_", suffix, "_ESD.asc.gz", sep=""))
        
        thslds <- c("UpperLeftROC")
        
        threshFile <- paste(outFolder, "/metrics/thresholds.csv", sep="")
        threshData <- read.csv(threshFile)
        
        thrNames <- names(threshData)
        thePos <- which(thrNames == thslds)
        theVal <- threshData[1,thePos]
        
        cat("Thresholding... \n")
        
        distMeanPR <- distMean
        distMeanPR[which(distMeanPR[] < theVal)] <- NA
        
        distMeanPA <- distMean
        distMeanPA[which(distMeanPA[] < theVal)] <- 0
        distMeanPA[which(distMeanPA[] != 0)] <- 1
        
        distStdvPR <- distStdv * distMeanPA
        
        #Now cut to native areas
        #Verify if the native area exists, else create one using the buffered convex hull
        
        NAGridName <- paste(NADir, "/", spID, "/narea.asc.gz", sep="")
        if (!file.exists(NAGridName)) {
          cat("The native area does not exist, generating one \n")
          NAGrid <- chullBuffer(inputDir, occFile, paste(NADir, "/", spID, sep=""), 500000)
        } else {
          cat("The native area exists, using it \n")
          NAGrid <- zipRead(paste(NADir, "/", spID, sep=""), "narea.asc.gz")
        }
        
        # Crop to Native Area extent
        distMeanPA <- crop(distMeanPA, NAGrid)
        distMeanPR <- crop(distMeanPR, NAGrid)
        distStdvPR <- crop(distStdvPR, NAGrid)
        
        distMeanPA <- distMeanPA * NAGrid
        distMeanPR <- distMeanPR * NAGrid
        distStdvPR <- distStdvPR * NAGrid
        
        #Writing these rasters
        
        #distMeanPA <- writeRaster(distMeanPA, paste(outFolder, "/projections/", spID, "_", suffix, "_EMN_PA.asc", sep=""), format='ascii', overwrite=T)
        #distMeanPR <- writeRaster(distMeanPR, paste(outFolder, "/projections/", spID, "_", suffix, "_EMN_PR.asc", sep=""), format='ascii', overwrite=T)
        #distStdvPR <- writeRaster(distStdvPR, paste(outFolder, "/projections/", spID, "_", suffix, "_ESD_PR.asc", sep=""), format='ascii', overwrite=T)
        
        distMeanPA <- writeRaster(distMeanPA, paste(outFolder, "/projections/", spID, "_", suffix, "_EMN_PA.asc", sep=""), overwrite=T)
        distMeanPR <- writeRaster(distMeanPR, paste(outFolder, "/projections/", spID, "_", suffix, "_EMN_PR.asc", sep=""), overwrite=T)
        distStdvPR <- writeRaster(distStdvPR, paste(outFolder, "/projections/", spID, "_", suffix, "_ESD_PR.asc", sep=""), overwrite=T)
        
        #Compressing everything within the projection dir
        
        ftoZIP <- list.files(paste(outFolder, "/projections/", sep=""), pattern=".asc")
        cat("Compressing... \n")
        for (fz in ftoZIP) {
          fName <- paste(outFolder, "/projections/", fz, sep="")
          if (OSys == "linux") {
            system(paste("gzip", fName))
          } else {
            system(paste("7z", "a", "-tgzip", paste(fName, ".gz", sep=""), fName),wait=T)
            file.remove(fName)
          }
        }
      }else{
        cat("Species with 0 datapoints, not to be modeled \n")
      }
      
    }else{
      cat("The occurrence file does not exist! \n")
    }
  } 
}

wrapper <- function(inputDir, OSys="LINUX") {
  spList <- list.files(paste(inputDir, "/occurrence_files", sep=""),pattern=".csv")
  for(sp in spList){
    spID <- unlist(strsplit(sp,".",fixed=T))[1]
    OSys <- tolower("linux")
    out <- VerificationProcess(spID, OSys, inputDir)
  }
}

wrapper(crop_dir,OSys="LINUX")