summarizeMetrics <- function(idir){

  #Setting the output directory
  odir <- paste(idir, "/summary-files", sep="")
  if (!file.exists(odir)) {dir.create(odir)}
  
  spList <- list.files(paste(idir,"/occurrence_files", sep=""),pattern=".csv")
  
  sppC <- 1
  for(sp in spList){
    sp <- unlist(strsplit(sp, ".", fixed=T))[1]
    
    cat(sp, "\n")
    
    metFile <- paste(idir,"/models/proj.",sp,"/metrics/metrics.csv",sep="")
    metFile <- read.csv(metFile)
    
    if (sppC == 1){
      #metFile <- paste(idir,"/models/proj.",sp,"/metrics/metrics.csv",sep="")
      metDummy <- metFile
      cat(sp, "is first", "\n")
    } else {
      metDummy <- rbind(metDummy, metFile)
      rm(metFile)
      cat(sp, "is other", "\n")
    } 
    sppC <- sppC + 1
  }

  # Writing the output
  cat("Writing the output!")
  oFile <- paste(odir, "/accuracy.csv", sep="")
  colnames(metDummy)[1] <- "spp"
  write.csv(metDummy, oFile, quote=F, row.names=F)
    
}

