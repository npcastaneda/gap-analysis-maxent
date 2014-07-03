require(rgdal)
require(raster)

source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.bufferPoints.R",sep=""))

#When presence surface not available, get all the populations and use the 50km buffer

speciesRichness_alt <- function(bdir) {
  idir <- paste(bdir, "/maxent_modeling", sep="")
  ndir <- paste(bdir, "/biomod_modeling/native-areas/asciigrids",sep="")
  mask <- raster(paste(bdir, "/masks/mask.asc", sep="")) # New line
  ddir <- paste(bdir, "/samples_calculations", sep="")
  
  outFolder <- paste(bdir, "/species_richness", sep="")
  if (!file.exists(outFolder)) {dir.create(outFolder)}
  
  spList <- read.csv(paste(bdir, "/summary-files/taxaForRichness.csv", sep=""))
  allOcc <- read.csv(paste(bdir, "/occurrences/",crop,".csv", sep=""))
  
  # Listing taxa without sampling buffer files
  spList_buffer=spList[spList$IS_VALID==0,]
  for (spp in spList_buffer$TAXON) {
    if(file.exists(paste(ddir, "/", spp, "/samples-buffer-na.asc.gz", sep=""))){
      spList_buffer=spList_buffer[spList_buffer$TAXON!=spp,]
    }else{
      spList_buffer=spList_buffer
    }    
  }
  
  # Creating samples buffer files for taxa requiring it
  if(dim(spList_buffer)[1]==0){
    print("Calculating samples buffer files is not necessary \n")
  }else{
    for (spp in spList_buffer$TAXON) {
      
      cat("Creating samples buffer for", spp, "\n")
      tallOcc <- allOcc[which(allOcc$Taxon == paste(spp)),]
      if (nrow(tallOcc) != 0) {
        exOcc <- T
        spOutFolder <- paste(ddir, "/", spp, sep="")
        if (!file.exists(spOutFolder)) {dir.create(spOutFolder)}
        tallOcc <- as.data.frame(cbind(as.character(tallOcc$Taxon), tallOcc$lon, tallOcc$lat))
        names(tallOcc) <- c("taxon", "lon", "lat")
        
        write.csv(tallOcc, paste(spOutFolder, "/samples.csv", sep=""), quote=F, row.names=F)
        rm(tallOcc)
        pagrid <- createBuffers(paste(spOutFolder, "/samples.csv", sep=""), spOutFolder, "samples-buffer.asc", 50000, paste(bdir, "/masks/mask.asc", sep=""))
        pagrid <- zipRead(spOutFolder,"samples-buffer.asc.gz") * zipRead(paste(ndir,"/",spp,sep=""),"narea.asc.gz")
        pagrid <- zipWrite(pagrid,spOutFolder,"samples-buffer-na.asc.gz")
      } else {
        cat("No occurrence points for ", spp,"\n")
      }
    }
  }
  
  spList_buffer=spList[spList$IS_VALID==0,]
  nrow(spList_buffer)
  for (spp in spList_buffer$TAXON) {
    if(!file.exists(paste(ddir, "/", spp, "/samples-buffer-na.asc.gz", sep=""))){
      spList_buffer=spList_buffer[spList_buffer$TAXON!=spp,]}
  }
  
  results_0=list()
  if(length(spList_buffer$TAXON)==0){
    cat("All species have valid niche models \n")
  }else{
    for(i in 1:length(spList_buffer$TAXON)){  
      results_0[[i]]=zipRead(paste(ddir, "/", spList_buffer$TAXON[i], sep=""),"samples-buffer-na.asc.gz")}
  }
  
  spList_buffer=spList[spList$IS_VALID==1,]
  results_1=list()
  if(length(spList_buffer$TAXON)==0){
    cat("None species have valid niche models \n")
  }else{
    for(spp in spList_buffer$TAXON){
      pos=which(spList_buffer$TAXON==spp)
      sppFolder <- paste(idir, "/models/", spp, sep="")
      projFolder <- paste(sppFolder, "/projections", sep="")
      pagrid <- paste(spp, "_worldclim2_5_EMN_PA.asc.gz", sep="")
      results_1[[pos]] <- zipRead(projFolder, pagrid)}
  }
  
  cat("Writing richness raster \n")
  if(length(results_0)==0){
    cat("No samples buffers to map \n")
  }else if(length(results_0)==1){
    results_sum_0=results_0[[1]]
  }else{
    results_sum_0=results_0[[1]]
    for(i in 2:length(results_0)){
      results_sum_0=extend(results_sum_0, mask) # New line
      results_0[[i]] = extend(results_0[[i]], mask) # New line
      results_sum_0=sum(results_sum_0,results_0[[i]], na.rm=T)
    }
  }
  
  if(length(results_1)==0){
    cat("No valid distribution models to map \n")
  }else if(length(results_1)==1){
    results_sum_1=results_1[[1]]
  }else{
    results_sum_1=results_1[[1]]
    for(i in 2:length(results_1)){
      results_sum_1=extend(results_sum_1, mask) # New line
      results_1[[i]] = extend(results_1[[i]], mask) # New line
      results_sum_1=sum(results_sum_1,results_1[[i]], na.rm=T)
    }
  }

    if(sum(ls()=="results_sum_0")){
      cat("results_sum_0 available \n")
      results_sum=extend(results_sum_0, mask)
      if(sum(ls()=="results_sum_1")){
        cat("results_sum_0 and results_sum_1 available \n")
        results_sum <- sum(results_sum, results_sum_1, na.rm=T)
      }
    }else{
      cat("only results_sum _1 available \n")
      results_sum <- results_sum_1
    }
    
    cat("Calculating mean sd raster \n")
    results_mean=results_sum/(length(results_0)+length(results_1))
    
    if(length(results_0)!=0){
      results_sd_0=(results_mean-results_0[[1]])^2
      
      if(length(results_0)==1){
        results_sd_0=(results_mean-results_0[[1]])^2
        
      }else{
        for(i in 2:length(results_0)){
          sum_sqrt=(results_mean-results_0[[i]])^2
          results_sd_0=sum(results_sd_0,sum_sqrt, na.rm=T)}
      }
    }
    
    if(length(results_1)!=0){
      results_sd_1=(results_mean-results_1[[1]])^2
      
      if(length(results_1)==1){
        results_sd_1=(results_mean-results_1[[1]])^2
        
      }else{
        for(i in 2:length(results_1)){
          sum_sqrt=(results_mean-results_1[[i]])^2
          results_sd_1=sum(results_sd_1,sum_sqrt, na.rm=T)}
      }
    }
  
    if(sum(ls()=="results_sd_0")){
      cat("results_sd_0 available \n")
      results_sd=extend(results_sd_0, mask)
      if(sum(ls()=="results_sd_1")){
        cat("results_sd_0 and results_sd_1 available \n")
        results_sd <- sum(results_sd, results_sd_1, na.rm=T)
      }
    }else{
      cat("only results_sd _1 available \n")
      results_sd <- results_sd_1
    }
 
    results_mean_sd=results_sd/(length(results_0)+length(results_1))
    
    zipWrite(results_sum, outFolder, "species-richness.asc.gz")
    zipWrite(results_mean_sd, outFolder, "species-richness-sdmean.asc.gz")
    
    cat("Done! \n")
}