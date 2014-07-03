#Take the list of species to model and create a new file with occurrences of selected species and climate data
# N. Castaneda 2012

require(raster)

createOccFilesBio <- function(occ, taxfield, outDir, env.dir) {
  
	cat("\n")
	
	if (!file.exists(outDir)) {
		dir.create(outDir)
	}
	
	cat("Loading occurrence file \n")
	
	occ <- read.csv(occ)
	
	spL <- unique(occ[,taxfield])
	
	spcounter <- 1
  
	# Stack climate variables
	for(i in 1:19){
	  if (i==1){
	    rs <- raster(paste(env.dir,"/bio_",i,".asc",sep=""))
	    bio.stk <- stack(rs)
	  }else{
	    rs <- raster(paste(env.dir,"/bio_",i,".asc",sep=""))
	    bio.stk <- addLayer(bio.stk,rs)
	  }
	}
	
	cat("Now printing \n")
  
  for (sp in spL) {
		
		nspp <- length(spL)
		cat("      ...", paste(round(spcounter/nspp*100,2),"% Completed", sep=""), "\n")
		
		spData <- occ[which(occ[,taxfield] == sp),]
    spData <- subset(spData,select=c("lon","lat"))
    
		# Points to Raster
		mask <- raster(paste(env.dir,"/bio_1.asc",sep=""))
		spData <- rasterize(spData, mask)
    spData[spData[]>=1]=1
    
    # Add occurrences layer to bio.stk
		bio.stk <- addLayer(bio.stk, spData)
    
    # Extract coordinates
		xy <- xyFromCell(mask,which(!is.na(mask[])))
    
    # Arrange matrix
		v <- getValues(bio.stk)
		i <- which(! is.na(v[,1]))
		spData <- v[i,]
		spData[is.na(spData)] <- 0 # Affects the species column
		
		spData <- cbind(spData,xy)
    #rm(xy)
		
		spData <- as.data.frame(spData)
		
		Sp.Env <- spData
		#head(Sp.Env)
		
    # Save files
    
    csvName <- paste(outDir, "/", sp, ".csv", sep="")
		write.csv(spData, csvName, row.names=F, quote=F)
		rm(spData)
		
		# Clean stack
    bio.stk <- dropLayer(bio.stk,c(20))
    
    spcounter <- spcounter+1
	}
	return("Done")
}
