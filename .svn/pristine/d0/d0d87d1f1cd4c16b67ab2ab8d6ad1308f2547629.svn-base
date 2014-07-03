# /* This program creates a SWD (sampling with data) file for a sample file in order to produce for maxent input
# /* Use this program to extract from bioclim datasets and produce a very reliable model and then you could project it
# /* to any geographic area, resolution or timeslice you may desire.
# 
# /* Written by Julian Ramirez
# /* CIAT, March 2012

library(raster)

#function to extract climate data
extractClimates <- function(input_dir,sample_file,env_dir,env_prefix,
                            env_ext,idfield,lonfield,latfield,taxfield,
                            output_dir) {
  occ <- read.csv(paste(input_dir,"/",sample_file, sep=""))
  xy <- data.frame(X=occ[,lonfield],Y=occ[,latfield])
  swd <- data.frame(cbind(paste(occ[,taxfield]),occ[,lonfield],occ[,latfield]))
  names(swd) <- c(taxfield,lonfield,latfield)
  
  #for (i in 1:19) {
  for (i in 1:20) {  
    cat("Reading environmental layer",i,"\n")
    #rs <- raster(paste(env_dir,"/",env_prefix,i,env_ext,sep=""))
    rs <- raster(paste(env_dir,"/",env_prefix,i,sep=""))
    swd$NEW <- extract(rs,xy)
    names(swd)[i+3] <- paste(env_prefix,i,sep="")
  }
  
  #cleaning the dataset and writing a separate dataset for those records with any missing record
  fun <- function(x) {
    nas <- which(is.na(x))
    if (length(nas)>0) {ok <- 0} else {ok <- 1}
    return(ok)
  }
  
  #get the cleaning field in the occurrence matrix
  swdin <- swd[,4:22]
  swd$IS_OK <- apply(swdin,1,fun)
  
  #write all the stuff
  swd_wrong <- swd[which(swd$IS_OK==0),]; swd_wrong$IS_OK <- NULL
  write.csv(swd_wrong,paste(output_dir,"/occurrences_swd_wrong.csv",sep=""),quote=F,row.names=F)
  
  swd_ok <- swd[which(swd$IS_OK==1),]; swd_ok$IS_OK <- NULL
  write.csv(swd_ok,paste(output_dir,"/occurrences_swd_ok.csv",sep=""),quote=F,row.names=F)
  
  swd_all <- swd
  write.csv(swd_all,paste(output_dir,"/occurrences_swd_all.csv",sep=""),quote=F,row.names=F)
  return("done!")
}

