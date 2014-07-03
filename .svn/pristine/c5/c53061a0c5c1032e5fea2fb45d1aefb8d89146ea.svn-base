###############################################
# BIOMOD
# N. Castaneda - 2012
###############################################


###############################################
# SETTINGS
source(paste(src.dir,"/000.zipRead.R",sep=""))
source(paste(src.dir,"/000.zipWrite.R",sep=""))
source(paste(src.dir,"/000.createChullBuffer.R",sep=""))

require(BIOMOD)
###############################################

###############################################
# TEST AREA
x <- ModelingProcess("Avena_abyssinica", "NT", crop_dir)
###############################################

ModelingProcess <- function(sp, OSys, inputDir){

  idir <- paste(inputDir,"/biomod_modeling/models",sep="")
  if (!file.exists(idir)) {dir.create(idir)}
  setwd(idir)
  
  mskDir <- paste(inputDir, "/masks", sep="")
  NADir <- paste(inputDir, "/biomod_modeling/native-areas/asciigrids", sep="")
    
    # Load species data
    Sp.Env <- read.csv(paste(inputDir,"/biomod_modeling/occurrence_files/",sp,".csv",sep=""))
    
    # Data organization
    LatLong <- Sp.Env[,21:22] # coordinates
    Expl.Var <- Sp.Env[,1:19] # bioclimatic variables
    Resp.Var <- Sp.Env[,20] # species occurrences
    
    # Prepare data for modelling
      Initial.State(Response = Resp.Var, Explanatory=Expl.Var, sp.name=sp,
                  IndependentResponse =NULL, IndependentExplanatory=NULL)
    
    # Running the models + Evaluation + Calibration
    Models(GLM = T, TypeGLM = "poly", Test = "AIC",
           GAM = F, Spline = 3, 
           CTA = F, CV.tree = 50,
           ANN = T, CV.ann = 3,
           SRE = F, quant=0.025,
           GBM = F, No.trees = 2000,
           RF = F,
           FDA = F,
           MARS = F,
           # Calibration
           # OJO:  - CHANGE NbRunEval accordingly (25?) / set nb.absences according to amount of records!!!!
           #NbRunEval = 5, DataSplit = 70, Yweights=NULL,
           NbRunEval = 2, DataSplit = 70, Yweights=NULL,
           NbRepPA=1, strategy="sre", coor=NULL, distance=2, nb.absences=10000,
           #Evaluation
           Roc = T, Optimized.Threshold.Roc = T, Kappa = F, TSS=T,
           KeepPredIndependent = T, VarImport=5)
    
    # Creating presence/absence binaries
    CurrentPred(GLM=F, GBM=F, GAM=F, CTA=F, ANN=T, SRE=F, FDA=F, MARS=T, RF=T, 
                BinRoc=T, BinKappa=F,BinTSS=T)
    
    # Projecting data
    Proj.name <- sp
    Proj <- Expl.Var
    Projection(Proj = Proj, Proj.name, GLM = T, GBM = F, GAM = T, CTA = F, ANN = T, 
               SRE = F, quant=0.025, FDA = F, MARS = F, RF = F, BinRoc = T, BinKappa = F, 
               BinTSS = T, FiltRoc = F, FiltKappa = F, FiltTSS = F, repetition.models=T, 
               compress="xz")
    
    # Ensamble Forecasting
    Ensemble.Forecasting(ANN = T, CTA = F, GAM = F, GBM = F, GLM = T, MARS = F, FDA = F,
                         RF = F, SRE = F, Proj.name = Proj.name, weight.method="Roc", 
                         decay = 1.6, PCA.median = F, binary = T, bin.method = "Roc", 
                         Test = T, repetition.models=T, final.model.out=F, 
                         qual.th=0, compress="xz")
    
    # Preparing raster files
    grid.dir <- paste(idir,"/proj.",sp,"/grdfiles",sep="")
    if (!file.exists(grid.dir)) {dir.create(grid.dir)}
    
    mask <- raster(paste(inputDir,"/masks/mask.asc",sep=""))
    
    # Binaries
    load(paste("proj.",sp,"/Total_consensus_",sp,"_Bin",sep=""))
    sp.array <- (get(paste("Total_consensus_",sp,"_Bin",sep="")))
    sp.array <- sp.array[,1,2] # Extracting weighted mean
    sp.array <- cbind(LatLong,sp.array)
    
    writeRaster(mask, paste(grid.dir,"/",sp,"_PA",sep=""), overwrite=T)
    rs <- raster(paste(grid.dir,"/",sp,"_PA",sep=""))
    
    coords <- sp.array[,1:2]
    cells <- cellFromXY(rs,coords)
    vals <- sp.array[,3]
    update(rs,vals,cells)
    writeRaster(rs, paste(grid.dir,"/",sp,"_PA.asc",sep=""), overwrite=T)
      
  #Now cut to native areas
  #Verify if the native area exists, else create one using the buffered convex hull
  
  NAGridName <- paste(NADir, "/", sp, "/narea.asc.gz", sep="")
  if (!file.exists(NAGridName)) {
    cat("The native area does not exist, generating one \n")
    occFile <- paste(inputDir,"/occurrence_files/",sp,".csv",sep="")
    NAGrid <- chullBuffer(inputDir, occFile, paste(NADir, "/", sp, sep=""), 500000)
  } else {
    cat("The native area exists, using it \n")
    NAGrid <- zipRead(paste(NADir, "/", sp, sep=""), "narea.asc.gz")
  }
  
  consensusBin <- rs
  consensusBinNA <- consensusBin * NAGrid

  # Writing Native Areas cut grids
  writeRaster(consensusBinNA, paste(grid.dir,"/",sp,"_PA_NA.asc",sep=""), overwrite=T)
  
  # Erase diva-grid format files
  a <- paste(grid.dir,"/",sp,"_PA.grd",sep=""); unlink(a)
  a <- paste(grid.dir,"/",sp,"_PA.gri",sep=""); unlink(a)
  rm(a)
  
  # Writing evaluation data
  metDir <- paste(inputDir,"/biomod_modeling/models/proj.",sp,"/metrics",sep="")
  if (!file.exists(metDir)) {dir.create(metDir)}
  
  load(paste("proj.",sp,"/consensus_",sp,"_results",sep=""))
  met.sp <- (get(paste("consensus_",sp,"_results",sep="")))
  met.sp <- lapply(met.sp[1],"[",4) # Selects all values from test.results
  met.sp <- as.data.frame(met.sp)
  
  # Calculate AUC standard deviation
  std.auc <- apply(met.sp,1,sd) 
  std.auc <- as.matrix(std.auc)
  colnames(std.auc) <- sp
  r.names <- rownames(std.auc)
  r.names <- paste("sd.",r.names,sep="")
  rownames(std.auc) <- r.names
  
  # Calculate average AUC among replicates
  met.sp <- rowMeans (met.sp, na.rm = TRUE, dims = 1) 
  met.sp <- as.data.frame(met.sp)
  colnames(met.sp) <- sp
  
  # Join STD and AUC and then transpose
  met.sp <- rbind(met.sp, std.auc)
  met.sp <- t(met.sp)
  
  write.csv(met.sp, paste(metDir,"/metrics.csv",sep=""))

    #NAGridName <- paste(NADir, "/narea.asc.gz", sep="")
    #if (!file.exists(NAGridName)) {
    #  cat("The native area does not exist, generating one \n")
      #NAGrid <- chullBuffer(inputDir, occFile, paste(NADir, "/", spID, sep=""), 500000)
    #} else {
    #  cat("The native area exists, using it \n")
    #  NAGrid <- zipRead(paste(NADir, "/", sep=""), "narea.asc.gz")
    #}
    
    #distMeanPA <- rs * NAGrid
    #writeRaster(distMeanPA, paste(grid.dir,"/",sp,"_PA_NA.asc",sep=""), overwrite=T)
    
    # Probabilities surfaces ***NEEDS TO BE THRESHOLDED***
    #load(paste("proj.",sp,"/Total_consensus_",sp,sep=""))
    #sp.array <- (get(paste("Total_consensus_",sp,sep="")))
    #sp.array <- sp.array[,1,2] # Extracting weighted mean
    #sp.array <- cbind(LatLong,sp.array)
    
    # Opción D. Arango
    #writeRaster(mask, paste(grid.dir,"/",sp,sep=""), overwrite=T)
    #rs <- raster(paste(grid.dir,"/",sp,sep=""))
    
    #coords <- sp.array[,1:2]
    #cells <- cellFromXY(rs,coords)
    #vals <- sp.array[,3]
    #update(rs,vals,cells)
}  

GapProcess <- function(inputDir, OSys="LINUX", ncpu){
  
  spL <- list.files(paste(inputDir,"/biomod_modeling/occurrence_files",sep=""), pattern=".csv")

  gap_wrapper <- function(i) {
    library(raster)
    library(SDMTools)
    library(BIOMOD)
    sp <- spL[i]
    sp <- unlist(strsplit(sp, ".", fixed=T))[1]
    cat("\n")
    cat("...Species", sp, "\n")
    out <- ModelingProcess(sp, OSys, inputDir)
    }
  
  library(snowfall)
  sfInit(parallel=T,cpus=ncpu)
  
  sfExport("gap_wrapper")
  sfExport("ModelingProcess")
  #sfExport("getMetrics")
  sfExport("zipRead")
  sfExport("zipWrite")
  sfExport("chullBuffer")
  sfExport("inputDir")
  sfExport("OSys")
  sfExport("spL")
  
  #run the control function
  system.time(sfSapply(as.vector(1:length(spL)), gap_wrapper))
  
  #stop the cluster
  sfStop()
  
  return("Done!")
  
}
