
# Select background points based in SSB correction (Hijmans, 2012)
# CIAT, 2014

require(rgdal)
require(dismo)
require(raster)
require(BIOMOD)

source(paste(src.dir,"/000.zipRead.R",sep=""))

### Function
selectBack <- function(taxon, occFile, outBackDir) {
  
  if(file.exists(occFile)) {
    
    cat("\n Selecting ...\n")
    
    # Name of taxon
    taxon <- strsplit(taxon, split=".csv")
    taxon <- unlist(taxon)
    
    # Occurrence data
    spData <- read.csv(occFile)
    
    # Load native area
    na_dir <- paste(ascis,"/",taxon,sep="")
    verf <- paste(na_dir, "/narea.asc.gz", sep="")
    
    if(file.exists(verf)){
      NArea <- zipRead(na_dir,"narea.asc.gz")
      rm(na_dir)
      NArea[which(NArea[]<1)] <- NA
    }; rm(verf)
    
    # Cleaning occurrence data
    spData <- spData[!duplicated(spData[]),]
    aux    <- extract(NArea,spData[,c("lon","lat")])
    spData$aux <- aux; rm(aux)
    spData <- spData[!is.na(spData$aux),]
    spData$aux <- NULL
    
    # Option 1
    # require(SpatialTools)
    # Distances between all points
    # distances <- dist1(spData[,c("lon","lat")])
    # distances <- distances[upper.tri(distances)]
    # index     <- sum(distances)
    
    # Option 2
    require(fields)
    # Distances to the mean on coordinates
    xy <- colMeans(spData[,c("lon","lat")])
    xy <- as.data.frame(t(xy))
    
    distances <- rdist(xy, spData[,c("lon","lat")])
    index <- sum(distances); rm(xy); rm(distances)
    
    if(index > 2000){
      
      # Extract coordinates from native area
      paData <- xyFromCell(NArea,which(NArea[]==1))
      colnames(paData) <- c("lon","lat")
      
      paData <- as.data.frame(paData)
      paData$Taxon <- taxon
      paData <- paData[,c("Taxon","lon","lat")]
      
      # Points of Presences and background area
      all_data <- rbind(spData[,c("Taxon","lon","lat")],
                        paData[,c("Taxon","lon","lat")])
      
      # Pseudo-absences selection
      # 5 sets of 10000 points within native area
      for(i in 1:5){
        
        # Random strategy
        pseudo.abs(coor           = all_data[,c("lon","lat")],
                   status         = c(rep(1,dim(spData)[1]),rep(0,dim(paData)[1])),
                   strategy       = "random",
                   create.dataset = TRUE,
                   nb.points      = 10000,
                   species.name   = paste(taxon,"_",i,sep=""),
                   add.pres       = TRUE)
        
        lab <- paste(taxon,"_PA_",i," <- Dataset.",taxon,"_",i,".random.partial",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("rm(Dataset.",taxon,"_",i,".random.partial)",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("rownames(",taxon,"_PA_",i,") <- 1:nrow(",taxon,"_PA_",i,")",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("colnames(",taxon,"_PA_",i,") <- c('lon','lat','id_pa')",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste(taxon,"_PA_",i,"$Taxon <- ","taxon",sep="")
        eval(parse(text=lab)); rm(lab)
        
      }; rm(i)
      
      # SSB correction function
      ssb_correction <- function(data, bios, ...){
        
        require(dismo)
        data = data
        
        # Presences and pseudo-absences (ID's of rows)
        presences <- which(data$id_pa==1)
        absences  <- which(data$id_pa==0)
        
        # Training and testing data for presences
        # set.seed(.)
        pres.train <- sample(presences,size=round(0.7*length(presences),0))
        pres.test  <- setdiff(presences,pres.train)
        
        # Training and testing data for absences
        # set.seed(.)
        abs.train <- sample(absences,size=round(0.7*length(absences),0))
        abs.test  <- setdiff(absences,abs.train)
        
        # Create ID_train variable
        data$id_train <- 0
        data$id_train[c(pres.train,abs.train)] <- 1
        
        # Presences training
        p_train <- data.frame(lon = data$lon[data$id_pa==1 & data$id_train==1],
                              lat = data$lat[data$id_pa==1 & data$id_train==1])
        
        # Presences testing
        p_test  <- data.frame(lon = data$lon[data$id_pa==1 & data$id_train==0],
                              lat = data$lat[data$id_pa==1 & data$id_train==0])
        
        # Pseudo-absences training
        a_train <- data.frame(lon = data$lon[data$id_pa==0 & data$id_train==1],
                              lat = data$lat[data$id_pa==0 & data$id_train==1])
        
        # Pseudo-absences testing
        a_test  <- data.frame(lon = data$lon[data$id_pa==0 & data$id_train==0],
                              lat = data$lat[data$id_pa==0 & data$id_train==0])
        
        # Measure to SSB
        SSB <- ssb(p_test, a_test, p_train)
        ssb_pa_b <- SSB[,1]/SSB[,2]; rm(SSB)
        
        cat("The Spatial Sorting Bias (SSB) before correction is:",ssb_pa_b,"\n")
        
        # Count of number of presences sites within native area
        n_pres <- sum(data$id_pa,na.rm=T)
        
        # if(n_pres >= 10){}
        
        # Correction of SSB
        absences <- data[absences,]
        rownames(absences) <- 1:nrow(absences)
        pa_select <- pwdSample(fixed     = p_test[,c("lon","lat")],
                               sample    = absences[,c("lon","lat")],
                               reference = p_train[,c("lon","lat")],
                               # 5 veces el número de presencias
                               n         = ((5*n_pres)/(0.30*n_pres)))
        rm(n_pres)
        
        pa_select <- as.vector(na.exclude(as.vector(pa_select)))
        pa_select <- absences[pa_select,]
        rownames(pa_select) <- 1:nrow(pa_select)
        
        ### Sesgo de selección espacial corregido
        SSB <- ssb(p_test, pa_select[,c("lon","lat")], p_train)
        ssb_pa_a <- SSB[,1]/SSB[,2]
        
        pa_select$id_pa    <- NULL
        pa_select$id_train <- NULL
        
        pa_biovar <- extract(bios, pa_select[,c("lon","lat")])
        pa_select <- cbind(pa_select,pa_biovar)
        pa_select <- as.data.frame(pa_select)
        pa_select$Taxon <- taxon
        pa_select <- pa_select[,c("Taxon","lon","lat",paste("bio_",1:19,sep=""))]
        
        cat("The Spatial Sorting Bias after correction is:",ssb_pa_a,"\n \n")
        
        return(pa_select)
        
      }
      
      # Create list of data frames
      dfList <- list(0)
      for(i in 1:5){
        
        lab <- paste(taxon,"_PA_",i,sep="")
        dfList[[i]] <- eval(parse(text=lab))
        rm(lab)
        
        lab <- paste("rm(",taxon,"_PA_",i,")",sep="")
        eval(parse(text=lab)); rm(lab)
        
      }; rm(i)
      
      # Climatic information
      bio_var <- list.files(env_dir,pattern="\\.asc",full.names=T)
      bio_var <- lapply(X=bio_var,FUN=raster)
      bio_var <- stack(bio_var)
      
      # Apply correction
      dfList <- lapply(dfList, function(x) ssb_correction(data=x,bios=bio_var))
      names(dfList) <- paste(taxon,"_PA_",1:5,sep="")
      
      outDir <- paste(outBackDir,"/",taxon,sep="")
      if(!file.exists(outDir)){dir.create(outDir)}
      
      # Write data
      for(i in 1:5){
        
        lab <- paste("write.csv(dfList[[",i,"]],paste(outDir,'/',names(dfList)[",i,"],'.csv',sep=''),row.names=FALSE)",sep="")
        eval(parse(text=lab)); rm(lab)
        
      }; rm(i)
      
    } else {
      
      # Extract coordinates from native area
      paData <- xyFromCell(NArea,which(NArea[]==1))
      colnames(paData) <- c("lon","lat")
      
      # Climatic information
      bio_var <- list.files(env_dir,pattern="\\.asc",full.names=T)
      bio_var <- lapply(X=bio_var,FUN=raster)
      bio_var <- stack(bio_var)
      
      paData <- as.data.frame(paData)
      paData$Taxon <- taxon
      paData <- paData[,c("Taxon","lon","lat")]
      
      # Points of Presences and background area
      all_data <- rbind(spData[,c("Taxon","lon","lat")],
                        paData[,c("Taxon","lon","lat")])
      
      count <- nrow(spData)
      
      # Pseudo-absences selection
      # 5 veces el número de presencias del taxón
      for(i in 1:5){
        
        # Random strategy
        pseudo.abs(coor           = all_data[,c("lon","lat")],
                   status         = c(rep(1,dim(spData)[1]),rep(0,dim(paData)[1])),
                   strategy       = "random",
                   create.dataset = TRUE,
                   nb.points      = 5*count,
                   species.name   = paste(taxon,"_",i,sep=""),
                   add.pres       = TRUE)
        
        lab <- paste(taxon,"_PA_",i," <- Dataset.",taxon,"_",i,".random.partial",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("rm(Dataset.",taxon,"_",i,".random.partial)",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("rownames(",taxon,"_PA_",i,") <- 1:nrow(",taxon,"_PA_",i,")",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("colnames(",taxon,"_PA_",i,") <- c('lon','lat','id_pa')",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste(taxon,"_PA_",i,"$Taxon <- ","taxon",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("bios <- extract(bio_var,",taxon,"_PA_",i,"[,c('lon','lat')])",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste(taxon,"_PA_",i," <- cbind(",taxon,"_PA_",i,",bios)",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste(taxon,"_PA_",i," <- as.data.frame(",taxon,"_PA_",i,")",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("rm(bios)",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste("rownames(",taxon,"_PA_",i,") <- 1:nrow(",taxon,"_PA_",i,")",sep="")
        eval(parse(text=lab)); rm(lab)
        
        lab <- paste(taxon,"_PA_",i," <- ",taxon,"_PA_",i,"[,c('Taxon','lon','lat',paste('bio_',1:19,sep=''))]",sep="")
        eval(parse(text=lab)); rm(lab)
        
      }; rm(i)
      
      # Create list of data frames
      dfList <- list(0)
      for(i in 1:5){
        
        eq <- paste(taxon,"_PA_",i,sep="")
        dfList[[i]] <- eval(parse(text=eq))
        rm(eq)
        
        eq <- paste("rm(",taxon,"_PA_",i,")",sep="")
        eval(parse(text=eq)); rm(eq)
        
      }; rm(i)
      
      names(dfList) <- paste(taxon,"_PA_",1:5,sep="")
      
      outDir <- paste(outBackDir,"/",taxon,sep="")
      if(!file.exists(outDir)){dir.create(outDir)}
      
      # Write data
      for(i in 1:5){
        
        lab <- paste("write.csv(dfList[[",i,"]],paste(outDir,'/',names(dfList)[",i,"],'.csv',sep=''),row.names=FALSE)",sep="")
        eval(parse(text=lab)); rm(lab)
        
      }; rm(i)
      
    }
    
  }
  
  return(cat("Done! \n \n"))
  
}
