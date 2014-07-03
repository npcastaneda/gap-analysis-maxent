#Creating maps and figures from the gap analysis
#March 2013
# Improvements by M. Syrfert (NHM) and H. Achicanoy(CIAT)

library(raster); library(rasterVis); library(maptools); data(wrld_simpl)

source(paste(src.dir,"/000.zipRead.R",sep=""))

### Species richness map
sp_rich_dir <- paste(crop_dir,"/species_richness",sep="")
sp_rich <- zipRead(sp_rich_dir, "species-richness.asc.gz")

# cols = colorRampPalette(c("gray87", "dark green","yellow","orange","red"))(255) # LOVE IT!!!
cols = colorRampPalette(c("white", "dark green","yellow","orange","red"))(255) # LOVE IT!!!
z <- extent(sp_rich)
aspect <- (z@ymax-z@ymin)*1.4/(z@xmax-z@xmin)

tiff("./figures/spp_richness.tif",
     res=300,pointsize=5,width=1500,height=1500*aspect,units="px",compression="lzw")
plot(sp_rich, col=cols, useRaster=T, 
     main="Species richness",
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5, border="azure4")
grid()
dev.off()

### Gap richness map
priorityLevels <- c("HPS","MPS","LPS","NFCR")
for(p in priorityLevels){
  gap_dir <- paste(crop_dir, "/gap_richness/", p, sep="")
  if(!file.exists(paste(gap_dir,"/gap-richness.asc.gz",sep=""))){
    cat("No raster file available for category:",p,"\n")
  }else{
    gap_ras <- zipRead(gap_dir,"gap-richness.asc.gz")
    
    cols = colorRampPalette(c("gray87", "dark green","yellow","orange","red"))(255) # LOVE IT!!!
    z <- extent(gap_ras)
    aspect <- (z@ymax-z@ymin)*1.4/(z@xmax-z@xmin)
    
    tiff(paste("./figures/gap_richness_",p,".tif",sep=""),
         res=300,pointsize=5,width=1500,height=1500*aspect,units="px",compression="lzw")
    plot(gap_ras, col=cols, useRaster=T, 
         main=paste(p,"gap richness",sep=" "),
         horizontal=T,
         legend.width=1,
         legend.shrink=0.99)
    plot(wrld_simpl,add=T,lwd=0.5, border="azure4")
    grid()
    dev.off() 
  }
}

### Summary graphs

bdir <- crop_dir
dir <- paste(bdir,"/gap_richness", sep="")

summ <- c("gap-richness-dpmax","gap-richness")

for(s in summ){
  rasterList <- list.dirs(dir)[-1]
  path = rasterList
  fname = paste(s,".asc.gz",sep="")
  rs <- lapply(rasterList,zipRead,fname)
  
  rs <- stack(rs)
  
  names(rs) <- c("a", "c", "b", "d")
  rs <- rs[[order(names(rs))]]
  
  tiff(paste("./figures/all-",s,".tif",sep=""),
       res=300,width=1500,height=1500,pointsize=5,units="px",compression="lzw")            
  # Color blind safe options: PuOrTheme
  mypalette <- rasterTheme(region=brewer.pal(n=9, "YlOrRd")) #Cool!
  levelplot(rs, par.settings=mypalette, names.attr=c("HPS", "MPS", "LPS", "NFCR"))
  dev.off()  
}

### Plot the CA50 vs Potential coverage

# prior <- read.csv(paste(crop_dir,"/maxent_modeling/summary-files/areas.csv",sep=""))
# 
# 
# 
# prior <- read.csv("D:/_code/R/gap-analysis-maxent-cip/areas.csv")
# table_base <- read.csv("D:/_code/R/gap-analysis-maxent-cip/areas.csv")
# #fit <- lm(prior$CA50_G~prior$PD_COV)
# # fit <- lm(prior$GBSize~prior$DRSize)
# 
# head(table_base)
# table_base <- data.frame(Taxon=table_base$TAXON)
# table_base$HS <- NA; table_base$HS_RP <- NA
# table_base$GS <- NA; table_base$GS_RP <- NA
# table_base$TOTAL <- NA; table_base$TOTAL_RP <- NA
# table_base$ATAUC <- NA; table_base$STAUC <- NA; table_base$ASD15 <- NA; table_base$IS_VALID <- NA
# table_base$SRS <- NA; table_base$GRS <- NA; table_base$ERS <- NA
# table_base$ERTS <- NA; table_base$FPS <- NA; table_base$FPCAT <- NA
# 
# 
# 
# rsize <- prior
# 
# head(prior)
# for (spp in prior$taxon) {
#   cat("Preparing species",paste(spp),"\n")
# #   spp <- prior$taxon[1]
# #   prior$GBSize[which(rsize$taxon==paste(spp))]
#   if(prior$DRSize=="NA"){
#     size <- rsize$CHSize[which(rsize$taxon==paste(spp))]
#   } else {
#     size <- rsize$DRSize[which(rsize$taxon==paste(spp))]
#   } 
#   prior$fill[which(prior$taxon==paste(spp))] <- size
# }
# 
# head(prior)
# 
# 
# if (!is.na(grs)) {
#   if (grs>10) {grs <- 10}
# }
# 
# 
# 
# 
# prior$fill <- NA
# if(){
#   
# }
# 
#   
# fit <- lm(prior$GBSize/1000~fill/1000)
# #lims <- c(min(prior$PD_COV,prior$CA50_G),max(prior$CA50_G,prior$PD_COV))/1000
# # lims <- c(min(prior$DRSize,prior$GBSize),max(prior$GBSize,prior$DRSize))/1000
# 
# lims <- c(min(fill/1000,prior$GBSize/1000),max(fill/1000,prior$GBSize/1000))
# 
# # this section checks if DRSize is a NA, if yes, then use MCP of species (CHSize) else keep the value
# spp= prior$taxon  
# fill = numeric(length(spp))
# 
# for (i in 1:length(spp)){
#   if (is.na(prior$DRSize[i])){
#     prior$DRSize[i] <- prior$CHSize[i]
#   } else { 
#     prior$DRSize[i]       
#   } 
#   fill[i]= prior$DRSize[i]  
# }
# 
# 
# #do the plot
# tiff(paste(crop_dir,"/figures/geographic_coverage.tif",sep=""),
#      res=300,pointsize=12,width=1500,height=1000,units="px",compression="lzw")
# 
# par(mar=c(5,5,1,1),cex=0.8)
# plot(prior$DRSize/1000,prior$GBSize/1000,pch=20,cex=0.75,col="red",xlim=lims,ylim=c(0,1000),
#      xlab="Potential geographic coverage (sq-km * 1000)",
#      ylab="Genebank accessions CA50 (sq-km * 1000")
# 
# abline(0,1,lwd=0.75,lty=2)
# lines(prior$DRSize/1000,fit$fitted.values/1000)
# grid(lwd=0.75)
# 
# dev.off()
# 
# 
# 
# 
# 
# #do the plot
# tiff(paste(crop_dir,"/figures/geographic_coveragetest.tif",sep=""),
#      res=300,pointsize=12,width=1500,height=1000,units="px",compression="lzw")
# 
# par(mar=c(5,5,1,1),cex=0.8)
# plot(fill/1000,prior$GBSize/1000,pch=20,cex=0.75,col="red",xlim=lims,ylim=c(0,max(prior$GBSize/1000)+10),
#      xlab="Potential geographic coverage (sq-km * 1000)",
#      ylab="Genebank accessions CA50 (sq-km * 1000)")
# 
# abline(0,1,lwd=0.75,lty=2)
# lines(fill/1000,fit$fitted.values)
# grid(lwd=0.75)
# 
# dev.off()