# Gap Analysis package
# CIAT, 2009-2014
# Author: Nora Castaneda
# Contributors: Julian Ramirez-Villegas, Harold Achicanoy

# Splits the occurrence files into two files: one holding germplasm accessions and the one holding reference (mainly herbaria) sightings. It also enables the user to visualize samples density

library(raster); library(maptools); data(wrld_simpl)

occ <- read.csv(paste("./occurrences/",crop,"_all.csv",sep=""))

occ <- occ[which(!is.na(occ$lat)),]
occ <- occ[which(!is.na(occ$lon)),]

h <- occ[which(occ$H==1),]
g <- occ[which(occ$G==1),]

write.csv(h,paste("./occurrences/",crop,"_h.csv",sep=""),quote=F,row.names=F) #georeferenced herbarium records
write.csv(g,paste("./occurrences/",crop,"_g.csv",sep=""),quote=F,row.names=F) #georeferenced germplasm records
write.csv(occ,paste("./occurrences/",crop,".csv",sep=""),quote=F,row.names=F) #all georeferenced records


#==map H/G densities ==#
g_ras <- read.csv(paste("./occurrences/",crop,"_g.csv",sep=""))
h_ras <- read.csv(paste("./occurrences/",crop,"_h.csv",sep=""))


r <- raster()
res(r) <- 1

g_occ <- g_occ[,c(lon,lat)]
g_ras <- rasterize(g_occ,r,fun="count")

h_occ <- h_occ[,c(lon,lat)]
h_ras <- rasterize(h_occ,r,fun="count")

h_ras[which(h_ras[]==0)] <- NA; g_ras[which(g_ras[]==0)] <- NA

brks <- unique(quantile(c(h_ras[],g_ras[]),na.rm=T,probs=c(seq(0,1,by=0.05))))
cols <- colorRampPalette(c("dark green","yellow","orange","red"))(length(brks)-1)
brks.lab <- round(brks,0)

if (!file.exists("./figures")) {dir.create("./figures")}

z <- extent(h_ras)
aspect <- (z@ymax-z@ymin)*1.4/(z@xmax-z@xmin)

par(1,1)

#herbarium map
tiff("./figures/h_samples_count.tif",
     res=300,pointsize=5,width=1500,height=1500*aspect,units="px",compression="lzw")
par(mar=c(2.5,2.5,1,1),cex=0.8,lwd=0.8)
plot(h_ras,col=cols,zlim=c(min(brks),max(brks)), main = "Reference samples (H)",
     breaks=brks,lab.breaks=brks.lab,useRaster=F,
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5, border="azure4")
grid()
dev.off()

#germplasm map
tiff("./figures/g_samples_count.tif",
     res=300,pointsize=5,width=1500,height=1500*aspect,units="px",compression="lzw")
par(mar=c(2.5,2.5,1,1),cex=0.8, lwd=0.8)
plot(g_ras,col=cols,zlim=c(min(brks),max(brks)),useRaster=F, main="Genebank accessions",
     breaks=brks,lab.breaks=brks.lab,
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5, border="azure4")
grid()
dev.off()
