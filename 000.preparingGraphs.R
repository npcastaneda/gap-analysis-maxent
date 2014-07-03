#plot the CA50 vs Potential coverage thing
prior <- read.csv(paste(crop_dir,"/priorities/priorities.csv",sep=""))

fit <- lm(prior$CA50_G~prior$PD_COV)
lims <- c(min(prior$PD_COV,prior$CA50_G),max(prior$CA50_G,prior$PD_COV))/1000

#do the plot
tiff(paste(crop_dir,"/figures/geographic_coverage.tif",sep=""),
     res=300,pointsize=12,width=1500,height=1000,units="px",compression="lzw")
par(mar=c(5,5,1,1),cex=0.8)
plot(prior$PD_COV/1000,prior$CA50_G/1000,pch=20,cex=0.75,xlim=lims,ylim=c(0,1000),
     xlab="Potential geographic coverage (sq-km * 1000)",
     ylab="Genebank accessions CA50 (sq-km * 1000")
abline(0,1,lwd=0.75,lty=2)
lines(prior$PD_COV/1000,fit$fitted.values/1000)
grid(lwd=0.75)
#text(prior$PD_COV[which(prior$Taxon=="Lathyrus_sativus")]/1000+2000,
#     prior$CA50_G[which(prior$Taxon=="Lathyrus_sativus")]/1000,
#     "L. sativus",cex=0.5)
dev.off()

#plot the gap richness maps, uncertainty and related stuff
source(paste(src.dir,"/000.zipRead.R",sep=""))

gap_rich <- zipRead(paste(crop_dir,"/gap_richness/",sep=""),"gap-richness.asc.gz")
gap_dpmax <- zipRead(paste(crop_dir,"/gap_richness/",sep=""),"gap-richness-dpmax.asc.gz")
gap_sdmax <- zipRead(paste(crop_dir,"/gap_richness/",sep=""),"gap-richness-sdmax.asc.gz")

library(maptools); data(wrld_simpl)

z <- extent(gap_rich)
aspect <- (z@ymax-z@ymin)*1.2/(z@xmax-z@xmin)

grich_brks <- unique(gap_rich[!is.na(gap_rich[])])
grich_cols <- c("grey 80",colorRampPalette(c("yellow","orange","red"))(length(grich_brks)-2))

#gap richness map
tiff(paste(crop_dir,"/figures/gap_richness.tif",sep=""),
     res=300,pointsize=7,width=1500,height=1500*aspect,units="px",compression="lzw")
par(mar=c(2.5,2.5,1,1),cex=0.8)
plot(gap_rich,col=grich_cols,zlim=c(min(grich_brks),max(grich_brks)),useRaster=F,
     breaks=grich_brks,lab.breaks=grich_brks,
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5)
grid()
dev.off()


###############
#gsd_brks <- unique(quantile(gap_sdmax[],probs=seq(0,1,by=0.05),na.rm=T))
gap_sdmax[which(gap_rich[]==0)] <- NA
gsd_brks <- c(seq(0,max(gap_sdmax[],na.rm=T),by=0.05),max(gap_sdmax[],na.rm=T))
gsd_cols <- colorRampPalette(c("light green","green","light blue","blue"))(length(gsd_brks)-1)
gsd_labs <- round(gsd_brks,2)

#gap uncertainty map (standard deviation)
tiff(paste(crop_dir,"/figures/gap_richness_sd.tif",sep=""),
     res=300,pointsize=7,width=1500,height=1500*aspect,units="px",compression="lzw")
par(mar=c(2.5,2.5,1,1),cex=0.8)
plot(gap_sdmax,col=gsd_cols,zlim=c(min(gsd_brks),max(gsd_brks)),useRaster=F,
     breaks=gsd_brks,lab.breaks=gsd_labs,
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5)
grid()
dev.off()


gap_dpmax[which(gap_rich[]==0)] <- NA
gdp_brks <- unique(quantile(gap_dpmax[],probs=seq(0,1,by=0.05),na.rm=T))
gdp_cols <- colorRampPalette(c("yellow","green","blue"))(length(gdp_brks)-1)
gdp_labs <- round(gdp_brks,2)


#gap uncertainty map (popdist)
tiff(paste(crop_dir,"/figures/gap_richness_dp.tif",sep=""),
     res=300,pointsize=7,width=1500,height=1500*aspect,units="px",compression="lzw")
par(mar=c(2.5,2.5,1,1),cex=0.8)
plot(gap_dpmax,col=gdp_cols,zlim=c(min(gdp_brks),max(gdp_brks)),useRaster=F,
     breaks=gdp_brks,lab.breaks=gdp_labs,
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5)
grid()
dev.off()


#comparison of expert and gap scores
eps <- read.csv(paste(crop_dir,"/priorities/expert_gap_comparison.csv",sep=""))
spear <- cor(eps$FPS,eps$EPS,method="spearman")
eps$RD <- (eps$FPS-eps$EPS)*10

fit <- lm(eps$EPS~eps$FPS)


tiff(paste(crop_dir,"/figures/expert_evaluation.tif",sep=""),
     res=300,pointsize=12,width=1500,height=1000,units="px",compression="lzw")
par(mar=c(5,5,1,1),cex=0.8)
plot(eps$FPS,eps$EPS,xlab="Gap Analysis Final priority score",
     ylab="Expert priority score",pch=20,xlim=c(0,8),ylim=c(0,8))
lines(eps$FPS,fit$fitted.values)
abline(0,1,lty=2)
grid()
dev.off()


tiff(paste(crop_dir,"/figures/expert_evaluation_RD.tif",sep=""),
     res=300,pointsize=12,width=1500,height=1000,units="px",compression="lzw")
par(mar=c(5,5,1,1),cex=0.8)
hist(eps$RD,xlab="Relative difference (%)",
     ylab="Frequency (number of taxa)",
     breaks=20,xlim=c(-100,100),col="grey 70",main=NA)
abline(v=0,col="red")
grid()
dev.off()


#plot species richness
source(paste(src.dir,"/000.zipRead.R",sep=""))
sp_rich <- zipRead(paste(crop_dir,"/species_richness",sep=""),"species-richness.asc.gz")
sdmax <- zipRead(paste(crop_dir,"/species_richness",sep=""),"species-richness-sdmax.asc.gz")


library(maptools); data(wrld_simpl)

z <- extent(sp_rich)
aspect <- (z@ymax-z@ymin)*1.2/(z@xmax-z@xmin)

rich_brks <- unique(sp_rich[!is.na(sp_rich[])])
rich_cols <- c("grey 80",colorRampPalette(c("yellow","orange","red"))(length(rich_brks)-2))

#gap richness map
tiff(paste(crop_dir,"/figures/sp_richness.tif",sep=""),
     res=300,pointsize=7,width=1500,height=1500*aspect,units="px",compression="lzw")
par(mar=c(2.5,2.5,1,1),cex=0.8)
plot(sp_rich,col=rich_cols,zlim=c(min(rich_brks),max(rich_brks)),useRaster=F,
     breaks=rich_brks,lab.breaks=rich_brks,
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5)
grid()
dev.off()


###############
#gsd_brks <- unique(quantile(gap_sdmax[],probs=seq(0,1,by=0.05),na.rm=T))
sdmax[which(sp_rich[]==0)] <- NA
sd_brks <- c(seq(0,max(sdmax[],na.rm=T),by=0.05),max(sdmax[],na.rm=T))
sd_cols <- colorRampPalette(c("light green","green","light blue","blue"))(length(sd_brks)-1)
sd_labs <- round(sd_brks,2)


#gap uncertainty map (standard deviation)
tiff(paste(crop_dir,"/figures/sp_richness_sd.tif",sep=""),
     res=300,pointsize=7,width=1500,height=1500*aspect,units="px",compression="lzw")
par(mar=c(2.5,2.5,1,1),cex=0.8)
plot(gap_sdmax,col=sd_cols,zlim=c(min(sd_brks),max(sd_brks)),useRaster=F,
     breaks=sd_brks,lab.breaks=sd_labs,
     horizontal=T,
     legend.width=1,
     legend.shrink=0.99)
plot(wrld_simpl,add=T,lwd=0.5)
grid()
dev.off()

