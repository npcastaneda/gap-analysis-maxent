# J. Ramirez 
# 2012

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

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_sativus")]/1000+2000,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_sativus")]/1000,
     "L. sativus",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_cicera")]/1000+1500,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_cicera")]/1000,
     "L. cicera",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_aphaca")]/1000+2000,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_aphaca")]/1000,
     "L. aphaca",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_pratensis")]/1000-2000,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_pratensis")]/1000,
     "L. pratensis",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_ochrus")]/1000+1500,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_ochrus")]/1000,
     "L. ochrus",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_clymenum")]/1000+2000,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_clymenum")]/1000,
     "L. clymenum",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_inconspicuus")]/1000+2500,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_inconspicuus")]/1000,
     "L. inconspicuus",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_annuus")]/1000+1500,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_annuus")]/1000,
     "L. annuus",cex=0.5)

text(prior$PD_COV[which(prior$Taxon=="Lathyrus_pseudocicera")]/1000+2400,
     prior$CA50_G[which(prior$Taxon=="Lathyrus_pseudocicera")]/1000-10,
     "L. pseudocicera",cex=0.5)

dev.off()
