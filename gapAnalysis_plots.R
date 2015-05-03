# Gap analysis plots
# H. Achicanoy
# CIAT, 2014
# Usar parte del código principal del Gap Analysis

# ============================================================================= #
# Arial font in R plots
# ============================================================================= #

install.packages("extrafont")
library("extrafont")
font_import()

fonts()

# ============================================================================= #
# Representatividad de muestreo
# ============================================================================= #

sampAll <- read.csv("./sample_counts/sample_count_table.csv")
# sampAll <- read.csv("F:/CIAT/CWR_project/Others/Potato_case/sample_count_table.csv")
sampAll <- sampAll[order(sampAll$TAXON),]
rownames(sampAll) <- 1:nrow(sampAll)

labels <- as.character(sampAll$TAXON)
sampAll$TAXON <- gsub(pattern="olanum_", replacement=". ", labels)
rm(labels)

fit <- lm(sampAll$GNUM~sampAll$TOTAL)
lims <- c(min(sampAll$TOTAL,sampAll$GNUM, na.rm=T),max(sampAll$TOTAL, na.rm=T))

View(sampAll[,c("TAXON", "TOTAL", "GNUM")])

tiff("./figures/genebank_vs_total_PAPER.tif", res=300, pointsize=8, width=1024, height=1024, units="px", compression="lzw")
# tiff("F:/CIAT/CWR_project/Others/Potato_case/genebank_vs_total_PAPER2.tif", family="Arial", res=300, pointsize=8, width=1024, height=1024, units="px", compression="lzw")

par(mar=c(5,5,1,1),cex=0.8)
plot(sampAll$TOTAL, sampAll$GNUM, pch=20, col="red", cex=1.7,
     xlim=c(lims[1],lims[2]+200), ylim=c(lims[1],lims[2]+50),
     xaxs="i", yaxs="i", xlab="", ylab="")
title(xlab=toupper("Total number of samples"),
      ylab=toupper("Number of genebank accessions"),
      col.lab="black", cex.lab=1.7, font.lab=2)
abline(0,1, lwd=0.75, lty=2)
abline(fit)
grid(lwd=0.75)
box()

etiq <- rep(NA, length(sampAll$TAXON))
etiq[1] <- sampAll$TAXON[1] # S. acaule
etiq[10] <- sampAll$TAXON[10] # S. boliviense
etiq[12] <- sampAll$TAXON[12] # S. brevicaule
etiq[17] <- sampAll$TAXON[17] # S. candolleanum
etiq[68] <- sampAll$TAXON[68] # S. stoloniferum

x = sampAll$TOTAL
y = sampAll$GNUM
y[10] <- y[10] - 350
x[12] <- x[12] - 450; y[12] <- y[12] - 350
x[68] <- x[68] - 400; y[68] <- y[68] - 350

with(sampAll, text(y~x, labels=etiq, cex=1.1, pos=3, col=1, font=4))

dev.off()

# ============================================================================= #
# Cobertura geográfica
# ============================================================================= #

prior <- read.csv(paste(crop_dir,"/maxent_modeling/summary-files/areas.csv",sep=""))
# prior <- read.csv("F:/CIAT/CWR_project/Others/Potato_case/areas.csv")
prior[,c("DRSize","GBSize")] <- prior[,c("DRSize","GBSize")]/1000
prior <- prior[order(prior$taxon),]
rownames(prior) <- 1:nrow(prior)
fit <- lm(prior$GBSize~prior$DRSize)
lims <- c(min(prior$DRSize,prior$GBSize,na.rm=T),
          max(prior$GBSize,prior$DRSize,na.rm=T)) # Problem 1

prior$taxon <- gsub(pattern="olanum_", replacement=". ", as.character(prior$taxon))

View(prior[,c("taxon", "DRSize", "GBSize")])

tiff(paste(crop_dir,"/figures/geographic_coverage_PAPER.tif",sep=""), res=300, pointsize=12, width=1500, height=1500, units="px", compression="lzw")
# tiff("F:/CIAT/CWR_project/Others/Potato_case/geographic_coverage_PAPER2.tif", family="Arial", res=300, pointsize=12, width=1500, height=1500, units="px", compression="lzw")

par(mar=c(5,5,1,1),cex=0.8)
plot(prior$DRSize, prior$GBSize, pch=20, cex=1.7, col="red",
     xlim=c(lims[1],lims[2]+200), ylim=c(0,1000),
     xaxs="i", yaxs="i", xlab="", ylab="")
title(xlab=toupper("Potential distribution coverage (sq-km * 1000)"),
      ylab=toupper("Genebank accessions CA50 (sq-km * 1000)"),
      col.lab="black", cex.lab=1.2, font.lab=2)

abline(0,1,lwd=0.75,lty=2)
abline(fit) # Problem 2
grid(lwd=0.75)
box()

etiq <- rep(NA, length(prior$taxon))

etiq[1] <- prior$taxon[1]
etiq[12] <- prior$taxon[12]
etiq[17] <- prior$taxon[17]
etiq[19] <- prior$taxon[19]
etiq[24] <- prior$taxon[24]
etiq[25] <- prior$taxon[25]
etiq[67] <- prior$taxon[67]

x = prior$DRSize
y = prior$GBSize

x[19] <- x[19] - 100 # S. chacoense
x[24] <- x[24] + 270; y[24] <- y[24] - 25 # S. colombianum

with(prior, text(y~x, labels=etiq, cex=1.1, pos=3, col=1, font=4))

dev.off()

# ============================================================================= #
# Environmental representativeness score plots
# ============================================================================= #

edist_info <- read.csv(paste(crop_dir,"/maxent_modeling/summary-files/edist_wwf.csv",sep=""),header=T)
# edist_info <- read.csv("F:/CIAT/CWR_project/Others/Potato_case/edist_wwf.csv",header=T)
edist_info <- edist_info[order(edist_info$taxon),]
rownames(edist_info) <- 1:nrow(edist_info)
edist_info$taxon <- gsub(pattern="olanum_", replacement=". ", edist_info$taxon)

## Numerador
# GBDist.PC1

## Denominador
# DRDist.PC1 # Modelo
# CHDist.PC1 # Convex hull

den <- edist_info$DRDist.PC1
den[which(den==0)] <- edist_info$CHDist.PC1[which(den==0)]

num <- edist_info$GBDist.PC1

mdata <- data.frame(Taxon=edist_info$taxon,
                    Num=num,
                    Den=den)

mdata$NArea <- edist_info$NADist.PC1
mdata$Herbm <- edist_info$HBDist.PC1

View(mdata[,c("Taxon","Den","Num")])

tiff("./figures/environmental_coverage_PAPER.tif", res=300, pointsize=8, width=1500, height=1500, units="px", compression="lzw")
# tiff("F:/CIAT/CWR_project/Others/Potato_case/environmental_coverage_PAPER2.tif", family="Arial", res=300, pointsize=8, width=1500, height=1500, units="px", compression="lzw")

par(mar=c(6.5,6.5,1,1),cex=0.8)
plot(x=mdata$Den, y=mdata$Num,
     pch=20, col=2, cex=2.3,
     xlim=c(0, 60), ylim=c(0, 40),
     xaxs="i", yaxs="i",
     xlab="", ylab="")
grid()
title(#xlab="",
  ylab=toupper("Number of ecosystems represented in\ngenebank accessions CA50"),
  col.lab="black", cex.lab=2.0, font.lab=2)
mtext(toupper("Number of ecosystems represented in\npotential distribution coverage"),
      side=1, line=4.3, cex=1.5, font=2)
abline(0,1, lwd=0.75, lty=2)
fit <- lm(Num ~ Den, mdata)
abline(fit)
box()

etiq <- rep(NA, length(mdata$Taxon))
etiq[14] <- as.character(mdata$Taxon[14])
etiq[19] <- as.character(mdata$Taxon[19])
etiq[24] <- as.character(mdata$Taxon[24])
etiq[67] <- as.character(mdata$Taxon[67])

x = mdata$Den
y = mdata$Num

x[19] <- x[19] + 7; y[19] <- y[19] - 1

with(mdata, text(y~x, labels=etiq, cex=1.5, pos=3, col=1, font=4))

dev.off()
