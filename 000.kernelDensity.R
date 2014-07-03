##################################################
## PCA + graphs                                 ##
## Author: N. Castaneda                         ##
## Date: 2012-07-06                             ##
## Comments: this version calculates the PCA per##
## genepool! (good version)                     ##
##################################################

require(adehabitat)
require(plyr)
require(sp)

#------------------------------------------------------
# Preparing data (removing columns lat, lon and Taxon)
#------------------------------------------------------

crop.dir <- "D:/CWR-collaborations/potatoCIP/_process" #TEST!!!!
setwd(crop.dir)

taxfield <- "Taxon"
data.pca <- read.csv("./swd/occurrences_swd_ok_kernel.csv")

# Saving taxfield
tax <- data.pca[,taxfield]

# Selecting columns
data.pca <- data.pca[-c(1,2,3)]

#------------------------------------------------------
# PCA - using singular value decomposition (SVD)
#------------------------------------------------------
# !!!!TO DO!!!!: include conditional for cases where PCA can't be calculated

prc <- prcomp(data.pca, center=TRUE, scale=TRUE)  

# Get the scores
pca.scores <- (prc$x)
pca.scores <- as.data.frame(pca.scores)

# Add column with taxa name
pca.scores <- cbind(Taxon=tax, pca.scores)

# Saving the files  
pca.dir <- "./pca_scores"
if (!file.exists(pca.dir)) {dir.create(pca.dir)}
write.csv(pca.scores, "./pca_scores/pca_scores.csv", quote=F, row.names=F)

#------------------------------------------------------
# Prepare data for graphs
#------------------------------------------------------
# http://genetic-codes.blogspot.co.uk/2010/05/kernel-density-plots.html
scores <- read.csv("./pca_scores/pca_scores.csv")

countRecords <- count(scores, "Taxon")

plot.new()

for(spp in countRecords$Taxon){
  cat("Processing", spp, "\n")
  spdata <- scores[scores$Taxon == spp,]
  occ <- scores
  ch <- occ[chull(cbind(occ$PC1,occ$PC2)),2:3]
  ch <- rbind(ch, ch[1,])
  cat("Transforming to polygons \n")
  pol <- SpatialPolygons(list(Polygons(list(Polygon(ch)), 1)))
  plot(pol, add=T)
}



plot(pol)




X <- matrix(stats::rnorm(2000), ncol = 2)
chull(X)
## Not run: 
# Example usage from graphics package
plot(X, cex = 0.5)
hpts <- chull(X)
lines(hpts)
hpts <- c(hpts, hpts[1])
lines(X[hpts, ])
is(hpts)




xy <- scores[,c("PC1","PC2")]
id <- scores[,c("Taxon")]
colours <- rainbow(length(unique(id)))
means <- aggregate(xy, list(id), mean)

## Run the kernel density analysis, this uses the least-square cross-validation (LSCV) 
## and 70% threshhold

(hr <- kernelUD(xy, id))
image(hr)
ver <- getverticeshr(hr, lev=65)

layout(cbind(1, 2), width=c(8, 3))  # put legend on bottom 1/8th of the chart
par(xpd=FALSE)
#par(mfrow=c(1,2))
par(mar=c(2.5,2.5,1,1),cex=0.8,lwd=1)

plot(ver, colborder=colours, colpol="NA", sub="", lwd=2.5,
     xlab="PC 1", ylab="PC 2", cex.axis=0.8)

abline(h=0,lty=2,col="grey")
abline(v=0,lty=2,col="grey")

# plot.new()

par(xpd=TRUE)
#par(mar=c(1, 1, 1, 1)+0.1, xpd=TRUE)
#par(mar=c(2.5,2.5,1,1),cex=0.8,lwd=1)


legend("topright", names(ver),lty=c(1,1),lwd=c(3,3),col=colours, cex=0.7)

