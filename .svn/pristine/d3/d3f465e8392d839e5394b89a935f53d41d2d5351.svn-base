#############################
# PCA -light version
# N. Castaneda - 2012
#############################

require(raster)
require(stats)

rm(list=ls()); g=gc(); rm(g)

wd <- "/curie_data2/ncastaneda/geodata/bio_2_5m"
setwd(wd)

#stk <- stack(paste("bio_",1:19,".asc",sep="")) #all variables
stk <- stack(paste("bio_",c(1,2,4:6,8:19),".asc",sep="")) #excluding BIO3 y BIO7

rs <- stk[[1]]
xy <- xyFromCell(rs,which(!is.na(rs[])))
xy <- as.data.frame(xy)

bios <- extract(stk,xy)
bios <- as.data.frame(bios)
bios <- cbind(xy,bios)

#biovars <- bios[,3:21] # When all 19 biovars are included
biovars <- bios[,3:19] # OJO REVISAR ANTES DE SEGUIR CON ESTE PASO!
x <- prcomp(biovars, scale=T)
xs <- summary(x)
y <- as.data.frame(predict(x,biovars))

dir.create("pca_result_raw")

#for (i in 1:19) {
for (i in 1:5) {
  cat("Assign",i,"\n")
  rs <- raster(stk)
  wcells <- cellFromXY(rs,data.frame(X=xy$x,Y=xy$y))
  rs[] <- NA
  rs[wcells] <- y[,paste("PC",i,sep="")]
  rs <- writeRaster(rs,paste("./pca_result_raw/pc_",i,".asc",sep=""),format='ascii')
}

write.csv(xs$importance,"./pca_result_raw/pca_importance.csv",quote=F,row.names=T)
