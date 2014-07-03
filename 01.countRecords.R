#Gap analysis workshop
#March 2012
#stop("")

#read occurrences
occ <- read.csv(paste("./occurrences/",crop,"_all.csv",sep=""))
#gpSpp <- read.csv("/curie_data2/ncastaneda/gap-analysis/PsppofPcrops_template.csv")

#prepare files - OJO CORRER PASO POR PASO Y VERIFICAR NUMERO DE REGISTROS CON nrow(occ)
occ$H <- ifelse(occ$Type =="H",1,0)
occ$G <- ifelse(occ$Type =="G",1,0)
if(sum(occ$final_cult_stat != "cultivated",na.rm=T)!=0){
occ <- occ[which(occ$final_cult_stat != "cultivated"),]}
if(sum(occ$final_origin_stat != "non_native",na.rm=T)!=0){
  occ <- occ[which(occ$final_origin_stat != "non_native"),]}
occ$is_hybrid[which(is.na(occ$is_hybrid))] <- 0
if(sum(occ$is_hybrid == 1, na.rm = T)!=0){
  occ <- occ[which(occ$is_hybrid != 1),]}

#template <- gpSpp$Taxon_name
#occ <- occ[occ$Taxon %in% template,] #Parece que si

#occ <- read.csv(paste("./occurrences/",crop,"_all.csv",sep=""))

if(sum(names(occ)=="taxon")==0){
  pos <- which(names(occ)=="taxon")
  names(occ)[pos] <- "Taxon"
}else{
  cat("No need to change taxon field name")
}

taxField <- "Taxon"
hField <- "H"
gField <- "G"

write.csv(occ,paste("./occurrences/",crop,"_all.csv",sep=""),quote=F,row.names=F)

#taxon unique values
taxNames <- unique(occ[,taxField])

for (tax in taxNames) {
  cat("Counting",paste(tax),"\n")
  
  #subselect the data for this taxon
  taxData <- occ[which(occ[,taxField]==tax),]
  taxDataQ <- data.frame(paste(taxData[,taxField]),taxData[,"lon"],taxData[,"lat"],taxData[,hField],taxData[,gField])
  names(taxDataQ) <- c(taxField,"lon","lat",hField,gField)
  
  allData <- taxDataQ; allData[,hField] <- NULL; allData[,gField] <- NULL
  hData <- taxDataQ[which(taxDataQ[,hField]==1),]
  hData[,hField] <- NULL; hData[,gField] <- NULL
  gData <- taxDataQ[which(taxDataQ[,gField]==1),]
  gData[,hField] <- NULL; gData[,gField] <- NULL
  
  #count herbarium recs and germplasm samples (totals regardless of populations)
  hc <- sum(taxData[,hField])
  gc <- sum(taxData[,gField])
  total <- nrow(taxData)
  
  #count h and g samples (totals but only considering populations -unique coordinates)
  hc_u <- nrow(unique(hData)[complete.cases(unique(hData)),]) # New addition
  gc_u <- nrow(unique(gData)[complete.cases(unique(gData)),]) # New addition
  total_u <- nrow(unique(allData)[complete.cases(unique(allData)),]) # New addition
  
  row_out <- data.frame(TAXON=paste(tax),HNUM=hc,GNUM=gc,HNUM_RP=hc_u,GNUM_RP=gc_u,TOTAL_RP=total_u)
  
  if (tax==taxNames[1]) {
    sampAll <- row_out
  } else {
    sampAll <- rbind(sampAll,row_out)
  }
}

if (!file.exists("./sample_counts")) {dir.create("./sample_counts")}
sampAll$TOTAL <- sampAll$HNUM+sampAll$GNUM
write.csv(sampAll,"./sample_counts/sample_count_table.csv",row.names=F,quote=F)

sampAll <- read.csv("./sample_counts/sample_count_table.csv")

#making the plot.
#1. do a regression between TOTAL(x) and GNUM(y)
fit <- lm(sampAll$GNUM~sampAll$TOTAL)

lims <- c(min(sampAll$TOTAL,sampAll$GNUM),max(sampAll$TOTAL))

#do the plot
tiff("./figures/genebank_vs_total.tif",
         res=300,pointsize=8, width=1000,height=1000,units="px",compression="lzw")
par(mar=c(5,5,1,1),cex=0.8)
plot(sampAll$TOTAL,sampAll$GNUM,pch=20, col="red",cex=1,xlim=lims,ylim=lims,
     xlab="Total number of samples",
     ylab="Number of genebank accessions")
abline(0,1,lwd=0.75,lty=2)
#abline(h=500,lwd=0.75,lty=1,col="red")
#abline(h=100,lwd=0.75,lty=2,col="red")
lines(sampAll$TOTAL,fit$fitted.values)
grid(lwd=0.75)

# NOTE: Personalize this according to the crop you're working with!
#
#text(sampAll$TOTAL[which(sampAll$TAXON=="Solanum_pimpinellifolium")]-45,
#     sampAll$GNUM[which(sampAll$TAXON=="Solanum_pimpinellifolium")]-15,
#     "S. pimpinellifolium",cex=0.55, font=3)

dev.off()
