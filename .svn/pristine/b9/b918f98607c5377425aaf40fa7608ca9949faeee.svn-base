# Prepare taxa name before the analysis. This code removes subsp. and var. from the taxonomic names
# April 2013

# Read occurrences
crop <- paste(toupper(substr(crop,1,1)),substr(crop,2,nchar(crop)),sep="")
occ <- read.csv(paste("./occurrences/",crop,"_all.csv",sep=""))
#occ <- read.csv("D:/CWR/_inputs/occurrences/Sorghum_all.csv")

names(occ)
# Erasing var. from text
tNames <- unique(occ$Taxon)
tNames<-as.vector(tNames)
fout <- tNames[grep("var._",tNames,fixed=T)]
occ$Taxon <- as.vector(occ$Taxon)
for (f in fout) {
  true_name <- strsplit(paste(f),"var._")[[1]]
  true_name3 <- paste(true_name[1],true_name[2],sep="")
  occ$Taxon[which(occ$Taxon==f)] <- array(true_name3,dim=length(which(occ$Taxon==f)))  
}

# Erasing subsp. from text
tNames <- unique(occ$Taxon)
tNames<-as.vector(tNames)
fout <- tNames[grep("subsp.",tNames,fixed=T)]
occ$Taxon <- as.vector(occ$Taxon)
for (f in fout) {
  true_name <- strsplit(paste(f),"subsp._")[[1]]
  true_name3 <- paste(true_name[1],true_name[2],sep="")
  occ$Taxon[which(occ$Taxon==f)] <- array(true_name3,dim=length(which(occ$Taxon==f)))  
}

crop <- tolower(crop)
write.csv(occ,paste("./occurrences/",crop,"_all.csv",sep=""),quote=F,row.names=F)
