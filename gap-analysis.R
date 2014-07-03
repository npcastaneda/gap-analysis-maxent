#Gap analysis
#Based on Phaseolus study (Ramirez-Villegas et al., 2010)
#J. Ramirez
#CIAT
#March 2012

#-------------------------------------------------
# Run outside linux, only when new code is available
# cd /curie_data2/ncastaneda/code/gap-analysis-cwr/gap-analysis/gap-code
# cp * /curie_data2/ncastaneda/gap-analysis/gap_[crop_name]/_scripts
# cd /curie_data2/ncastaneda/gap-analysis
#-------------------------------------------------

stop("Warning: do not run the whole thing")

crop <- "avena" #change according to the coded name
#basic stuff - where is the code
src.dir <- paste("/curie_data2/ncastaneda/gap-analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!
gap.dir <-"/curie_data2/ncastaneda/gap-analysis" # !!! change accordingly !!!

#crop details
crop_dir <- paste(gap.dir,"/gap_",crop,sep="")

if (!file.exists(crop_dir)) {dir.create(crop_dir)}
setwd(crop_dir)

#basic stuff - creating folders
biomod <- paste(crop_dir,"/biomod_modeling",sep=""); if (!file.exists(biomod)) {dir.create(biomod)}

figs <- paste(crop_dir,"/figures",sep=""); if (!file.exists(figs)) {dir.create(figs)}

msks <- paste(crop_dir,"/masks",sep=""); if (!file.exists(msks)) {dir.create(msks)}
rm(msks)

narea <- paste(biomod,"/native-areas",sep=""); if (!file.exists(narea)) {dir.create(narea)}

# polys <- paste(narea,"/polyshps",sep=""); if (!file.exists(polys)) {dir.create(polys)}

ascis <- paste(narea,"/asciigrids",sep=""); if (!file.exists(ascis)) {dir.create(ascis)}

#== prepare taxonomic names for the analysis ==#
# source(paste(src.dir,"/000.fixTaxNames.R",sep=""))

#== compare germplasm vs. total records==#
source(paste(src.dir,"/01.countRecords.R",sep=""))

#== split H/G occurrences ==#
occ <- read.csv(paste("./occurrences/",crop,"_all.csv",sep=""))
occ <- occ[which(!is.na(occ$lat)),]
occ <- occ[which(!is.na(occ$lon)),]

h <- occ[which(occ$H==1),]
g <- occ[which(occ$G==1),]

write.csv(h,paste("./occurrences/",crop,"_h.csv",sep=""),quote=F,row.names=F)
write.csv(g,paste("./occurrences/",crop,"_g.csv",sep=""),quote=F,row.names=F)
write.csv(occ,paste("./occurrences/",crop,".csv",sep=""),quote=F,row.names=F)

#== samples densities comparison ==#
source(paste(src.dir,"/02.splitHG.R",sep=""))

#== prepare masks ==#
# source(paste(src.dir,"/000.prepareMasks.R",sep=""))
#set climate dir
env_dir <- "/curie_data2/ncastaneda/geodata/bio_2_5m" # !!! change accordingly !!!
# msks <- paste(crop_dir,"/masks",sep="")
# x <- createMasks(msks,env_dir)

#== crop climate data to extent of interest==#
eco_dir <- "/curie_data2/ncastaneda/geodata/wwf_eco_terr" # !!! change accordingly !!!
#source(paste(src.dir,"/000.ExtractVariables.R",sep=""))
#x <- maskVariables(crop_dir,env_dir,eco_dir)

#== create SWD occurrence files ==#
source(paste(src.dir,"/001.createSWD.R",sep=""))

occ_dir <- paste(crop_dir,"/occurrences",sep="")
swd_dir <- paste(crop_dir,"/swd",sep="")
if (!file.exists(swd_dir)) {dir.create(swd_dir)}

sample_file = paste(crop,".csv", sep="")

x <- extractClimates(input_dir=occ_dir,sample_file=sample_file,env_dir=env_dir,
                     env_prefix="bio_",env_ext="",lonfield="lon",
                     latfield="lat",taxfield="Taxon",output_dir=swd_dir)

#== preparing kernel density files ==#
# source(paste(src.dir,"/000.kernelDensity.R", sep="")) NEEDS TO BE FIXED!

#== splitting the occurrence files for biomod==# THIS NEEDS TIME!
# source(paste(src.dir,"/003.createOccurrenceFilesBiomod.R",sep=""))
# oDir <- paste(crop_dir,"/biomod_modeling/occurrence_files",sep="")
# if (!file.exists(oDir)) {dir.create(oDir)}
# x <- createOccFilesBio(occ=paste(crop_dir,"/swd/occurrences_swd_ok.csv",sep=""),
# taxfield="Taxon", outDir=oDir, env.dir=paste(crop_dir, "/biomod_modeling/current-clim", sep=""))

#== splitting the occurrence files ==#
source(paste(src.dir,"/003.createOccurrenceFiles.R",sep=""))
oDir <- paste(crop_dir,"/occurrence_files",sep="")
if (!file.exists(oDir)) {dir.create(oDir)}
x <- createOccFiles(occ=paste(crop_dir,"/swd/occurrences_swd_ok_kernel.csv",sep=""),
                    taxfield="Taxon", outDir=oDir)

#== prepare native areas ==#
source(paste(src.dir,"/004.createNARasters.R",sep=""))

#== making the pseudo-absences ==#
source(paste(src.dir,"/002.selectBackgroundArea.R",sep=""))
fList <- list.files("./occurrence_files",pattern=".csv")

MaxModDir <- paste(crop_dir,"/maxent_modeling",sep="")
if (!file.exists(MaxModDir)) {dir.create(MaxModDir)}

bkDir <- paste(crop_dir,"/maxent_modeling/background",sep="")
if (!file.exists(bkDir)) {dir.create(bkDir)}

for (f in fList) {
  cat("Processing",paste(f),"\n")
  iFile <- paste("./occurrence_files/",f,sep="")
  oFile <- paste("./maxent_modeling/background/",f,sep="")
  x <- selectBack(occFile=iFile, outBackName=oFile,
                  msk=paste(gap.dir,"/_backgroundFiles/backselection.asc",sep=""),
                  backFilesDir=paste(gap.dir,"/_backgroundFiles/",sep=""))
}

#== perform the maxent modelling in parallel ==#
source(paste(src.dir,"/005.modelingApproach.R",sep=""))
GapProcess(inputDir=crop_dir, OSys="linux", ncpu=3, j.size="-mx8192m")

#== summarise the models metrics ==#
source(paste(src.dir,"/006.summarizeMetricsThresholds.R",sep=""))
x <- summarizeMetrics(idir=crop_dir)

#== calculate area with SD<0.15 (aSD15) ==#
source(paste(src.dir,"/007.calcASD15.R",sep=""))
x <- summarizeASD15(idir=crop_dir)

#select which taxa are of use for species richness
#get the following modelling metrics:
# a. 25-fold average test AUC (ATAUC)
# b. 25-fold stdev of test AUC (STAUC)
# c. proportion of potential distribution with SD>15 (ASD15)

#== isValid==1 if ATAUC>0.7, STAUC<0.15, ASD15<10% ==#
acc <- read.csv(paste(crop_dir,"/maxent_modeling/summary-files/accuracy.csv",sep=""))
asd <- read.csv(paste(crop_dir,"/maxent_modeling/summary-files/ASD15.csv",sep=""))

for (spp in acc$SPID) {
  cat("Processing taxon",paste(spp),"\n")
  
  #getting the quality metrics
  atauc <- acc$TestAUC[which(acc$SPID==spp)]
  stauc <- acc$TestAUCSD[which(acc$SPID==spp)]
  asd15 <- asd$rateThresholded[which(asd$taxon==paste(spp))]
  
  #putting everything onto a row for appending
  row_res <- data.frame(Taxon=paste(spp),ATAUC=atauc,STAUC=stauc,ASD15=asd15,ValidModel=NA)
  
  #checking if any is na and correcting consequently
  if (is.na(atauc)) {atauc <- 0}
  if (is.na(stauc)) {stauc <- 1}
  if (is.na(asd15)) {asd15 <- 100}
  
  #reporting model quality
  if (atauc>=0.7 & stauc<=0.15 & asd15<=10) {
    row_res$ValidModel <- 1
  } else {
    row_res$ValidModel <- 0
  }
  
  #appending everything
  if (spp == acc$SPID[1]) {
    res_all <- row_res
  } else {
    res_all <- rbind(res_all,row_res)
  }
  
}
write.csv(res_all,paste(crop_dir,"/maxent_modeling/summary-files/modelsMets.csv",sep=""),quote=F,row.names=F)

#== create taxa for spp richness table ==#
table_base <- read.csv(paste(crop_dir,"/sample_counts/sample_count_table.csv",sep=""))
table_base <- data.frame(Taxon=table_base$TAXON)
table_base$IS_VALID <- NA

#== reading tables ==#
samples <- read.csv(paste(crop_dir,"/sample_counts/sample_count_table.csv",sep=""))
model_met <- read.csv(paste(crop_dir,"/maxent_modeling/summary-files/modelsMets.csv",sep=""))

names(model_met)[1]="TAXON"
table_base <- samples

for (spp in table_base$TAXON) {
  cat("Processing species",paste(spp),"\n")
  
  #modelling metrics
  if(sum(model_met$TAXON==paste(spp))==0){atauc <- NA
                                          stauc <- NA
                                          asd15 <- NA
                                          isval <- NA}else{
                                            atauc <- model_met$ATAUC[which(model_met$TAXON==paste(spp))]
                                            stauc <- model_met$STAUC[which(model_met$TAXON==paste(spp))]
                                            asd15 <- model_met$ASD15[which(model_met$TAXON==paste(spp))]
                                            isval <- model_met$ValidModel[which(model_met$TAXON==paste(spp))]}
  
  table_base$ATAUC[which(table_base$TAXON==paste(spp))] <- atauc
  table_base$STAUC[which(table_base$TAXON==paste(spp))] <- stauc
  table_base$ASD15[which(table_base$TAXON==paste(spp))] <- asd15
  table_base$IS_VALID[which(table_base$TAXON==paste(spp))] <- isval
  
}

table_base=table_base[c(1,11)]

table_base$IS_VALID[which(is.na(table_base$IS_VALID))]<-0

write.csv(table_base,paste(crop_dir,"/summary-files/taxaForRichness.csv",sep=""),row.names=F,quote=F)

rm(table_base, samples, model_met)

#== Filtering occurrences by native area ==#
source(paste(src.dir,"/000.filterOccFilesByNArea.R",sep=""))
x <- OccFilterNArea(crop_dir)

#== calculate size of distributional range ==#
source(paste(src.dir,"/008.sizeDR2.R",sep=""))
sizeDRProcess(inputDir=crop_dir, ncpu=2, crop=crop)

#== summarise area files ==#
source(paste(src.dir,"/008.summarizeDR.R",sep=""))
summarizeDR(crop_dir)

#== calculate environmental distance of distributional range ==#
source(paste(src.dir,"/009.edistDR.R",sep=""))
x <- summarizeDR_env(crop_dir)

#== calculate species richness ==#
source(paste(src.dir,"/010.speciesRichness2.R",sep=""))
x <- speciesRichness_alt(bdir=crop_dir)

#== create the priorities table ==#
source(paste(src.dir,"/000.prioritiesTable.R",sep=""))
x <- priTable(crop_dir)

#== calculate distance to populations ==#
source(paste(src.dir,"/011.distanceToPopulations2.R",sep=""))
ParProcess(crop_dir, ncpu=3)

#== calculate final gap richness ==#
source(paste(src.dir,"/012.gapRichness2.R",sep=""))
x <- gapRichness(bdir=crop_dir)

#== verify if gap map is available for each taxon ==#
priFile <- read.csv(paste(crop_dir, "/priorities/priorities.csv", sep=""))
newpriFile <- priFile
newpriFile$MAP_AVAILABLE <- NA

spList <- priFile$TAXON

for (spp in spList){
  fpcat <- priFile$FPCAT[which(priFile$TAXON==paste(spp))]
  if(file.exists(paste(crop_dir, "/gap_spp/", fpcat, "/", spp, ".asc.gz", sep=""))){
    newpriFile$MAP_AVAILABLE[which(newpriFile$TAXON==paste(spp))] <- 1
  } else {
    newpriFile$MAP_AVAILABLE[which(newpriFile$TAXON==paste(spp))] <- 0
  }
}

write.csv(newpriFile, paste(crop_dir, "/priorities/priorities.csv", sep=""), row.names=F, quote=F)
rm(newpriFile)

#== calculate distances to known populations ==#
source(paste(src.dir, "/000.maxpdistance.R", sep=""))
x <- maxpdist(crop_dir)

#== getting maps and figures ==#
source(paste(src.dir,"/013.mapsAndFigures.R",sep=""))

#== ensuring access to folders ==#
system(paste("chmod", "-R", "777", crop_dir))