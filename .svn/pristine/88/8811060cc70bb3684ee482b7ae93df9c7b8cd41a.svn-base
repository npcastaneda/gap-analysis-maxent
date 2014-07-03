#Extract the AUC evaluation statistics to select species with relatively high AUC (i.e. >.65 or .7)
#Test and train AUCs should be extracted from cross validated runs

cat(" \n")
cat("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n")
cat("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n")
cat("XXXXXXXX SUMMARIZE METRICS SCRIPT XXXXXXXXXX \n")
cat("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n")
cat("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \n")

cat(" \n")
cat(" \n")

summarizeMetrics <- function(idir) {
  
	#Setting the output directory
	odir <- paste(idir, "/maxent_modeling/summary-files", sep="")

	if (!file.exists(odir)) {
		dir.create(odir)
	}
	
	spList <- list.files(paste(idir, "/occurrence_files", sep=""))
	nspp <- nrow(spList)
	
	sppC <- 1
	sppCC <- 1
	for (spp in spList) {
		spp <- unlist(strsplit(spp, ".", fixed=T))[1]
		fdName <- spp #paste("sp-", spp, sep="")
		spFolder <- paste(idir, "/maxent_modeling/models/", fdName, sep="")
		
		#Performing only for existing folders
		if (file.exists(spFolder)) {
			
			cat("Processing species", spp, paste("...",round(sppCC/length(spList)*100,2),"%",sep=""), "\n")
			
			#Loading metrics files and adding one more field (SPID)
			metricsFile <- paste(spFolder, "/metrics/metrics.csv", sep="")
      
			#Creating Metrics and Threshold files
			Met <- data.frame(NSamples=rep(NA,1),TrainSamples=rep(NA,1),TestSamples=rep(NA,1),
			                  TrainAUC=rep(NA,1),TrainAUCSD=rep(NA,1),TestAUC=rep(NA,1),
			                  TestAUCSD=rep(NA,1),TrainR=rep(NA,1),TrainRSD=rep(NA,1),
			                  TestR=rep(NA,1),TestRSD=rep(NA,1),AllR=rep(NA,1),AllRSD=rep(NA,1),
			                  TrainLogD=rep(NA,1),TrainLogDSD=rep(NA,1),TestLogD=rep(NA,1),
			                  TestLogDSD=rep(NA,1),AllLogD=rep(NA,1),AllLogDSD=rep(NA,1),
			                  TrainRMSQD=rep(NA,1),TrainRMSQDSD=rep(NA,1),TestRMSQD=rep(NA,1),
			                  TestRMSQDSD=rep(NA,1),AllRMSQD=rep(NA,1),AllRMSQDSD=rep(NA,1))
			
			Thres <- data.frame(TenPercentile=rep(NA,1),TenPercentileSD=rep(NA,1),
			                    Prevalence=rep(NA,1),PrevalenceSD=rep(NA,1),FixedValue=rep(NA,1),
			                    MaxTrainSensSpec=rep(NA,1),MaxTrainSesnsSpecSD=rep(NA,1),
			                    EqualTrainSensSpec=rep(NA,1),EqualTrainSensSpecSD=rep(NA,1),
			                    BalanceTrainOmission=rep(NA,1),BalanceTrainOmissionSD=rep(NA,1),
			                    UpperLeftROC=rep(NA,1),UpperLeftROCSD=rep(NA,1))  
			
			write.csv(Met,paste(idir,"/maxent_modeling/models/metricsDummy.csv",sep=""),row.names=FALSE)
			write.csv(Thres,paste(idir,"/maxent_modeling/models/threshDummy.csv",sep=""),row.names=FALSE)
      
      dumMetFile <- paste(idir,"/maxent_modeling/models/metricsDummy.csv",sep="")
      dumThreshFile <- paste(idir,"/maxent_modeling/models/threshDummy.csv",sep="")
      if (file.exists(metricsFile)) {
			  metrics <- read.csv(metricsFile)
        threshFile <- paste(spFolder, "/metrics/thresholds.csv", sep="")
        thresholds <- read.csv(threshFile)
      } else {
        metrics <- read.csv(dumMetFile)
        thresholds <- read.csv(dumThreshFile)
      }
			
			#Adding one more field (SPID)
      metrics <- cbind(SPID=spp, metrics)
			thresholds <- cbind(SPID=spp, thresholds)
		  
			#Comprising everything onto a matrix
			if (sppC == 1) {
				finRes <- metrics
				finThr <- thresholds
				rm(thresholds)
				rm(metrics)
			} else {
				finRes <- rbind(finRes, metrics)
				finThr <- rbind(finThr, thresholds)
				rm(thresholds)
				rm(metrics)
			}
			
			sppC <- sppC + 1
		} else {
			cat("The species", spp, "was not modeled", paste("...",round(sppCC/length(spList)*100,2),"%",sep=""), "\n")
		}
		sppCC <- sppCC + 1
	}

	#Now writing the outputs
	cat("\n")
	cat("Writing outputs... \n")
	oFile <- paste(odir, "/accuracy.csv", sep="")
	write.csv(finRes, oFile, quote=F, row.names=F)
	
	oFile <- paste(odir, "/thresholds.csv", sep="")
	write.csv(finThr, oFile, quote=F, row.names=F)
	
	#Return the metrics data-frame
	return(finRes)
}
