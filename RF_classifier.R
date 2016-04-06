#!/usr/bin/env Rscript
library(randomForest)

# perform the classification

args <- commandArgs(trailingOnly = TRUE)

# only for testing
# args <- c("data/NKTCL.Jiang.L01.snp.tab")

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("Usage: Rscript RF_classifier.R <model_dir> <in.tab> <tum_num> [out.prefix]", call.=FALSE)
} else if (length(args) == 3) {
  # default output file
  args[4] <- args[2]
}

data <- read.table(args[2],header=T)

# clean up the data
snp_mut_type <- c('nonsynonymous_SNV','splicing',
                  'stopgain','stoploss')

# output predicted the input variants with un-recognized mutation types
outdata <- data[!data$MutType %in% snp_mut_type,]
out <- paste0(args[4],".UnknownMutationType.txt")
write.table(outdata,out,row.names = F, sep = "\t" )
data <- data[data$MutType %in% snp_mut_type,]

# output predicted the input variants with VAF < 0.05 or low Var read depth 
outdata <- data[!((data$VAF >= 0.05 & data$VarDP >=10) | data$VAF >= 0.2),]
out <- paste0(args[4],".LowQualityMutations.txt")
write.table(outdata,out,row.names = F, sep = "\t" )
data <- data[(data$VAF >= 0.05 & data$VarDP >=10) | data$VAF >= 0.2,]

# output predicted the input common SNPs
outdata <- data[data$snp138 == 1 & data$PopFreqMax >= 0.01,]
out <- paste0(args[4],".CommonSNPs.txt")
write.table(outdata,out,row.names = F, sep = "\t" )

# remove common SNPs with freq>0.01
data <- data[data$snp138 == 0 | data$PopFreqMax < 0.01,]

# Correcting the missing data with logical ways
attach(data)
for (i in (1:nrow(data))){
  if (MutType[i] %in% snp_mut_type[3:4]){
    data$LOFMut[i] <- 1
  }else{
    data$LOFMut[i] <- 0
  }
  
  # Normalize the values into 0~1 range
  # The min and max values were based on a randomly sampling 
  # of 1,000,000 records from the ANNOVAR database file
  data$SIFT[i] <- 1 - SIFT[i]
  data$FATHMM[i] <- 1 - (FATHMM[i] + 12.45) / (9.95 + 12.45)
  data$GERP[i] <- (GERP[i] + 12.3) / (6.17 + 12.3)
  data$MutationAssessor[i] <- (MutationAssessor[i] + 5.545) / (5.975 + 5.545)
  data$MetaSVM[i] <- (MetaSVM[i] - 1.629) / (1.557 + 1.629)
  data$phyloP7way_vertebrate[i] <- (phyloP7way_vertebrate[i] + 4.552) / (1.602 + 4.552)
  data$phyloP20way_mammalian[i] <- (phyloP20way_mammalian[i] + 9.336) / (1.199 + 9.336)
  data$CADD_phred[i] <- CADD_phred[i] / 60
  data$integrated_fitCons[i] <- integrated_fitCons[i] / 0.84
  data$SiPhy_29way_logOdds[i] <- SiPhy_29way_logOdds[i] / 30.947
  data$PROVEAN[i] <- 1 - (data$PROVEAN[i] + 14) / (13.15 + 14)
}
detach(data)
nlof <- nrow(data[which(data$LOFMut==1),])
for (score in 
     c('SIFT', 'Polyphen2_HDIV', 'MutationAssessor', 'FATHMM',
       'PROVEAN','VEST3','MetaSVM','MetaLR')){
  data[which(data$LOFMut==1),score] <- rep(max(data[,score]),nlof)
}
med_score <- c()
for (score in 
     c("SIFT","Polyphen2_HDIV","LRT","MutationTaster",
       "MutationAssessor","FATHMM","PROVEAN","VEST3",
       "CADD_phred","DANN","fathmm_MKL_coding","MetaSVM",
       "MetaLR","integrated_fitCons","GERP",
       "phyloP7way_vertebrate","phyloP20way_mammalian",
       "phastCons7way_vertebrate","phastCons20way_mammalian",
       "SiPhy_29way_logOdds","dbscSNV_ADA","dbscSNV_RF")){
  med_score[score] <- median(data[,score],na.rm = T)
  for (i in 1:nrow(data)){
    if(is.na(data[i,score])) data[i,score] <- med_score[score]
  }
}

for (i in 1:nrow(data)){
  # if the 40M baf info is missing, firstly use the info of the other
  # tumors as they always correlate well; if still NA, use the overall
  # median value 0.275
  if(is.na(data$BAF_s40M[i])) data$BAF_s40M[i] <- data$T2_BAF_s40M[i]
  if(is.na(data$BAF_s40M[i])) data$BAF_s40M[i] <- 0.275
  if(is.na(data$T2_BAF_s40M[i])) data$T2_BAF_s40M[i] <- data$BAF_s40M[i]

  # Caculating binominal p-values and adding a very small number (1e-100) to prevent log errors for p=0
  data$logBinomP[i] <- 
    -log10(binom.test(data$VarDP[i],
                      (data$RefDP[i]+data$VarDP[i]))$p.value + 1e-100) 
  if(data$logBinomP[i] > 100) data$logBinomP[i] <- 100
  if(!is.na(data$T2_RefDP[i] + data$T2_VarDP[i])){
    data$T2_logBinomP[i] <-
      -log10(binom.test(data$T2_VarDP[i],
                        (data$T2_RefDP[i]+data$T2_VarDP[i]))$p.value + 1e-100)
    if(data$T2_logBinomP[i] > 100) data$T2_logBinomP[i] <- 100
    data$T2_VAFdiff[i] <- data$VAF[i] - data$T2_VAF[i]
  }else{
	  data$T2_logBinomP[i] <- NA
	  data$T2_VAFdiff[i] <- NA
  }
}  

# cut the VAF in half if surrounding-40M mean BAF of
# common SNP < 0.15 (LOH)
data$aVAF <- data$VAF
for (i in 1:nrow(data)){
  if (is.na(data$T2_VAF[i])){
    data$T2_aVAF[i] <- NA
  }else{
    data$T2_aVAF[i] <- data$T2_VAF[i]
  }
  if(data$BAF_s40M[i] <= 0.15){
    data$aVAF[i] <- data$VAF[i] / 2
  }
  if(!is.na(data$T2_VAF[i])){
    if(data$T2_BAF_s40M[i] <= 0.15){
      data$T2_aVAF[i] <- data$T2_VAF[i] / 2
    }
  }
}
if (args[3]==1){
	feat <- c(
	  'PopFreqMax','snp138','SIFT','Polyphen2_HDIV',
	  'MutationAssessor','PROVEAN','VEST3','CADD_phred',
	  'DANN','fathmm_MKL_coding','MetaSVM','MetaLR','GERP',
	  'phyloP20way_mammalian', 'phastCons7way_vertebrate',
	  'phastCons20way_mammalian','SiPhy_29way_logOdds',
	  'VarDP','VAF','BAF_s40M','LOFMut','logBinomP','aVAF'
	  )
}else{
	feat <- c(
      'PopFreqMax','snp138','SIFT','Polyphen2_HDIV',
	  'MutationAssessor','PROVEAN','VEST3','CADD_phred',
	  'DANN','fathmm_MKL_coding','MetaSVM','MetaLR','GERP',
	  'phyloP20way_mammalian', 'phastCons7way_vertebrate',
	  'phastCons20way_mammalian','SiPhy_29way_logOdds',
	  'VarDP','VAF','BAF_s40M','T2_VarDP','T2_VAF',
	  'T2_BAF_s40M','LOFMut','logBinomP','T2_logBinomP',
	  'T2_VAFdiff','aVAF','T2_aVAF'
	  )
	sg_d <- data[is.na(data$T2_VAF),]  # for sites with missing T2 inforation, use the single-tumor model
	data <- data[!is.na(data$T2_VAF),]	
}

## RF
td <- data[,feat]
require(randomForest)
# load the model
if(args[3]==1){
	rf_model <- paste0(args[1],"/RF.single.rda")
}else{
	rf_model <- paste0(args[1],"/RF.multiple.rda")
}
if(file.exists(rf_model)){
  load(rf_model)
}else{
  stop(paste0("Missing model file: ",rf_model), call.=FALSE)
}
somatic_cutoff <- 0.5  # the default value
td_prob <- predict(x.rf, type="prob", newdata=td)
data$somatic_prob <- td_prob[,2]
if(args[3] > 1){
	feat <- c(
	  'PopFreqMax','snp138','SIFT','Polyphen2_HDIV',
	  'MutationAssessor','PROVEAN','VEST3','CADD_phred',
	  'DANN','fathmm_MKL_coding','MetaSVM','MetaLR','GERP',
	  'phyloP20way_mammalian', 'phastCons7way_vertebrate',
	  'phastCons20way_mammalian','SiPhy_29way_logOdds',
	  'VarDP','VAF','BAF_s40M','LOFMut','logBinomP','aVAF'
	  )
	rf_model <- paste0(args[1],"/RF.single.rda")
	td <- sg_d[,feat]
	if(file.exists(rf_model)){
	  load(rf_model)
	}else{
	  stop(paste0("Missing model file: ",rf_model), call.=FALSE)
	}
	td_prob <- predict(x.rf, type="prob", newdata=td)
	sg_d$somatic_prob <- td_prob[,2]
	data <- rbind(data,sg_d)
}

# output predicted somatic variants
out <- paste0(args[4],".somatic.txt")
write.table(data[data$somatic_prob > somatic_cutoff,],out,row.names = F, sep = "\t" )

# output predicted non-somatic variants
out <- paste0(args[4],".non-somatic.txt")
write.table(data[data$somatic_prob <= somatic_cutoff,],out,row.names = F, sep = "\t" )



# test if an evaluation is needed
nvaf_na <- summary(is.na(data$N_VAF))
if(is.na(nvaf_na['TRUE'])){
  nvaf_count <- nrow(data)
}else{
  nvaf_count <- nrow(data) - is.na(nvaf_na['TRUE'])
}

# evalutaton dataset
if(nvaf_count >= 10){
  data$Prediction <- data$somatic_prob > somatic_cutoff
  evd <- data[!is.na(data$N_VAF),]
  evd <- evd[evd$N_RefDP + evd$N_VarDP >= 10,]
}
nvaf_count <- nrow(evd)
if (nvaf_count >= 10){
  for (i in 1:nrow(evd)){
    if(evd$N_VAF[i] < 0.02){
      evd$Somatic[i] <- 1
    }else{
      evd$Somatic[i] <- 0
    }
  }
  evd$Somatic <- factor(evd$Somatic)
  evd.tab <- table(evd$Somatic,evd$Prediction)
  out <- paste0(args[4],".evaluation.txt")
  perf_rep <- rbind(c(args[4], evd.tab[1,1], evd.tab[1,2], 
                      evd.tab[2,1], evd.tab[2,2]))
  perf_rep <- as.data.frame(perf_rep)
  colnames(perf_rep) <- c('ID','TN','FP','FN','TP')
  write.table(perf_rep, out, row.names = F, sep = "\t")
}

# Saving the image for future inquiry 
# save.image("Rdata/somaticID.v01.rdata")


