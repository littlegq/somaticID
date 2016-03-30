#!/usr/bin/Rscript
library(randomForest)

# Train novel models using data with matched normal

files <- commandArgs(trailingOnly = TRUE)
if (length(files) < 1){
  stop("Usage: Rscript RF_trainer.R <file1.tab> <file2.tab> ...", call.=FALSE)
}

# Train novel models using data with matched normal
data <- read.table(files[1],header=T)
if(length(files) > 1){
  for (i in 2:length(files)){
    newdata <- read.table(files[i],header=T)
    data <- rbind(data,newdata)
  }
  rm(newdata)
}


# clean up data
snp_mut_type <- c('nonsynonymous_SNV','splicing',
                  'stopgain','stoploss');
data <- data[data$N_RefDP + data$N_VarDP >= 10,]
data <- data[data$MutType %in% snp_mut_type,]
data <- data[(data$VAF >= 0.05 & data$VarDP >=10) | data$VAF >= 0.2,]

# remove common SNPs with freq>0.01
data <- data[data$snp138 == 0 | data$PopFreqMax < 0.01,]
# remove data without normal information
data <- data[!is.na(data$N_VAF),]
if(nrow(data)==0){
  stop("Lack of Normal data. Quit.", call.=FALSE)
}

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
  }
}  

for (i in 1:nrow(data)){
  if(is.na(data$T2_RefDP[i] + data$T2_VarDP[i])){
    data$Specific[i] <- NA
  }else{
    if(data$T2_VAF[i] < 0.02 || data$T2_VarDP <= 1){
      data$Specific[i] <- 1
      }else{
        data$Specific[i] <- 0
      }
  }
  if(data$N_VAF[i] < 0.02){
    data$Somatic[i] <- 1
  }else{
    data$Somatic[i] <- 0
  }
}
data$Specific <- factor(data$Specific)
data$Somatic <- factor(data$Somatic)
rm(i,score,med_score,nlof,snp_mut_type,files)

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

feat_and_class <- c(
  'PopFreqMax','snp138','SIFT','Polyphen2_HDIV',
  'MutationAssessor','PROVEAN','VEST3','CADD_phred',
  'DANN','fathmm_MKL_coding','MetaSVM','MetaLR','GERP',
  'phyloP20way_mammalian', 'phastCons7way_vertebrate',
  'phastCons20way_mammalian','SiPhy_29way_logOdds', 
  'VarDP','VAF','BAF_s40M','LOFMut','logBinomP','aVAF',
  'Somatic'
)
mldata <- data[,feat_and_class]
## create model using random forest
x.rf <- randomForest(Somatic ~ ., data=mldata)
save(x.rf, file = "self_trained_RDA/RF.single.rda")

# remove data without T2 information
data <- data[!is.na(data$T2_VAF),]
if (nrow(data) >= 10){                                 # with T2
  feat_and_class <- c(
    'PopFreqMax','snp138','SIFT','Polyphen2_HDIV',
    'MutationAssessor','PROVEAN','VEST3','CADD_phred',
    'DANN','fathmm_MKL_coding','MetaSVM','MetaLR','GERP', 
    'phyloP20way_mammalian', 'phastCons7way_vertebrate',
    'phastCons20way_mammalian','SiPhy_29way_logOdds', 
    'VarDP','VAF','BAF_s40M','T2_VarDP','T2_VAF',
    'T2_BAF_s40M','LOFMut','logBinomP','T2_logBinomP',
    'T2_VAFdiff','aVAF','T2_aVAF','Somatic'
  )
  mldata <- data[,feat_and_class]
  ## create model using random forest
  x.rf <- randomForest(Somatic ~ ., data=mldata)
  save(x.rf, file = "self_trained_RDA/RF.multiple.rda")
}




