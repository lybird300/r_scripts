#03/01/2016 
# script to perform analyses for metabolic syndrome

rm(list=ls())
library(buRlo)

############################## STEP 1: load phenotype data #####################
#load all the phenotypes from the merged genotype set
fvg_all_gen <- load.burlo.data(pop = "FVG", phen = "general", data.dir = "/nfs/servizio/") 

#inclusion criteria: call rate > 95%
fvg_all_summary <- perid.summary(fvg_all_gen)
# fvg_all_summary_call_95 <- fvg_all_summary[which(fvg_all_summary$CallPP > 0.95),]
fvg_all_summary_call_95 <- row.names(fvg_all_summary[which(fvg_all_summary$CallPP > 0.95),])

#inclusion criteria
# fvg_all_metabolic <- fvg_all_gen[which(fvg_all_gen@phdata$id %in% fvg_all_summary_call_95),]
# fvg_all_metabolic <- fvg_all_metabolic[which(fvg_all_metabolic@phdata$age >= 18),]

# selected_cols <- c("id","sex","village","age","BMI","VITA","PSMEAN","PDMEAN","GLUCOSIO","COL_HDL","TRIGLICERIDI","F_lip","F_iper","F_diab","diabete")
#select only relevant columns to calculate the score
#conversion from our pheno names to the GENERAL used names
# id=fvg_all_gen@phdata$id
# sex=fvg_all_gen@phdata$sex
# village=fvg_all_gen@phdata$village
# age=fvg_all_gen@phdata$age
# BMI=fvg_all_gen@phdata$BMI
# WC=fvg_all_gen@phdata$VITA
# SBP=fvg_all_gen@phdata$PSMEAN
# DBP=fvg_all_gen@phdata$PDMEAN
# GLU=fvg_all_gen@phdata$GLUCOSIO
# HDL=fvg_all_gen@phdata$COL_HDL
# TRIG=fvg_all_gen@phdata$TRIGLICERIDI
# LLM=fvg_all_gen@phdata$F_lip
# AHM=fvg_all_gen@phdata$F_iper
# ADM=fvg_all_gen@phdata$F_diab
# DIABETES=fvg_all_gen@phdata$diabete

fvg_all_metabolic <- data.frame(id=fvg_all_gen@phdata$id,sex=fvg_all_gen@phdata$sex,village=fvg_all_gen@phdata$village,age=fvg_all_gen@phdata$age,WC=fvg_all_gen@phdata$VITA,BMI=fvg_all_gen@phdata$BMI,SBP=fvg_all_gen@phdata$PSMEAN,DBP=fvg_all_gen@phdata$PDMEAN,GLU=fvg_all_gen@phdata$GLUCOSIO,HDL=fvg_all_gen@phdata$COL_HDL,TRIG=fvg_all_gen@phdata$TRIGLICERIDI,LLM=fvg_all_gen@phdata$F_lip,AHM=fvg_all_gen@phdata$F_iper,ADM=fvg_all_gen@phdata$F_diab,DIABETES=fvg_all_gen@phdata$diabete)

#add data retrieved from the new anamnesis
fvg_all_metabolic_missing <- read.table("/home/cocca/analyses/MetabolicSyndrome/FVG/fvg_all_metabolic_missing.txt",header=T)
#remove useless column
fvg_all_metabolic_missing$MetS_score <- NULL
#fix sex
fvg_all_metabolic_missing_males <- fvg_all_metabolic_missing[which(fvg_all_metabolic_missing$sex == 1),]
fvg_all_metabolic_missing_females <- fvg_all_metabolic_missing[which(fvg_all_metabolic_missing$sex == 2),]
fvg_all_metabolic_missing_males$sex <- 1
fvg_all_metabolic_missing_females$sex <- 0

fvg_all_metabolic_missing_fix <- rbind(fvg_all_metabolic_missing_males,fvg_all_metabolic_missing_females)

#inclusion criteria
fvg_all_metabolic <- fvg_all_metabolic[which(fvg_all_metabolic$id %in% fvg_all_summary_call_95),]
fvg_all_metabolic <- fvg_all_metabolic[which(fvg_all_metabolic$age >= 18),]
fvg_all_metabolic <- fvg_all_metabolic[-which(fvg_all_metabolic$id %in% fvg_all_metabolic_missing_fix$id),]
fvg_all_metabolic <- rbind(fvg_all_metabolic,fvg_all_metabolic_missing_fix)

fvg_all_metabolic_males <- fvg_all_metabolic[which(fvg_all_metabolic$sex == 1),]
fvg_all_metabolic_females <- fvg_all_metabolic[which(fvg_all_metabolic$sex == 0),]

dim(fvg_all_metabolic_males)
# [1] 610  15
dim(fvg_all_metabolic_females)
# [1] 834  15

#set to 0 missing parameters for the IDF calculation
fvg_all_metabolic_males[is.na(fvg_all_metabolic_males$LLM),]$LLM <- 0
fvg_all_metabolic_males[is.na(fvg_all_metabolic_males$AHM),]$AHM <- 0
fvg_all_metabolic_males[is.na(fvg_all_metabolic_males$ADM),]$ADM <- 0
fvg_all_metabolic_males[is.na(fvg_all_metabolic_males$DIABETES),]$DIABETES <- 0

fvg_all_metabolic_females[is.na(fvg_all_metabolic_females$LLM),]$LLM <- 0
fvg_all_metabolic_females[is.na(fvg_all_metabolic_females$AHM),]$AHM <- 0
fvg_all_metabolic_females[is.na(fvg_all_metabolic_females$ADM),]$ADM <- 0
fvg_all_metabolic_females[is.na(fvg_all_metabolic_females$DIABETES),]$DIABETES <- 0

#Recode sex as 1 for Males and 2 for Females
# fvg_all_metabolic_males$sex <- 1
# fvg_all_metabolic_females$sex <- 2

dim(fvg_all_metabolic_males)
dim(fvg_all_metabolic_females)

#calculate the MetS_score (gender specific)
fvg_all_metabolic_males$MetS_score <- (0.645*fvg_all_metabolic_males$WC + 0.933*fvg_all_metabolic_males$BMI + 0.059*fvg_all_metabolic_males$SBP + 0.087*fvg_all_metabolic_males$DBP + 0.011*fvg_all_metabolic_males$GLU - 0.022*fvg_all_metabolic_males$HDL+ 0.003*fvg_all_metabolic_males$TRIG - 63.0)
fvg_all_metabolic_females$MetS_score <- (0.342*fvg_all_metabolic_females$WC+ 0.636*fvg_all_metabolic_females$BMI + 0.133*fvg_all_metabolic_females$SBP+ 0.146*fvg_all_metabolic_females$DBP + 0.021*fvg_all_metabolic_females$GLU - 0.027*fvg_all_metabolic_females$HDL+ 0.009*fvg_all_metabolic_females$TRIG - 44.4)

head(fvg_all_metabolic_males,20)
head(fvg_all_metabolic_females,20)

#write tables
write.table(fvg_all_metabolic_males,file=paste("fvg_all_metabolic_males_MetS_score.csv",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(fvg_all_metabolic_females,file=paste("fvg_all_metabolic_females_MetS_score.csv",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#put all back together to calculate stats and to have the dataset for the analyses
fvg_all_metabolic_all_score <- rbind(fvg_all_metabolic_males,fvg_all_metabolic_females)

# MetS prevalence based on IDF definition calculated by sex:
# Waist circumference (WC) ≥ 94 cm in men, plus two of the following:
# triglycerides ≥150 mg/dl or lipid-lowering medication,
# HDL-cholesterol <40 mg/dl in men or lipid-lowering medication,
# systolic blood pressure (SBP) ≥130 mmHg or diastolic blood pressure (DBP) ≥85 mmHg, or antihypertensive medication and,
# fasting blood glucose ≥100 mg/dL or previously diagnosed type 2 diabetes.

#remove all NA data
fvg_all_metabolic_males <- fvg_all_metabolic_males[which(!is.na(fvg_all_metabolic_males$BMI)),]
fvg_all_metabolic_males <- fvg_all_metabolic_males[which(!is.na(fvg_all_metabolic_males$WC)),]
fvg_all_metabolic_males <- fvg_all_metabolic_males[which(!is.na(fvg_all_metabolic_males$SBP)),]
fvg_all_metabolic_males <- fvg_all_metabolic_males[which(!is.na(fvg_all_metabolic_males$DBP)),]
fvg_all_metabolic_males <- fvg_all_metabolic_males[which(!is.na(fvg_all_metabolic_males$GLU)),]
fvg_all_metabolic_males <- fvg_all_metabolic_males[which(!is.na(fvg_all_metabolic_males$HDL)),]
fvg_all_metabolic_males <- fvg_all_metabolic_males[which(!is.na(fvg_all_metabolic_males$TRIG)),]

fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$BMI)),]
fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$WC)),]
fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$SBP)),]
fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$DBP)),]
fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$GLU)),]
fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$HDL)),]
fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$TRIG)),]
fvg_all_metabolic_females <- fvg_all_metabolic_females[which(!is.na(fvg_all_metabolic_females$TRIG)),]

dim(fvg_all_metabolic_males)
# [1] 363  16
dim(fvg_all_metabolic_females)
# [1] 497  16

fvg_all_metabolic_males$affected <- 0
# fvg_all_metabolic_males$affected_status <- 0
for(i in 1:length(fvg_all_metabolic_males$id)){
  affected_status <- 0
  if (fvg_all_metabolic_males[i,]$WC >= 94){
    if (fvg_all_metabolic_males[i,]$TRIG >= 150 | (!is.na(fvg_all_metabolic_males[i,]$LLM) & fvg_all_metabolic_males[i,]$LLM == 1) ) {
      affected_status <- affected_status +1
    }
    if(fvg_all_metabolic_males[i,]$HDL < 40 | (!is.na(fvg_all_metabolic_males[i,]$LLM) & fvg_all_metabolic_males[i,]$LLM == 1)) {
      affected_status <- affected_status +1
    }
    if(fvg_all_metabolic_males[i,]$SBP >= 130 | fvg_all_metabolic_males[i,]$DBP >= 85 | (!is.na(fvg_all_metabolic_males[i,]$AHM) & fvg_all_metabolic_males[i,]$AHM == 1)) {
      affected_status <- affected_status +1
    }
    if(fvg_all_metabolic_males[i,]$GLU >= 100 | (!is.na(fvg_all_metabolic_males[i,]$DIABETES) & fvg_all_metabolic_males[i,]$DIABETES == 1)) {
      affected_status <- affected_status +1
    }
    if (affected_status >= 2){
    fvg_all_metabolic_males[i,]$affected <- 1
    }
  } else {
    fvg_all_metabolic_males[i,]$affected <- 0
  }
  # fvg_all_metabolic_males[i,]$affected_status <- affected_status

}

# Waist circumference (WC) ≥ 80 cm in women, plus two of the following:
# triglycerides ≥150 mg/dl or lipid-lowering medication,
# HDL-cholesterol <50 mg/dl in women or lipid-lowering medication,
# systolic blood pressure (SBP) ≥130 mmHg or diastolic blood pressure (DBP) ≥85 mmHg, or antihypertensive medication and,
# fasting blood glucose ≥100 mg/dL or previously diagnosed type 2 diabetes.
fvg_all_metabolic_females$affected <- 0
for(i in 1:length(fvg_all_metabolic_females$id)){
  if (fvg_all_metabolic_females[i,]$WC >= 80){
    affected_status <- 0
    if (fvg_all_metabolic_females[i,]$TRIG >= 150 | (!is.na(fvg_all_metabolic_females[i,]$LLM) & fvg_all_metabolic_females[i,]$LLM == 1)) {
      affected_status <- affected_status +1
    }
    if(fvg_all_metabolic_females[i,]$HDL < 50 | (!is.na(fvg_all_metabolic_females[i,]$LLM) & fvg_all_metabolic_females[i,]$LLM == 1)) {
      affected_status <- affected_status +1
    }
    if(fvg_all_metabolic_females[i,]$SBP >= 130 | fvg_all_metabolic_females[i,]$DBP >= 85 | (!is.na(fvg_all_metabolic_females[i,]$AHM) & fvg_all_metabolic_females[i,]$AHM == 1)) {
      affected_status <- affected_status +1
    }
    if(fvg_all_metabolic_females[i,]$GLU >= 100 | (!is.na(fvg_all_metabolic_females[i,]$DIABETES) & fvg_all_metabolic_females[i,]$DIABETES == 1)) {
      affected_status <- affected_status +1
    }
    if (affected_status >= 2){
    fvg_all_metabolic_females[i,]$affected <- 1
    }
  } else {
    fvg_all_metabolic_females[i,]$affected <- 0
  }
}

dim(fvg_all_metabolic_males[fvg_all_metabolic_males$affected == 1,])
dim(fvg_all_metabolic_females[fvg_all_metabolic_females$affected == 1,])
# > dim(fvg_all_metabolic_males[fvg_all_metabolic_males$affected == 1,])
# [1] 126  17
# > dim(fvg_all_metabolic_females[fvg_all_metabolic_females$affected == 1,])
# [1] 135  17

#put all back together to compute stats
fvg_all_metabolic_cleaned <- rbind(fvg_all_metabolic_males,fvg_all_metabolic_females)
# > dim(fvg_all_metabolic_cleaned)
# [1] 859  17

#descriptive stats
library(pastecs)
stat.desc(fvg_all_metabolic)
stat.desc(fvg_all_metabolic_all_score)

sd(fvg_all_metabolic_cleaned$age)
sd(fvg_all_metabolic_cleaned$MetS_score)

stat.desc(fvg_all_metabolic_cleaned)
library(moments)
skewness(fvg_all_metabolic_cleaned$MetS_score)
kurtosis(fvg_all_metabolic_cleaned$MetS_score)
# > skewness(fvg_all_metabolic_cleaned$MetS_score)
# [1] 0.6240077
# > kurtosis(fvg_all_metabolic_cleaned$MetS_score)
# [1] 4.016046

#need to plot the istogram for MetS_score
# for males, females and all together
library(ggplot2)
all_datasets <- c("fvg_all_metabolic_cleaned","fvg_all_metabolic_males","fvg_all_metabolic_females")
for (dataset in all_datasets){
  current_dataset <- get(dataset)
  pl <- ggplot(current_dataset)
  pl <- pl + aes(MetS_score)
  pl <- pl + geom_histogram(binwidth=0.5)
  pl <- pl + labs(title = paste(dataset,"MetS score",sep=" "))
  pl <- pl + theme_bw()
  pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
  pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
  pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
  # pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
  # ggsave(filename=paste(getwd(),"/complete_19102015_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
  ggsave(filename=paste(getwd(),"/",dataset,".jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
  
}

#MetS_score>36.33 on overal sample
dim(fvg_all_metabolic_cleaned[fvg_all_metabolic_cleaned$MetS_score > 36.33,])
fvg_all_metabolic_cleaned$MetS_affected <- ifelse(fvg_all_metabolic_cleaned$MetS_score > 36.33, 1, 0)
# 349
dim(fvg_all_metabolic_cleaned)
# 859  17

#Calculate MetS prevalence as defined in the summary sheet:
# Prevalence= (number of diseased subjects)/(total number of subjects)
dim(fvg_all_metabolic_cleaned[fvg_all_metabolic_cleaned$affected == 1,])
# [1] 261  17

# Cohen's K:
# Cohen's kappa coefficient is a statistic which measures inter-rater agreement for qualitative (categorical) items.
##Cohen's K computation using an R package
# install.packages("irr")
library(irr)

#After creating a data frame containing only the 2 columns with the 2 diagnosis (ratings):
ratings <- (fvg_all_metabolic_cleaned[,c("MetS_affected","affected")])
kappa2(ratings)
# Cohen's Kappa for 2 Raters (Weights: unweighted)

#  Subjects = 859 
#    Raters = 2 
#     Kappa = 0.558 

#         z = 16.8 
#   p-value = 0 

#########################genotype data preparation for genabel
phenofile <- "/home/cocca/analyses/MetabolicSyndrome/FVG/fvg_all_metabolic_ALL_MetS_score.csv"
phenotypes <- read.table(phenofile,header=T,sep="\t")

all_gen <- load.burlo.data(pop = "FVG", phen = "general", data.dir = "/nfs/servizio/") 

all_to_analyze <- all_gen[which(all_gen@phdata$id %in% phenotypes$id),]

export.plink(all_to_analyze,filebasename="FVG", phenotypes=c("sex"), transpose=TRUE)
convert.snp.tped(tpedfile="/home/cocca/analyses/MetabolicSyndrome/FVG/FVG.tped",tfamfile="/home/cocca/analyses/MetabolicSyndrome/FVG/FVG.tfam",outfile="FVG_out")


