#03/01/2016 
# script to perform analyses for metabolic syndrome

rm(list=ls())
library(buRlo)

############################## STEP 1: load phenotype data #####################
#load all the phenotypes from the merged genotype set
carl_all_gen <- load.burlo.data(pop = "CARL", phen = "general", data.dir = "/nfs/servizio/") 

#inclusion criteria: call rate > 95%
carl_all_summary <- perid.summary(carl_all_gen)
# carl_all_summary_call_95 <- carl_all_summary[which(carl_all_summary$CallPP > 0.95),]
carl_all_summary_call_95 <- row.names(carl_all_summary[which(carl_all_summary$CallPP > 0.95),])


#select only relevant columns to calculate the score
#conversion from our pheno names to the GENERAL used names
# id=carl_all_gen@phdata$id
# sex=carl_all_gen@phdata$sex
# village=carl_all_gen@phdata$village
# age=carl_all_gen@phdata$age
# BMI=carl_all_gen@phdata$BMI
# WC=carl_all_gen@phdata$VITA
# SBP=carl_all_gen@phdata$PSMEAN
# DBP=carl_all_gen@phdata$PDMEAN
# GLU=carl_all_gen@phdata$GLUCOSIO
# HDL=carl_all_gen@phdata$COL_HDL
# TRIG=carl_all_gen@phdata$TRIGLICERIDI

carl_all_metabolic <- data.frame(id=carl_all_gen@phdata$id,sex=carl_all_gen@phdata$sex,village=carl_all_gen@phdata$village,age=carl_all_gen@phdata$age,WC=carl_all_gen@phdata$VITA,BMI=carl_all_gen@phdata$BMI,SBP=carl_all_gen@phdata$PSMEAN,DBP=carl_all_gen@phdata$PDMEAN,GLU=carl_all_gen@phdata$GLUCOSIO,HDL=carl_all_gen@phdata$COL_HDL,TRIG=carl_all_gen@phdata$TRIGLICERIDI)

#inclusion criteria
carl_all_metabolic <- carl_all_metabolic[which(carl_all_metabolic$id %in% carl_all_summary_call_95),]
carl_all_metabolic <- carl_all_metabolic[which(carl_all_metabolic$age >= 18),]


carl_all_metabolic_males <- carl_all_metabolic[which(carl_all_metabolic$sex == 1),]
carl_all_metabolic_females <- carl_all_metabolic[which(carl_all_metabolic$sex == 0),]

dim(carl_all_metabolic_males)
dim(carl_all_metabolic_females)

#remove all NA data
carl_all_metabolic_males <- carl_all_metabolic_males[which(!is.na(carl_all_metabolic_males$BMI)),]
carl_all_metabolic_males <- carl_all_metabolic_males[which(!is.na(carl_all_metabolic_males$WC)),]
carl_all_metabolic_males <- carl_all_metabolic_males[which(!is.na(carl_all_metabolic_males$SBP)),]
carl_all_metabolic_males <- carl_all_metabolic_males[which(!is.na(carl_all_metabolic_males$DBP)),]
carl_all_metabolic_males <- carl_all_metabolic_males[which(!is.na(carl_all_metabolic_males$GLU)),]
carl_all_metabolic_males <- carl_all_metabolic_males[which(!is.na(carl_all_metabolic_males$HDL)),]
carl_all_metabolic_males <- carl_all_metabolic_males[which(!is.na(carl_all_metabolic_males$TRIG)),]

carl_all_metabolic_females <- carl_all_metabolic_females[which(!is.na(carl_all_metabolic_females$BMI)),]
carl_all_metabolic_females <- carl_all_metabolic_females[which(!is.na(carl_all_metabolic_females$WC)),]
carl_all_metabolic_females <- carl_all_metabolic_females[which(!is.na(carl_all_metabolic_females$SBP)),]
carl_all_metabolic_females <- carl_all_metabolic_females[which(!is.na(carl_all_metabolic_females$DBP)),]
carl_all_metabolic_females <- carl_all_metabolic_females[which(!is.na(carl_all_metabolic_females$GLU)),]
carl_all_metabolic_females <- carl_all_metabolic_females[which(!is.na(carl_all_metabolic_females$HDL)),]
carl_all_metabolic_females <- carl_all_metabolic_females[which(!is.na(carl_all_metabolic_females$TRIG)),]

#Recode sex as 1 for Males and 2 for Females
carl_all_metabolic_males$sex <- 1
carl_all_metabolic_females$sex <- 2

dim(carl_all_metabolic_males)
dim(carl_all_metabolic_females)

#calculate the MetS_score (gender specific)
carl_all_metabolic_males$G_males <- 0.645*carl_all_metabolic_males$WC + 0.933*carl_all_metabolic_males$BMI + 0.059*carl_all_metabolic_males$SBP + 0.087*carl_all_metabolic_males$DBP + 0.011*carl_all_metabolic_males$GLU - 0.022*carl_all_metabolic_males$HDL+ 0.003*carl_all_metabolic_males$TRIG - 63.0
carl_all_metabolic_females$G_females <- 0.342*carl_all_metabolic_females$WC+ 0.636*carl_all_metabolic_females$BMI + 0.133*carl_all_metabolic_females$SBP+ 0.146*carl_all_metabolic_females$DBP + 0.021*carl_all_metabolic_females$GLU - 0.027*carl_all_metabolic_females$HDL+ 0.009*carl_all_metabolic_females$TRIG - 44.4

head(carl_all_metabolic_males,20)
head(carl_all_metabolic_females,20)

#write tables
write.table(carl_all_metabolic_males,file=paste("carl_all_metabolic_males_MetS_score.csv",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(carl_all_metabolic_females,file=paste("carl_all_metabolic_females_MetS_score.csv",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
