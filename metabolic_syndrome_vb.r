#03/01/2016 
# script to perform analyses for metabolic syndrome

rm(list=ls())
library(buRlo)

############################## STEP 1: load phenotype data #####################
#load all the phenotypes from the merged genotype set
#/home/cocca/analyses/MetabolicSyndrome/VBI/tratti_VB.csv

# vb_all_gen <- load.burlo.data(pop = "vb", phen = "general", data.dir = "/nfs/servizio/") 
vb_all_gen <- read.table("/home/cocca/analyses/MetabolicSyndrome/VBI/tratti_VB_2.csv",sep="\t",header=T)

#inclusion criteria: call rate > 95%
# vb_all_summary <- perid.summary(vb_all_gen)
# vb_all_summary_call_95 <- vb_all_summary[which(vb_all_summary$CallPP > 0.95),]
# vb_all_summary_call_95 <- row.names(vb_all_summary[which(vb_all_summary$CallPP > 0.95),])


#select only relevant columns to calculate the score
#conversion from our pheno names to the GENERAL used names
# id=vb_all_gen$clinic_id
# sex=vb_all_gen$sex
# age=vb_all_gen$age..years.
# WC=vb_all_gen$waist..cm.
# BMI=vb_all_gen$bmi..kg.m.2.
# SBP=vb_all_gen$average.pao.sist..mmHg.
# DBP=vb_all_gen$average.pao.diast..mmHg.
# GLU=vb_all_gen$glucose..mg.dl.
# HDL=vb_all_gen$HDL..mg.dl.
# TRIG=vb_all_gen$tryglicerids..mg.dl.
# LLM=vb_all_gen$drug.ipolipemizzanti
# AHM=vb_all_gen$drug.antipertensivi
# ADM=vb_all_gen$drug.antidiabetici
# DIABETES=vb_all_gen$self.reported.diabetes

vb_all_metabolic <- data.frame(id=vb_all_gen$clinic_id,sex=vb_all_gen$sex,age=vb_all_gen$age..years.,WC=vb_all_gen$waist..cm.,BMI=vb_all_gen$bmi..kg.m.2.,SBP=vb_all_gen$average.pao.sist..mmHg.,DBP=vb_all_gen$average.pao.diast..mmHg.,GLU=vb_all_gen$glucose..mg.dl.,HDL=vb_all_gen$HDL..mg.dl.,TRIG=vb_all_gen$tryglicerids..mg.dl.,LLM=vb_all_gen$drug.ipolipemizzanti,AHM=vb_all_gen$drug.antipertensivi,ADM=vb_all_gen$drug.antidiabetici,DIABETES=vb_all_gen$self.reported.diabetes)

#inclusion criteria
# vb_all_metabolic <- vb_all_metabolic[which(vb_all_metabolic$id %in% vb_all_summary_call_95),]
vb_all_metabolic <- vb_all_metabolic[which(vb_all_metabolic$age >= 18),]

#read genotyped samples 
vb_all_genotyped <- read.table("/home/cocca/analyses/MetabolicSyndrome/VBI/VBI_1785_genotyped.sample",header=T)

#keep only genotyped samples
vb_all_metabolic <- vb_all_metabolic[which(vb_all_metabolic$id %in% vb_all_genotyped$ID),]

#split by gender
vb_all_metabolic_males <- vb_all_metabolic[which(vb_all_metabolic$sex == "M"),]
vb_all_metabolic_females <- vb_all_metabolic[which(vb_all_metabolic$sex == "F"),]

dim(vb_all_metabolic_males)
dim(vb_all_metabolic_females)

summary(vb_all_metabolic_males)
summary(vb_all_metabolic_females)


#remove all NA data
vb_all_metabolic_males <- vb_all_metabolic_males[which(!is.na(vb_all_metabolic_males$BMI)),]
vb_all_metabolic_males <- vb_all_metabolic_males[which(!is.na(vb_all_metabolic_males$WC)),]
vb_all_metabolic_males <- vb_all_metabolic_males[which(!is.na(vb_all_metabolic_males$SBP)),]
vb_all_metabolic_males <- vb_all_metabolic_males[which(!is.na(vb_all_metabolic_males$DBP)),]
vb_all_metabolic_males <- vb_all_metabolic_males[which(!is.na(vb_all_metabolic_males$GLU)),]
vb_all_metabolic_males <- vb_all_metabolic_males[which(!is.na(vb_all_metabolic_males$HDL)),]
vb_all_metabolic_males <- vb_all_metabolic_males[which(!is.na(vb_all_metabolic_males$TRIG)),]

vb_all_metabolic_females <- vb_all_metabolic_females[which(!is.na(vb_all_metabolic_females$BMI)),]
vb_all_metabolic_females <- vb_all_metabolic_females[which(!is.na(vb_all_metabolic_females$WC)),]
vb_all_metabolic_females <- vb_all_metabolic_females[which(!is.na(vb_all_metabolic_females$SBP)),]
vb_all_metabolic_females <- vb_all_metabolic_females[which(!is.na(vb_all_metabolic_females$DBP)),]
vb_all_metabolic_females <- vb_all_metabolic_females[which(!is.na(vb_all_metabolic_females$GLU)),]
vb_all_metabolic_females <- vb_all_metabolic_females[which(!is.na(vb_all_metabolic_females$HDL)),]
vb_all_metabolic_females <- vb_all_metabolic_females[which(!is.na(vb_all_metabolic_females$TRIG)),]

#Recode sex as 1 for Males and 2 for Females
vb_all_metabolic_males$sex <- 1
vb_all_metabolic_females$sex <- 2

dim(vb_all_metabolic_males)
dim(vb_all_metabolic_females)

#calculate the MetS_score (gender specific)
vb_all_metabolic_males$MetS_score <- 0.645*vb_all_metabolic_males$WC + 0.933*vb_all_metabolic_males$BMI + 0.059*vb_all_metabolic_males$SBP + 0.087*vb_all_metabolic_males$DBP + 0.011*vb_all_metabolic_males$GLU - 0.022*vb_all_metabolic_males$HDL+ 0.003*vb_all_metabolic_males$TRIG - 63.0
vb_all_metabolic_females$MetS_score <- 0.342*vb_all_metabolic_females$WC+ 0.636*vb_all_metabolic_females$BMI + 0.133*vb_all_metabolic_females$SBP+ 0.146*vb_all_metabolic_females$DBP + 0.021*vb_all_metabolic_females$GLU - 0.027*vb_all_metabolic_females$HDL+ 0.009*vb_all_metabolic_females$TRIG - 44.4

head(vb_all_metabolic_males,20)
head(vb_all_metabolic_females,20)

#write tables
write.table(vb_all_metabolic_males,file=paste("vb_all_metabolic_males_MetS_score.csv",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(vb_all_metabolic_females,file=paste("vb_all_metabolic_females_MetS_score.csv",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

# MetS prevalence based on IDF definition calculated by sex:
# Waist circumference (WC) ≥ 94 cm in men, plus two of the following:
# triglycerides ≥150 mg/dl or lipid-lowering medication,
# HDL-cholesterol <40 mg/dl in men or lipid-lowering medication,
# systolic blood pressure (SBP) ≥130 mmHg or diastolic blood pressure (DBP) ≥85 mmHg, or antihypertensive medication and,
# fasting blood glucose ≥100 mg/dL or previously diagnosed type 2 diabetes.
vb_all_metabolic_males$affected <- 0
for(i in 1:length(vb_all_metabolic_males$id)){
	if (vb_all_metabolic_males[i,]$WC >= 94){
		affected_status <- 0
		if (vb_all_metabolic_males[i,]$TRIG >= 150 | vb_all_metabolic_males[i,]$LLM == 1 ) {
			affected_status <- affected_status +1
		}
		if(vb_all_metabolic_males[i,]$HDL < 40 | vb_all_metabolic_males[i,]$LLM == 1) {
			affected_status <- affected_status +1
		}
		if(vb_all_metabolic_males[i,]$SBP >= 130 | vb_all_metabolic_males[i,]$DBP >= 85 | vb_all_metabolic_males[i,]$AHM == 1) {
			affected_status <- affected_status +1
		}
		if(vb_all_metabolic_males[i,]$GLU >= 100 | vb_all_metabolic_males[i,]$DIABETES == 1) {
			affected_status <- affected_status +1
		}
		if (affected_status >= 2){
		vb_all_metabolic_males[i,]$affected <- 1
		}
	} else {
		vb_all_metabolic_males[i,]$affected <- 0
	}
}

# Waist circumference (WC) ≥ 80 cm in women, plus two of the following:
# triglycerides ≥150 mg/dl or lipid-lowering medication,
# HDL-cholesterol <50 mg/dl in women or lipid-lowering medication,
# systolic blood pressure (SBP) ≥130 mmHg or diastolic blood pressure (DBP) ≥85 mmHg, or antihypertensive medication and,
# fasting blood glucose ≥100 mg/dL or previously diagnosed type 2 diabetes.
vb_all_metabolic_females$affected <- 0
for(i in 1:length(vb_all_metabolic_females$id)){
	if (vb_all_metabolic_females[i,]$WC >= 80){
		affected_status <- 0
		if (vb_all_metabolic_females[i,]$TRIG >= 150 | vb_all_metabolic_females[i,]$LLM == 1 ) {
			affected_status <- affected_status +1
		}
		if(vb_all_metabolic_females[i,]$HDL < 50 | vb_all_metabolic_females[i,]$LLM == 1) {
			affected_status <- affected_status +1
		}
		if(vb_all_metabolic_females[i,]$SBP >= 130 | vb_all_metabolic_females[i,]$DBP >= 85 | vb_all_metabolic_females[i,]$AHM == 1) {
			affected_status <- affected_status +1
		}
		if(vb_all_metabolic_females[i,]$GLU >= 100 | vb_all_metabolic_females[i,]$DIABETES == 1) {
			affected_status <- affected_status +1
		}
		if (affected_status >= 2){
		vb_all_metabolic_females[i,]$affected <- 1
		}
	} else {
		vb_all_metabolic_females[i,]$affected <- 0
	}
}

dim(vb_all_metabolic_males[vb_all_metabolic_males$affected == 1,])
dim(vb_all_metabolic_females[vb_all_metabolic_females$affected == 1,])
# > dim(vb_all_metabolic_males[vb_all_metabolic_males$affected == 1,])
# [1] 176  16
# > dim(vb_all_metabolic_females[vb_all_metabolic_females$affected == 1,])
# [1] 269  16


#put all back together to compute stats
vb_all_metabolic_cleaned <- rbind(vb_all_metabolic_males,vb_all_metabolic_females)

#descriptive stats
summary(vb_all_metabolic_cleaned)
sd(vb_all_metabolic_cleaned$age)
sd(vb_all_metabolic_cleaned$MetS_score)
library(pastecs)
stat.desc(vb_all_metabolic_cleaned)
library(moments)
skewness(vb_all_metabolic_cleaned$MetS_score)
kurtosis(vb_all_metabolic_cleaned$MetS_score)

#need to plot the istogram for MetS_score
# for males, females and all together
library(ggplot2)
all_datasets <- c("vb_all_metabolic_cleaned","vb_all_metabolic_males","vb_all_metabolic_females")
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

# > skewness(vb_all_metabolic_cleaned$MetS_score)
# [1] 0.419236
# > kurtosis(vb_all_metabolic_cleaned$MetS_score)
# [1] 3.410815
#MetS_score>36.33 on overal sample
dim(vb_all_metabolic_cleaned[vb_all_metabolic_cleaned$MetS_score > 36.33,])
vb_all_metabolic_cleaned$MetS_affected <- ifelse(vb_all_metabolic_cleaned$MetS_score > 36.33, 1, 0)
# 644

dim(vb_all_metabolic_cleaned)
# 1741

#Calculate MetS prevalence as defined in the summary sheet:
# Prevalence= (number of diseased subjects)/(total number of subjects)
dim(vb_all_metabolic_cleaned[vb_all_metabolic_cleaned$affected == 1,])
# [1] 445  16

# Cohen's K:
# Cohen's kappa coefficient is a statistic which measures inter-rater agreement for qualitative (categorical) items.
##Cohen's K computation using an R package
# install.packages("irr")
library(irr)

#After creating a data frame containing only the 2 columns with the 2 diagnosis (ratings):
ratings <- (vb_all_metabolic_cleaned[,c("MetS_affected","affected")])
kappa2(ratings)
 




