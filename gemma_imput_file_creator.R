#Script to create a correct gemma imput formatted file
rm(list=ls())

#set arguments
# trait <- "HDL"
# trait <- "LDL"
# trait <- "TC"
# trait <- "TG"
# trait <- "TotalFM"
# trait <- "TotalLM"
# trait <- "DiastolicBP"
# trait <- "hemoglobin"
# trait <- "SystolicBP"
# trait <- "mch"
# trait <- "mchc"
# trait <- "mcv"
# trait <- "packedcellvolume"
# trait <- "plateletcount"
# trait <- "rbc"
# trait <- "whitebloodcell"
# trait <- "CREA"
# trait <- "FOSFAT"
# trait <- "SODIUM"
# trait <- "UREA"
# trait <- "URICACID"
# trait <- "HIPadjBMI"
# trait <- "WAISTadjBMI"
# trait <- "WHRadjBMI"
# trait <- "ALK_F"
# trait <- "BIL"
# trait <- "GGT"
# trait <- "DBP"
# trait <- "SBP"
trait <- "HR"


# basepath <- "/nfs/users/nfs_m/mc14/Work/SANGER/FVG/PHENO/"
basepath <- "/nfs/users/nfs_m/mc14/Work/SANGER/CARLANTINO/PHENO/"
phenopath <- paste(basepath,"CARDIO/",sep="")
pop <- "CARL"
# pop <- "FVG"

#if there is gender separation
male_table <- paste(phenopath,"res_",trait,"_",pop,"_males.txt",sep="")
female_table <- paste(phenopath,"res_",trait,"_",pop,"_females.txt",sep="")
#for complete table
#complete_table <- paste(phenopath,"res_",trait,"_",pop,"_all.txt",sep="")

sorted_file <- paste(basepath,"imputation_order_id.samples",sep="")
out_file <- paste(phenopath,"res_",trait,"_",pop,"_all.phen",sep="")

#read males and females files
males <- read.table(male_table,header=F,sep="\t",skip=1)
females <- read.table(female_table,header=F,sep="\t",skip=1)
#for complete table
complete <- read.table(complete_table,header=F,sep="\t",skip=1)

#check data
str(males)
str(females)

#format colnames
colnames(males) <- c("ID",trait)
colnames(females) <- c("ID",trait)
#for complete table
colnames(complete) <- c("ID",trait)

#format ID column
males$ID <- as.character(males$ID)
females$ID <- as.character(females$ID)
#for complete table
complete$ID <- as.character(complete$ID)

#read the sorted sample list
sorted_ids <- read.table(sorted_file,header=T,sep=" ")
sorted_ids$ID_2 <- as.character(sorted_ids$ID_2)

#merge males and females datasets
pheno_merged_males <- merge(sorted_ids,males,all=T,by.x="ID_2",by.y="ID")
pheno_merged_males_females <- merge(pheno_merged_males,females,all=T,by.x="ID_2",by.y="ID")
#for complete table
#pheno_merged_males_females <- merge(sorted_ids,complete,all=T,by.x="ID_2",by.y="ID")

#create the merged column
source("~/Work/r_scripts/add_column.r")
pheno_merged_males_females$TRAIT <- apply(pheno_merged_males_females,1,function(row) pheno_add(row[3],row[4]))
pheno_merged_males_females$TRAIT <- as.numeric(as.character(pheno_merged_males_females$TRAIT))
#for complete table
#pheno_merged_males_females$TRAIT <- as.numeric(as.character(pheno_merged_males_females[,3]))

#create the final output file
pheno_merged <- data.frame(pheno_merged_males_females$ID_2,pheno_merged_males_females$TRAIT)
str(pheno_merged)

colnames(pheno_merged) <- c("ID",trait)
pheno_merged$ID <- as.character(pheno_merged$ID)

write.table(pheno_merged,file=out_file,row.names=F,quote=F,col.names=T,sep="\t")