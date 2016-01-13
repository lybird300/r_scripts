#define finction for pvalue calculation
snp_chi_p <- function(data,col1,col2,col3,col4,col_name){
  for(i in 1:dim(data)[1]){
    A <- matrix(c(col1[i], col2[i], col3[i], col4[i]), nrow=2, ncol=2, byrow=T)
    pippo <- chisq.test(A)
    data$P_new[i] <- pippo$p.value
  }
  colnames(data)[colnames(data) %in% "P_new"] <- col_name
  return (data)
}

# define function to flip alleles
flip_all <- function(all){
  if(all == "A"){
    swap <- "T"
  }
  if(all == "T"){
    swap <- "A"
  }
  if(all == "C"){
    swap <- "G"
  }
  if(all == "G"){
    swap <- "C"
  }
  return(swap)
}

################################ FUNCTION DEF END ######################################################
#
#script to assess significance on freq differencies
rm(list=ls())
library(buRlo)

############################## STEP 1: load data from tab files formatted as:
#                              CHROM\tPOS\tAN\tAC\tEAS_AF\tAMR_AF\tAFR_AF\tSAS_AF\tEUR_AF\tAF #####################
# fgrep -h -v CHROM *.tab | sort -g -k1,1 -k2,2 | awk '$5!="." && $6!="." && $7!="." && $8!="." && $9!="."' | awk '{if ($4!~",") print $0}' > CARL_all_freq.tab
# fgrep -h -v CHROM *.tab | sort -g -k1,1 -k2,2 | awk '$5!="." && $6!="." && $7!="." && $8!="." && $9!="."' | awk '{if ($4!~",") print $0}' > FVG_all_freq.tab
# fgrep -h -v CHROM *.tab | sort -g -k1,1 -k2,2 | awk '$5!="." && $6!="." && $7!="." && $8!="." && $9!="."' | awk '{if ($4!~",") print $0}' > VBI_all_freq.tab
cohorts <- c("VBI","FVG","CARL")
cohort <- "VBI"

ingi_cohort <- read.table(paste("/lustre/scratch113/projects/esgi-vbseq/12122015_NICOLA_PLOT/",cohort,"/",cohort,"_all_freq.tab",sep=""),header=F,sep="\t")
colnames(ingi_cohort) <- c("CHROM","POS","AN","AC","EAS_AF","AMR_AF","AFR_AF","SAS_AF","EUR_AF","AF")


# We need to add AN for each TGPphIII cohort
ingi_cohort$AN_AFR <- 661*2 
ingi_cohort$AC_AFR <- ingi_cohort$AN_AFR * ingi_cohort$AFR_AF
ingi_cohort$AN_AMR <- 347*2 
ingi_cohort$AC_AMR <- ingi_cohort$AN_AMR * ingi_cohort$AMR_AF
ingi_cohort$AN_EAS <- 504*2 
ingi_cohort$AC_EAS <- ingi_cohort$AN_EAS * ingi_cohort$EAS_AF
ingi_cohort$AN_EUR <- 503*2 
ingi_cohort$AC_EUR <- ingi_cohort$AN_EUR * ingi_cohort$EUR_AF
ingi_cohort$AN_SAS <- 489*2 
ingi_cohort$AC_SAS <- ingi_cohort$AN_SAS * ingi_cohort$SAS_AF

ingi_cohort$AN_TGP <- 661*2 + 347*2 + 504*2 + 503*2 + 489*2 
ingi_cohort$AC_TGP <- apply(ingi_cohort[,c("AC_AFR","AC_AMR","AC_EAS","AC_EUR","AC_SAS")],1,sum)

#calculate freq differences
ingi_cohort$DIFF_AFR <- ingi_cohort$AF - ingi_cohort$AFR_AF
ingi_cohort$DIFF_AMR <- ingi_cohort$AF - ingi_cohort$AMR_AF
ingi_cohort$DIFF_EAS <- ingi_cohort$AF - ingi_cohort$EAS_AF
ingi_cohort$DIFF_EUR <- ingi_cohort$AF - ingi_cohort$EUR_AF
ingi_cohort$DIFF_SAS <- ingi_cohort$AF - ingi_cohort$SAS_AF

ingi_cohort$TGP_AF <- ingi_cohort$AC_TGP / ingi_cohort$AN_TGP
ingi_cohort$DIFF_TGP <- ingi_cohort$AF - ingi_cohort$TGP_AF

save(ingi_cohort,file=paste(cohort,"_AF.RData",sep=""))
load(paste(cohort,"_AF.RData",sep=""))

#create subset of datasets to check for enrichment
# o_cohorts <- c("AFR","AMR","EAS","EUR","SAS","TGP")

# for(o_cohort in o_cohorts){
#   current_o_cohort <- 
# }

ingi_cohort_AFR <- ingi_cohort[which(ingi_cohort$DIFF_AFR > 0),]
ingi_cohort_AMR <- ingi_cohort[which(ingi_cohort$DIFF_AMR > 0),]
ingi_cohort_EAS <- ingi_cohort[which(ingi_cohort$DIFF_EAS > 0),]
ingi_cohort_EUR <- ingi_cohort[which(ingi_cohort$DIFF_EUR > 0),]
ingi_cohort_SAS <- ingi_cohort[which(ingi_cohort$DIFF_SAS > 0),]
ingi_cohort_TGP <- ingi_cohort[which(ingi_cohort$DIFF_TGP > 0),]

#select only differences bigger tha 1%
ingi_cohort_gt1_AFR <- ingi_cohort[which(ingi_cohort$DIFF_AFR > 0.01),]
ingi_cohort_gt1_AMR <- ingi_cohort[which(ingi_cohort$DIFF_AMR > 0.01),]
ingi_cohort_gt1_EAS <- ingi_cohort[which(ingi_cohort$DIFF_EAS > 0.01),]
ingi_cohort_gt1_EUR <- ingi_cohort[which(ingi_cohort$DIFF_EUR > 0.01),]
ingi_cohort_gt1_SAS <- ingi_cohort[which(ingi_cohort$DIFF_SAS > 0.01),]
ingi_cohort_gt1_TGP <- ingi_cohort[which(ingi_cohort$DIFF_TGP > 0.01),]

dim(ingi_cohort_gt1_AFR)
dim(ingi_cohort_gt1_AMR)
dim(ingi_cohort_gt1_EAS)
dim(ingi_cohort_gt1_EUR)
dim(ingi_cohort_gt1_SAS)
dim(ingi_cohort_gt1_TGP)
#calculate pvalue for each site for each difference
ingi_cohort_gt1_AFR <- snp_chi_p(ingi_cohort_gt1_AFR,ingi_cohort_gt1_AFR$AC_AFR, ingi_cohort_gt1_AFR$AN_AFR, ingi_cohort_gt1_AFR$AC, ingi_cohort_gt1_AFR$AN,"P_AFR")
ingi_cohort_gt1_AMR <- snp_chi_p(ingi_cohort_gt1_AMR,ingi_cohort_gt1_AMR$AC_AMR, ingi_cohort_gt1_AMR$AN_AMR, ingi_cohort_gt1_AMR$AC, ingi_cohort_gt1_AMR$AN,"P_AMR")
ingi_cohort_gt1_EAS <- snp_chi_p(ingi_cohort_gt1_EAS,ingi_cohort_gt1_EAS$AC_EAS, ingi_cohort_gt1_EAS$AN_EAS, ingi_cohort_gt1_EAS$AC, ingi_cohort_gt1_EAS$AN,"P_EAS")
ingi_cohort_gt1_EUR <- snp_chi_p(ingi_cohort_gt1_EUR,ingi_cohort_gt1_EUR$AC_EUR, ingi_cohort_gt1_EUR$AN_EUR, ingi_cohort_gt1_EUR$AC, ingi_cohort_gt1_EUR$AN,"P_EUR")
ingi_cohort_gt1_SAS <- snp_chi_p(ingi_cohort_gt1_SAS,ingi_cohort_gt1_SAS$AC_SAS, ingi_cohort_gt1_SAS$AN_SAS, ingi_cohort_gt1_SAS$AC, ingi_cohort_gt1_SAS$AN,"P_SAS")
ingi_cohort_gt1_TGP <- snp_chi_p(ingi_cohort_gt1_TGP,ingi_cohort_gt1_TGP$AC_TGP, ingi_cohort_gt1_TGP$AN_TGP, ingi_cohort_gt1_TGP$AC, ingi_cohort_gt1_TGP$AN,"P_TGP")

#we want to assign a category to each increment
# we can define 1%, 2%,3%,5%,10%,>10%
o_cohorts <- c("AFR","AMR","EAS","EUR","SAS","TGP")

freqs <- c(0.01,0.02,0.05,0.10,0.10)

for(o_cohort in o_cohorts){
  # o_cohort <- "EUR"
  print(o_cohort)
  current_ingi <- get(paste("ingi_cohort_",o_cohort,sep=""))
  colname_diff <- paste("DIFF_",o_cohort,sep="")

  current_ingi_all_bin <- NULL
  for(i in 1:length(freqs)){
    if(i==1){
      current_ingi_current_bin <- current_ingi[which(current_ingi[,colname_diff] <= freqs[i]),]
    }else if(i==length(freqs)){
      current_ingi_current_bin <- current_ingi[which(current_ingi[,colname_diff] > freqs[i]),]
    }else{
      current_ingi_current_bin <- current_ingi[which(current_ingi[,colname_diff] > freqs[i-1] & current_ingi[,colname_diff] <= freqs[i]),]
    }
  current_ingi_current_bin_resume <- as.data.frame(cbind(bin=freqs[i],tot_ingi=length(current_ingi_current_bin$POS)))

  current_ingi_all_bin <- rbind(current_ingi_all_bin,current_ingi_current_bin_resume)
  assign(paste("ingi_cohort_",o_cohort,"_bins",sep=""),current_ingi_all_bin)

  }
}



write.table(ingi_cohort_all,file=paste(cohort,"_diff_TGP_pval.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")

############################ LOAD EVERYTHING FROM EXOME CHIP WITHOUT SPLITTING BY VILLAGE ###############
#now load only the Exome chip data (80177 markers filtered by frequency/monomorphic, ids updated)
load('/nfs/servizio/FVG_exome.data')

##extract plink format data for exome
export.plink(dati,"whole_exome_chip",transpose=TRUE)

#calculate frequency outside R, and upload data ONLY FOR ALREADY SIGNIFICATIVE SITES FROM PREVIOUS ANALYSIS
chip <- "all_exome"
FVG_data <- read.table("/home/cocca/breast_cancer/ANTI_TUMORAL_DRUGS/exome_snp_list.txt",header=TRUE,sep="\t")
##############################################################################################################

########################### LOAD EVERYTHING FROM CARLANTINO ###############
cohort <- "CARL"
#load chip data from CARLANTINO cohort
load('/nfs/servizio/CARL.data')
chip <- "all_carl_complete"


#now load only the Exome chip data (80177 markers filtered by frequency/monomorphic, ids updated)
load('/nfs/servizio/CARL_exome.data')
chip <- "all_carl_exome_complete"

############################################################################################################
# extract data matching given rsIDs
dati2 <- (dati[,which(snpnames(dati) %in% breast$ID)])
##extract plink format data
export.plink(dati2,paste(cohort,chip,sep="_"),transpose=TRUE)
export.plink(dati,paste(cohort,chip,sep="_"),transpose=TRUE)

#calculate frequency outside R, and upload data ONLY FOR ALREADY SIGNIFICATIVE SITES FROM PREVIOUS ANALYSIS
FVG_data <- read.table("/home/cocca/genotypes/CARL/ARRAY/CARL_chip.freq.frq",header=TRUE,sep="\t")
##############################################################################################################


#calculate frequency outside R, and upload data FOR ALL EXOME SITES 
chip <- "all_exome_complete"
FVG_data <- read.table("/home/cocca/genotypes/CARL/EXOME/CARL_chip.freq.frq",header=TRUE,sep="\t")
#Upload ALL data from BREAST CANCER TARGET SMPS with TGP info
breast <- read.table("/home/cocca/breast_cancer/ANTI_TUMORAL_DRUGS/breast_target_snp_TGP_data.txt",header=TRUE,sep=" ")
# add AC count for EUR
breast$EUR_AC <- (breast$EUR_AN*breast$EUR_AF)
##############################################################################################################

#merge dataset
data_merged <- merge(breast,FVG_data,by.x='ID',by.y='SNP',sort=F)

# cast allele type into string
data_merged$REF <- as.character(data_merged$REF)
data_merged$ALT <- as.character(data_merged$ALT)
data_merged$A1 <- as.character(data_merged$A1)
data_merged$A2 <- as.character(data_merged$A2)

# we need to work for each row
for(n in 1:length(data_merged$ID)){
  if (data_merged$ALT[n] == data_merged$A1[n]){
    data_merged$VILL_CHK_ALL[n] <- data_merged$A1[n]
    data_merged$VILL_CHK_AF[n] <- data_merged$MAF[n]
  }else if (data_merged$ALT[n] == data_merged$A2[n]) {
    data_merged$VILL_CHK_ALL[n] <- data_merged$A2[n]
    data_merged$VILL_CHK_AF[n] <- 1-data_merged$MAF[n]
  }else if (data_merged$ALT[n] == flip_all(data_merged$A1[n])){
    data_merged$VILL_CHK_ALL[n] <- flip_all(data_merged$A1[n])
    data_merged$VILL_CHK_AF[n] <- data_merged$MAF[n]
  }else if (data_merged$ALT[n] == flip_all(data_merged$A2[n])){
    data_merged$VILL_CHK_ALL[n] <- flip_all(data_merged$A2[n])
    data_merged$VILL_CHK_AF[n] <- 1-data_merged$MAF[n]
  }
}

# add correct allele count
data_merged$VILL_CHK_AC <- data_merged$VILL_CHK_AF*data_merged$NCHROBS

# now we can calculate and add pvalue column
data_merged_pval <- snp_chi_p(data_merged,data_merged$EUR_AC, data_merged$EUR_AN, data_merged$VILL_CHK_AC, data_merged$NCHROBS,"data_P")

#write table
write.table(data_merged_pval,file=paste("tabella_geni",cohort,chip,"pval.txt",sep="_"),quote=F,row.names=F,col.names=T,sep="\t")

############################################################################################
# Now upload the data from different chips splitted by village and create datasets for each village, ordering with the corrected allele

# 1) 700K
chip <- "omni_700"
FVG_data <- read.table("/home/cocca/genotypes/VILLAGES/700K/All_villages.sig_freq.list",header=TRUE,sep="\t")
# 2) 300K+700K
chip <- "join_300_700"
FVG_data <- read.table("/home/cocca/genotypes/VILLAGES/JOIN_300_700/All_villages.sig_freq.list",header=TRUE,sep="\t")
# 3) EXOME
chip <- "exome"
FVG_data <- read.table("/home/cocca/genotypes/VILLAGES/EXOME/All_villages.sig_freq.list",header=TRUE,sep="\t")
# 3) EXOME ALL
chip <- "exome_all"
FVG_data <- read.table("/home/cocca/genotypes/VILLAGES/EXOME/All_villages.exome.sig_freq.list",header=TRUE,sep="\t")

##########################################################################################

# create different dataset for each village
villages <- c("clauzetto","smc","erto","illegio","resia","sauris")
base_folder <- "/home/cocca/genotypes/EXOME_CHIP/VILLAGES/"
# villages <- unique(FVG_data$village)

for (village in villages){
  # current_village <- FVG_data[which(FVG_data$village == village),]
  current_village <- read.table(paste(base_folder,village,".freq.frq",sep=""),header=T)
  village_dataset <- paste(chip,village,sep="_")
  # merge data with EUR/TGP data
  current_merged <- merge(breast,current_village,by.x='ID',by.y='SNP',sort=F)
  # cast allele type into string
  current_merged$REF <- as.character(current_merged$REF)
  current_merged$ALT <- as.character(current_merged$ALT)
  current_merged$A1 <- as.character(current_merged$A1)
  current_merged$A2 <- as.character(current_merged$A2)
  # we need to work for each row
  for(n in 1:length(current_merged$ID)){
    if (current_merged$ALT[n] == current_merged$A1[n]){
      current_merged$VILL_CHK_ALL[n] <- current_merged$A1[n]
      current_merged$VILL_CHK_AF[n] <- current_merged$MAF[n]
    }else if (current_merged$ALT[n] == current_merged$A2[n]) {
      current_merged$VILL_CHK_ALL[n] <- current_merged$A2[n]
      current_merged$VILL_CHK_AF[n] <- 1-current_merged$MAF[n]
    }else if (current_merged$ALT[n] == flip_all(current_merged$A1[n])){
      current_merged$VILL_CHK_ALL[n] <- flip_all(current_merged$A1[n])
      current_merged$VILL_CHK_AF[n] <- current_merged$MAF[n]
    }else if (current_merged$ALT[n] == flip_all(current_merged$A2[n])){
      current_merged$VILL_CHK_ALL[n] <- flip_all(current_merged$A2[n])
      current_merged$VILL_CHK_AF[n] <- 1-current_merged$MAF[n]
    }
  }
  # add correct allele count
  current_merged$VILL_CHK_AC <- current_merged$VILL_CHK_AF*current_merged$NCHROBS

  # now we can calculate and add pvalue column
  current_merged_pval <- snp_chi_p(current_merged,current_merged$EUR_AC, current_merged$EUR_AN, current_merged$VILL_CHK_AC, current_merged$NCHROBS,paste(village,"P",sep="_"))
  # sort by pval
  # current_merged_pval <- current_merged_pval[order(current_merged_pval$P),]

  write.table(current_merged_pval,file=paste("/home/cocca/breast_cancer/RESULTS/20140108_RESULTS/tabella_geni",village,chip,"pval.txt",sep="_"),quote=F,row.names=F,col.names=T,sep="\t")
  assign(village_dataset,current_merged_pval)

}
