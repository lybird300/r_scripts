#script to estract statisics fron vcf-stats results	and from others for FVG pop
#TODO: parametrize all this!
rm(list=ls())

#read commandline args
args <- commandArgs(trailing=TRUE)

#Args list:
#args[[1]]:list with individuals and sex (absolute path) (i.e '/nfs/users/nfs_m/mc14/Work/SANGER/DOCS/ValBorbera/sequenced_id_conversion_sex.txt')
#args[[2]]:the path of vcf-stats files to process (absolute path)
#args[[3]]: prefix used to assign names to vcf-stats output files

#Assign args
#indiv_file <- args[[1]]
#vcfstats_filepath <- args[[2]]
#vcfstats_file_prefix <- args[[3]]
indiv_file <- '/nfs/users/nfs_m/mc14/Work/SANGER/DOCS/ValBorbera/sequenced_id_conversion_sex.txt'
vcfstats_filepath <- 'vcf_stats'
vcfstats_file_prefix <- 'generic_stats'

#upload list with sex
individuals <- read.table(indiv_file, header=F)

#set colnames based on col number
if ( dim(individuals)[2] == 2 ){
colnames(individuals) <- c('SEQ_ID','sex')
}

if ( dim(individuals)[2] == 3 ){
colnames(individuals) <- c('SEQ_ID','PRIVATE_ID','sex')
}
#We assume by default that the last column is the sex column,
# and it is coded as follow:
#MALE:1
#FEMALE:2

#separate male/female
male_inds <- individuals[which(individuals$sex == 1),]
female_inds <- individuals[which(individuals$sex == 2),]

#Open a file stream to write some outputs
#sink('stats_resumes.log')

#INDELS
indels <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'indels',sep="."),sep="/"), header=F, skip=2,col.names=c('count','sample'))
indels$sample <- gsub("samples/","",indels$sample)

indels_summary <- summary(indels)

#cat('INDELS\n')
#cat(indels_summary,labels=names(indels_summary))
#cat('\n\n')

write.table(indels_summary, file='indels_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("indels_density.jpg",width=1000, height=1000)
plot(density(indels$count),main ="Indels density distribution in VBSEQ", xlab="Indels count")
dev.off()

jpeg("indels_boxplot.jpg", width=1000, height=1000)
indels_boxplot <- boxplot(indels$count,main ="Indels in VBSEQ", xlab="All individuals", ylab="Indels count")
dev.off()

#cat('SNPs Outliers\n')
#cat(indels[which(indels$count %in%  indels_boxplot$out),]$sample)
#cat('\n\n')

#to retrieve the outliers from boxplot
indels_outliers <- indels[which(indels$count %in% indels_boxplot$out),]$sample

#male/female boxplot
jpeg("indels_boxplot_male_female.jpg", width=1000, height=1000)
indels_boxplot_mf <- boxplot(indels[which(indels$sample %in% male_inds$SEQ_ID),]$count,indels[which(indels$sample %in% female_inds$SEQ_ID),]$count,main ="Indels in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Indels count")
dev.off()

#to retrieve the outliers from boxplot
#male_indels_outliers <- indels[which(indels$sample %in% male_inds$SEQ_ID && indels$count %in% indels_boxplot$out),]$sample
#female_indels_outliers <- indels[which(indels$count %in% indels_boxplot$out),]$sample

#######PRIVATE SNPS
private_s <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'private',sep="."),sep="/"), header=F, skip=1,col.names=c('count','sample'))

private_s_summary <- summary(private_s)

#cat('PRIVATE SNPS\n')
#cat(private_s_summary)
#cat('\n\n')

write.table(private_s_summary, file='private_snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("private_snps_density.jpg", width=1000, height=1000)
plot(density(private_s$count),main ="Private SNPs density distribution in VBSEQ", xlab="Private SNPs count")
dev.off()

jpeg("private_snps_boxplot.jpg", width=1000, height=1000)
private_s_boxplot <- boxplot(private_s$count,main ="Private SNPs in VBSEQ", xlab="All individuals", ylab="Private SNPs count")
dev.off()

#to retrieve the outliers from boxplot
private_outliers <- private_s[which(private_s$count %in%  private_s_boxplot$out),]$sample

#cat('PRIVATE SNPs Outliers\n')
#cat(private_s[which(private_s$count %in%  private_s_boxplot$out),]$sample)
#cat('\n\n')
#sink()

#male/female boxplot
jpeg("private_snps_boxplot_male_female.jpg", width=1000, height=1000)
private_s_boxplot_mf <- boxplot(private_s[which(private_s$sample %in% male_inds$SEQ_ID),]$count,private_s[which(private_s$sample %in% female_inds$SEQ_ID),]$count,main ="Private SNPs in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Private SNPs count")
dev.off()

################SHARED SNPS
#Number of sites having a non-reference allele in 0,1,2,etc samples

shared_s <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'shared',sep="."),sep="/"), header=F, skip=1,col.names=c('ind_count','freq'))

shared_s_summary <- summary(shared_s$freq)

#write.table(shared_s_summary, file='shared_snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("shared_snps_density.jpg", width=1000, height=1000)
plot(density(shared_s$freq),main ="Density distribution of the number of sites having a non-reference allele in VBSEQ",xlab="SNPs count")
dev.off()

jpeg("shared_snps_hist.jpg", width=1000, height=1000)
barplot(shared_s$freq,names.arg=shared_s$ind_count,main ="Number of sites having a non-reference allele in VBSEQ", xlab="# Individuals", ylab="Number of sites having a non-reference allel (count)")
dev.off()

#SNPS
#Number of positions with SNPs

snps <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'snps',sep="."),sep="/"), header=F, skip=2,col.names=c('count','sample'))
snps$sample <- gsub("samples/","",snps$sample)

snps_summary <- summary(snps$count)

#write.table(snps_summary, file='snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("snps_density.jpg", width=1000, height=1000)
plot(density(snps$count),main ="SNPs density distribution in VBSEQ",xlab="SNPs count")
dev.off()

jpeg("snps_boxplot.jpg", width=1000, height=1000)
snps_boxplot <- boxplot(snps$count,main ="SNPs in VBSEQ", xlab="All individuals", ylab="SNPs count")
dev.off()

#to retrieve the outliers from boxplot
snps_outliers <- snps[which(snps$count %in%  snps_boxplot$out),]$sample

#male/female boxplot
jpeg("snps_boxplot_male_female.jpg", width=1000, height=1000)
snps_boxplot_mf <- boxplot(snps[which(snps$sample %in% male_inds$SEQ_ID),]$count,snps[which(snps$sample %in% female_inds$SEQ_ID),]$count,main ="SNPs in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="SNPs count")
dev.off()

################SITES by SAMPLE stats
persample_file_path <- "vcf_check/per_sample_count.chk"
persample <- read.table(persample_file_path,header=F,skip=2)

colnames(persample) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions","nIndels","average_depth","nSingletons")

#read sample resume file
sample_resume_file <- "~/Work/SANGER/FVG/SEQ_CALLING/DOCS/samples_resume.txt"
sample_resume <- read.table(sample_resume_file,header=T)
sample_stats <- read.table("~/Work/SANGER/FVG/SEQ_CALLING/DOCS/samples_sites_stats.txt",header=T)

#merge together
sample_complete_resume <- merge(sample_stats,sample_resume, by.x="sample",by.y="called_id") #no "all" option because we want the merged on the sequenced only
# sample_complete_resume <- merge(persample,sample_resume, by.x="sample",by.y="called_id") #no "all" option because we want the merged on the sequenced only

# sample_complete_resume$nSites <- apply(sample_complete_resume,c("nRefHom","nNonRefHom","nHets"),sum)
# sample_complete_resume$nSites <- sum(sample_complete_resume$nRefHom,sample_complete_resume$nNonRefHom,sample_complete_resume$nHets)

##Plot snp Number
#plot stratifiyng by village
jpeg(paste(pop_name,"snps_count_plot_by_village.jpg",sep="_"),width=800, height=600)
plot(sample_complete_resume[which(sample_complete_resume$village == "Erto"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Erto"),]$n_snps,main="N SNPs by samples",xlab="Samples",ylab="N Snps",col="purple",ylim=c(min(sample_complete_resume$n_snps),3350000),xlim=c(0,300),pch=19)
points(sample_complete_resume[which(sample_complete_resume$village == "Illegio"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Illegio"),]$n_snps,col="blue",pch=19)
points(sample_complete_resume[which(sample_complete_resume$village == "Resia"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Resia"),]$n_snps,col="red",pch=19)
points(sample_complete_resume[which(sample_complete_resume$village == "Sauris"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Sauris"),]$n_snps,col="green",pch=19)
legend(250,3350000,legend=c("Erto","Illegio","Resia","Sauris"),fill=c("purple","blue","red","green"))
dev.off()

#stratify by sex
jpeg(paste(pop_name,"snps_count_plot_by_sex.jpg",sep="_"),width=800, height=600)
plot(sample_complete_resume[which(sample_complete_resume$sex == "Male"),]$index,sample_complete_resume[which(sample_complete_resume$sex == "Male"),]$n_snps,main="N SNPs by samples by sex",xlab="Samples",ylab="N SNPs",col="darkblue",ylim=c(min(sample_complete_resume$n_snps),3350000),xlim=c(0,300),pch=19)
points(sample_complete_resume[which(sample_complete_resume$sex == "Female"),]$index,sample_complete_resume[which(sample_complete_resume$sex == "Female"),]$n_snps,col="darkred",pch=19)
legend(250,3350000,legend=c("Males","Females"),fill=c("darkblue","darkred"))
dev.off()


##Plot indels number
#plot stratifiyng by village
jpeg(paste(pop_name,"indels_count_plot_by_village.jpg",sep="_"),width=800, height=600)
plot(sample_complete_resume[which(sample_complete_resume$village == "Erto"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Erto"),]$n_indels,main="N indels by samples",xlab="Samples",ylab="N indels",col="purple",ylim=c(min(sample_complete_resume$n_indels),570000),xlim=c(0,300),pch=19)
points(sample_complete_resume[which(sample_complete_resume$village == "Illegio"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Illegio"),]$n_indels,col="blue",pch=19)
points(sample_complete_resume[which(sample_complete_resume$village == "Resia"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Resia"),]$n_indels,col="red",pch=19)
points(sample_complete_resume[which(sample_complete_resume$village == "Sauris"),]$index,sample_complete_resume[which(sample_complete_resume$village == "Sauris"),]$n_indels,col="green",pch=19)
legend(250,570000,legend=c("Erto","Illegio","Resia","Sauris"),fill=c("purple","blue","red","green"))
dev.off()

#stratify by sex
jpeg(paste(pop_name,"indels_site_count_plot_by_sex.jpg",sep="_"),width=800, height=600)
plot(sample_complete_resume[which(sample_complete_resume$sex == "Male"),]$index,sample_complete_resume[which(sample_complete_resume$sex == "Male"),]$n_indels,main="N indels by samples by sex",xlab="Samples",ylab="N indels",col="darkblue",ylim=c(min(sample_complete_resume$n_indels),570000),xlim=c(0,300),pch=19)
points(sample_complete_resume[which(sample_complete_resume$sex == "Female"),]$index,sample_complete_resume[which(sample_complete_resume$sex == "Female"),]$n_indels,col="darkred",pch=19)
legend(250,570000,legend=c("Males","Females"),fill=c("darkblue","darkred"))
dev.off()

###############SAMPLE TS/TV
#samples-tstv

s_tstv <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'samples-tstv',sep="."),sep="/"), header=F,skip=1,col.names=c('Transitions','Transversions','ts_tv_ratio','sample'))

s_tstv_summary <- summary(s_tstv)

#write.table(s_tstv_summary, file='s_tstv_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("s_ts_density.jpg", width=1000, height=1000)
plot(density(s_tstv$Transitions),col='red',main ="Transitions density distribution in VBSEQ",xlab="Transition count")
dev.off()

jpeg("s_tv_density.jpg", width=1000, height=1000)
plot(density(s_tstv$Transversions),col='green',main ="Transversions density distribution in VBSEQ",xlab="Transversions count")
dev.off()

jpeg("s_tstv_density.jpg", width=1000, height=1000)
plot(density(s_tstv$ts_tv_ratio),col='blue',main ="Ts/Tv ratio density distribution in VBSEQ",xlab="Ts/Tv ratio")
dev.off()

#boxplots
jpeg("s_ts_boxplot.jpg", width=1000, height=1000)
s_ts_boxplot <- boxplot(s_tstv$Transitions,main ="Transition in VBSEQ",ylab="Transition count", xlab="All individuals", names=c("Transitions"))
dev.off()

jpeg("s_tv_boxplot.jpg", width=1000, height=1000)
s_tv_boxplot <- boxplot(s_tstv$Transversions,main ="Transversion in VBSEQ",ylab="Transversion count", xlab="All individuals", names=c("Transversions"))
dev.off()

jpeg("s_tstv_boxplot.jpg", width=1000, height=1000)
s_tstv_boxplot <- boxplot(s_tstv$ts_tv_ratio,main ="Ts/Tv ratio in VBSEQ",ylab="Ts/Tv ratio", xlab="All individuals", names=c("Ts/Tv"))
dev.off()

#male/female boxplot
jpeg("s_ts_boxplot_male_female.jpg", width=1000, height=1000)
s_ts_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$Transitions,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$Transitions,main ="Transitions in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Transitions count")
dev.off()

jpeg("s_tv_boxplot_male_female.jpg", width=1000, height=1000)
s_tv_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$Transversions,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$Transversions,main ="Transversions in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Transversions count")
dev.off()

jpeg("s_tstv_boxplot_male_female.jpg", width=1000, height=1000)
s_tstv_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$ts_tv_ratio,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$ts_tv_ratio,main ="Ts/Tv ratio in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Ts/Tv ratio")
dev.off()

#to retrieve the outliers from boxplot
s_ts_outliers <- s_tstv[which(s_tstv$Transitions %in%  s_ts_boxplot$out),]$sample
s_tv_outliers <- s_tstv[which(s_tstv$Transversions %in%  s_tv_boxplot$out),]$sample
s_tstv_outliers <- s_tstv[which(s_tstv$ts_tv_ratio %in%  s_tstv_boxplot$out),]$sample

###write a file for outliers
sink('Outliers.txt')
cat('Outliers for INDELS\n')
cat(indels_outliers)
cat('\n\nOutliers for PRIVATE SNPs\n')
cat(as.character(private_outliers))
cat('\n\nOutliers for SNPs\n')
cat(snps_outliers)
cat('\n\nOutliers for Transitions\n')
cat(as.character(s_ts_outliers))
cat('\n\nOutliers for Transversions\n')
cat(as.character(s_tv_outliers))
cat('\n\nOutliers for Ts/Tv ratio\n')
cat(as.character(s_tstv_outliers))
cat('\n')
sink()

######## STATS and plots for post calling QC and reports ###################
#SAMPLE LEVEL:
#
#first remove remove multiallelic sites
#read table with some info:
# all_variants_table_path <- "/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vqslod.fixed.no_geno.vcf.gz.csv"
# all_variants_table_path <- "/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vqslod.fixed.noX.UNRELATED.no_geno.vcf.gz.maf_table.tab"
all_variants_table_path <- "/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vqslod.fixed.noX.UNRELATED.no_geno.vcf.gz.maf_table_snp.tab"
variants <- read.table(all_variants_table_path,sep=" ", header=T)
# colnames(variants) <- c("CHROM","POS","REF","ALT","TYPE","AC","AN","TGP_AF","TGP_EUR","DP","IMP2","QBD","HWE","ICF","VQSLOD","AF","MAF","TGP_MAF","TGP_EUR_MAF")

variants$AC <- as.integer(as.character(variants$AC))
variants$TGP_AF <- as.numeric(as.character(variants$TGP_AF))
variants$TGP_EUR <- as.numeric(as.character(variants$TGP_EUR))
variants$DP <- as.integer(as.character(variants$DP))
variants$IMP2 <- as.character(variants$IMP2)
variants$QBD <- as.numeric(as.character(variants$QBD))
variants$HWE <- as.numeric(as.character(variants$HWE))
variants$ICF <- as.numeric(as.character(variants$ICF))
variants$VQSLOD <- as.numeric(as.character(variants$VQSLOD))
variants$AF <- as.numeric(as.character(variants$AF))
variants$MAF <- as.numeric(as.character(variants$MAF))
variants$TGP_MAF <- as.numeric(as.character(variants$TGP_MAF))
variants$TGP_EUR_MAF <- as.numeric(as.character(variants$TGP_EUR_MAF))

#summary for data
summary(variants)

#plot maf seq vs maf TGP/EUR
jpeg("AF_FVGvsTGP_EUR.jpg",width=1754, height=1024)
  par(cex=3)
  plot(variants$TGP_EUR,variants$AF, xlab="AF EUR", ylab="AF FVG")
dev.off()

#split by variant type (if you didn't do it at the beginning)
snps <- variants[which(variants$TYPE == "SNP"),]
indels <- variants[which(variants$TYPE == "INDEL"),]

#summaries
summary(snps)
summary(indels)

#Generate bar plots/density for VQSLOD,depth,maf,HWE,impute2 info score
#

#Generate manhattan plots by chr for VQSLOD,depth,maf,HWE,impute2 info score
maf_seq <-read.table("/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1_stats/PLINK/GWAS_OVERLAP/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.sex.GWAS.sorted.filter_right_sex_cleaned_freq.frq", header=T)
maf_gwas <- read.table("/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/MERGED/AUTO/GW_SEQ_OVER/FVG_seq_subset_MERGED_sorted_auto_overlap_right_sex.sorted_filter_cleaned_freq.frq",header=T)
#screenshot for the picture...
plot(maf_seq$MAF,maf_gwas$MAF,col="red")
#merge sequence and gwas to calculate DMAF
merged_maf <- merge(maf_seq,maf_gwas,by=c("CHR","SNP"),all=T)
merged_maf$DMAF <- abs(merged_maf$MAF.x - merged_maf$MAF.y)

#load a file with variants positions
pos_ref <- read.table("/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1_stats/PLINK/GWAS_OVERLAP/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.sex.GWAS.sorted.filter_right_sex_cleaned.bim",header=F)
colnames(pos_ref) <- c("CHR","SNP","cmPOS","POS","A1","A2")

#merge together the last two to have positions...
merged_maf_pos <- merge(merged_maf,pos_ref,by=c("CHR","SNP"),all=T,sort=F)


#plot DMAF by chr
for(i in 1:22){
	jpeg(paste("DMAF_chr",i,".jpg",sep=""),width=877, height=512,pointsize = 20)
		current_merged_maf_pos <- merged_maf_pos[which(merged_maf_pos$CHR == i),]
		plot(current_merged_maf_pos$POS,current_merged_maf_pos$DMAF,col="darkgreen",xlab="position",ylab="DMAF",main=paste("CHR",i,sep=" "))
	dev.off()
}

#upload non nreference discordance tables by sample
discordance_chr_path <- "/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1_stats/DISCORDANCE/QC_OUT"
sample_conc_table_file <- paste(discordance_chr_path,"all_sample_concordance_discordance.txt",sep="/")
sample_conc <- read.table(sample_conc_table_file,header=T)
sample_conc$SAMPLE_ID <- as.character(sample_conc$SAMPLE_ID)

#print some stats
jpeg("NRD_by_sample.jpg",width=1754, height=1024)
  par(cex=3)
  plot(sample_conc$NRD, xlab="Samples", ylab="NRD")
dev.off()

jpeg("CONC_by_sample.jpg",width=1754, height=1024)
  par(cex=3)
  plot(sample_conc$OGC, xlab="Samples", ylab="OGC")
dev.off()

conc_box <- boxplot(sample_conc$OGC)
nrd_box <- boxplot(sample_conc$NRD)

 #   SAMPLE_ID        OGC              NRD          
 # 582105 :  1   Min.   :0.5120   Min.   :0.009645  
 # 582106 :  1   1st Qu.:0.9922   1st Qu.:0.013253  
 # 582109 :  1   Median :0.9932   Median :0.015637  
 # 582111 :  1   Mean   :0.9880   Mean   :0.024342  
 # 582114 :  1   3rd Qu.:0.9943   3rd Qu.:0.017807  
 # 582115 :  1   Max.   :0.9958   Max.   :0.768571  
 # (Other):243                            

#upload all non nreference discordance tables by chromosomes
for(j in 1:22){
  current_chr_path <- paste(discordance_chr_path,"/CHR",j,sep="")
  current_chr <- read.table(paste(current_chr_path,"/site_concordance_discordance_table_chr",j,".txt",sep=""),header=T)

  jpeg(paste("NRD_chr",j,".jpg",sep=""),width=800, height=800)
    plot(current_chr$POS,current_chr$NRD,col="darkgreen",xlab="position",ylab="NRD",main=paste("CHR",j,sep=" "))
  dev.off()
}

#upload het data
seq_het <- read.table("/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1_stats/PLINK/GWAS_OVERLAP/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.sex.GWAS.sorted.filter_right_sex_cleaned_het.het",header=T)

gwas_het <- read.table("/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/MERGED/AUTO/GW_SEQ_OVER/FVG_seq_subset_MERGED_sorted_auto_overlap_right_sex.sorted_filter_cleaned_het.het",header=T)

#work with chr20
gwas_het_chr20 <- read.table("/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/MERGED/AUTO/GW_SEQ_OVER/FVG_seq_subset_MERGED_sorted_auto_overlap_right_sex.sorted_filter_cleaned_het_chr20.het",header=T)
sample_conc_20 <- read.table("/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1_stats/DISCORDANCE/QC_OUT/CHR20/sample_concordance_discordance_table_chr20.txt",header=T)
sample_conc_20 <- sample_conc_20[,1:4]

#merge het data with concordance
seq_het_conc<- merge(seq_het,sample_conc,by.x="FID",by.y="SAMPLE_ID",all=T)
gwas_het_conc <- merge(gwas_het,sample_conc,by.x="FID",by.y="SAMPLE_ID",all=T)
gwas_het_conc_chr20 <- merge(gwas_het_chr20,sample_conc_20,by.x="FID",by.y="SAMPLE",all=T)

#plot inbreeding coeff vs NRD
jpeg("FvsNRD_by_sample.jpg",width=877, height=877)
  par(cex=3)
  plot(gwas_het_conc$F,gwas_het_conc$NRD, xlab="F inbreeding coefficient", ylab="NR Discordance",col="darkred")
dev.off()

jpeg("FvsNRD_by_sample_chr20.jpg",width=877, height=877)
  par(cex=3)
  plot(gwas_het_conc_chr20$F,gwas_het_conc_chr20$NRD, xlab="F inbreeding coefficient", ylab="NR Discordance",col="darkred")
dev.off()

#plot maf distribution
jpeg("MAF_density.jpg",width=877, height=877)
	plot(density(variants[which(variants$MAF != 0),]$MAF),xlab="MAF")
dev.off()


#work by village
for(i in 1:length(villages)){
current_village <- villages[i]
current_list <- assign(paste(current_village,"list",sep="_"), sample_complete_resume[which(sample_complete_resume$village == current_village),]$sample)
write.table(current_list, file=paste(current_village,"list.txt",sep="_"),append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = FALSE)
}

########################################################################################################
########################### EXTRACT FREQ SPECTRUM STATS ################################################

#read village specific files
#make this parametric
pop_name <- "FVG_auto"
# all_variants_table_path <- "/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vqslod.fixed.no_geno.vcf.gz.csv"
all_variants_table_path <- "/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vqslod.fixed.noX.UNRELATED.no_geno.vcf.gz.maf_table_snp.tab"
current_seq <- read.table(all_variants_table_path,header=TRUE,stringsAsFactors=F, comment.char="",sep=" ")
#remove multiallelic sites
# current_seq_1 <- current_seq[-grep(",",current_seq$ALT),]

#remove X chr
# current_seq <- current_seq[-which(current_seq$CHROM == "X"),]
# colnames(variants) <- c("CHROM","POS","REF","ALT","TYPE","AC","AN","TGP_AF","TGP_EUR","DP","IMP2","QBD","HWE","ICF","VQSLOD","AF","MAF","TGP_MAF","TGP_EUR_MAF")

#convert to numeric
# variants$TGP_AF <- as.numeric(as.character(variants$TGP_AF))
# variants$TGP_EUR <- as.numeric(as.character(variants$TGP_EUR))
# variants$DP <- as.integer(as.character(variants$DP))
# variants$QBD <- as.numeric(as.character(variants$QBD))
# variants$HWE <- as.numeric(as.character(variants$HWE))
# variants$ICF <- as.numeric(as.character(variants$ICF))
# variants$TGP_MAF <- as.numeric(as.character(variants$TGP_MAF))
# variants$TGP_EUR_MAF <- as.numeric(as.character(variants$TGP_EUR_MAF))
current_seq$AC <- as.integer(as.character(current_seq$AC))
current_seq$IMP2 <- as.character(current_seq$IMP2)
current_seq$VQSLOD <- as.numeric(as.character(current_seq$VQSLOD))
current_seq$AF <- as.numeric(as.character(current_seq$AF))
current_seq$MAF <- as.numeric(as.character(current_seq$MAF))
current_seq$TGP_AF <- as.numeric(as.character(current_seq$TGP_AF))
current_seq$TGP_AMR_AF <- as.numeric(as.character(current_seq$TGP_AMR_AF))
current_seq$TGP_ASN_AF <- as.numeric(as.character(current_seq$TGP_ASN_AF))
current_seq$TGP_AFR_AF <- as.numeric(as.character(current_seq$TGP_AFR_AF))
current_seq$TGP_EUR_AF <- as.numeric(as.character(current_seq$TGP_EUR_AF))


#count by maf
# maf_classes <- c(0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
maf_classes <- c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5)

current_village_classes <- paste(pop_name,"maf_classes",sep="_")

#we removed the monomorphic sites
assign(current_village_classes,split_bins(maf_classes,current_seq[which(current_seq$MAF > 0),],pop_name))

# rm(current_seq)
# gc()
#write a cute output
sink(paste(pop_name,'maf_bin_resume.txt',sep="_"))
	print(get(current_village_classes))
sink()

jpeg(paste(pop_name,"site_count_plot.jpg",sep="_"),width=1000, height=1000)
	barplot(as.matrix(get(current_village_classes)),names.arg=maf_classes[2:length(maf_classes)],col=colors()[72])
dev.off()


############################################################################################
############################## PSEUDO-ENRICHMENT CHECK #####################################

#define this only once
EUR_maf_bin <- NULL
TGP_maf_bin <- NULL


#read table for each bin
for(bin in maf_classes[2:length(maf_classes)]){

#select a threshold to define enrichment or depletion based on dmaf percentage of each bin
enrich_thr <- (0.5*bin)
deplet_thr <- -(0.5*bin)

current_table <- paste("/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1_stats/",pop_name,"_maf_lte_",bin,"_table.txt",sep="")
current_table_maf <- read.table(current_table,header=T,stringsAsFactors=F,comment.char="") #sauris_auto_maf_lt3

#convert to numeric
current_table_maf$TGP_AF <- as.numeric(as.character(current_table_maf$TGP_AF))
current_table_maf$AMR_AF <- as.numeric(as.character(current_table_maf$AMR_AF))
current_table_maf$ASN_AF <- as.numeric(as.character(current_table_maf$ASN_AF))
current_table_maf$AFR_AF <- as.numeric(as.character(current_table_maf$AFR_AF))
current_table_maf$EUR_AF <- as.numeric(as.character(current_table_maf$EUR_AF))
current_table_maf$EUR_MAF <- ifelse(current_table_maf$EUR_AF > 0.5,1-current_table_maf$EUR_AF,current_table_maf$EUR_AF)
current_table_maf$EUR_MINOR <- ifelse(current_table_maf$EUR_AF > 0.5,as.character(current_table_maf$REF),as.character(current_table_maf$ALT))
current_table_maf$TGP_MAF <- ifelse(current_table_maf$TGP_AF > 0.5,1-current_table_maf$TGP_AF,current_table_maf$TGP_AF)
current_table_maf$TGP_MINOR <- ifelse(current_table_maf$TGP_AF > 0.5,as.character(current_table_maf$REF),as.character(current_table_maf$ALT))


# calculate DMAF with european
current_table_maf$DMAF_EUR <- (current_table_maf$MAF - current_table_maf$EUR_MAF)
current_table_maf$DMAF_TGP <- (current_table_maf$MAF - current_table_maf$TGP_MAF)


#count for each population, how many sites we have above and below the treshold EUR/TGP
EUR_enriched <- length(current_table_maf[which(current_table_maf$DMAF_EUR >= enrich_thr),]$EUR_MINOR)
EUR_depleted <- length(current_table_maf[which(current_table_maf$DMAF_EUR <= deplet_thr),]$EUR_MINOR)
EUR_stable <- length(current_table_maf[which(current_table_maf$DMAF_EUR > deplet_thr & current_table_maf$DMAF_EUR < enrich_thr),]$EUR_MINOR)
EUR_private <- length(current_table_maf[which(is.na(current_table_maf$EUR_MAF)),]$EUR_MINOR)
EUR_enriched + EUR_stable + EUR_depleted + EUR_private

TGP_enriched <- length(current_table_maf[which(current_table_maf$DMAF_TGP >= enrich_thr),]$TGP_MINOR)
TGP_depleted <- length(current_table_maf[which(current_table_maf$DMAF_TGP <= deplet_thr),]$TGP_MINOR)
TGP_stable <- length(current_table_maf[which(current_table_maf$DMAF_TGP > deplet_thr & current_table_maf$DMAF_TGP < enrich_thr),]$TGP_MINOR)
TGP_private <- length(current_table_maf[which(is.na(current_table_maf$TGP_MAF)),]$TGP_MINOR)
TGP_enriched + TGP_stable + TGP_depleted + TGP_private

dim(current_table_maf)

#calculate median with respect to overlapping sites
my_pop_median_EUR <- median(current_table_maf[-which(is.na(current_table_maf$EUR_MAF)),]$MAF)
my_pop_median_TGP <- median(current_table_maf[-which(is.na(current_table_maf$TGP_MAF)),]$MAF)

EUR_median <- median(current_table_maf[-which(is.na(current_table_maf$EUR_MAF)),]$EUR_MAF)
TGP_median <- median(current_table_maf[-which(is.na(current_table_maf$TGP_MAF)),]$TGP_MAF)

#calculate a ks test for overlapping sites in this bin


#create table with proportion of sites
EUR_current_maf_bin <- c(EUR_enriched,EUR_stable,EUR_depleted,EUR_private)
TGP_current_maf_bin <- c(TGP_enriched,TGP_stable,TGP_depleted,TGP_private)

#define a matrix to plot everything
EUR_maf_bin <- as.matrix(rbind(EUR_maf_bin,EUR_current_maf_bin))

TGP_maf_bin <- as.matrix(rbind(TGP_maf_bin,TGP_current_maf_bin))



}

colnames(EUR_maf_bin) <- c("enriched","stable","depleted","private")
colnames(TGP_maf_bin) <- c("enriched","stable","depleted","private")

rownames(EUR_maf_bin) <- maf_classes[2:length(maf_classes)]
rownames(TGP_maf_bin) <- maf_classes[2:length(maf_classes)]

sink(paste(pop_name,'EUR_maf_bin_resume.txt',sep="_"))
print(t(EUR_maf_bin))
sink()

sink(paste(pop_name,'TGP_maf_bin_resume.txt',sep="_"))
print(t(TGP_maf_bin))
sink()

jpeg(paste(pop_name,"EUR_site_count_plot.jpg",sep="_"),width=1000, height=1000)
	barplot(t(EUR_maf_bin),col=c("green","darkgrey","red","white"))
dev.off()

jpeg(paste(pop_name,"TGP_site_count_plot.jpg",sep="_"),width=1000, height=1000)
	barplot(t(TGP_maf_bin),col=c("green","darkgrey","red","white"))
dev.off()


#19/02/2014
#work with exome chip data for genotype concordance check on rare variants
exome_tab <- read.table("/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/EXOME_CHIP/SEQ_SAMPLES/FVG_seq_sample_exome_freq.frq",header=T)
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
# maf_classes <- c(0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
maf_classes <- c(0,0.02,0.05,0.1,0.2,0.3,0.4,0.5)


#2014/03/21
# Check relatedness between samples on the first released set: most stringend VQSLOD threshold for FVG
#create a matrix from a 3 column dataframe
relatedness2_fvg <- read.csv('/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vqslod.fixed.noX.vcf.pop_info.relatedness2',sep="\t" ,header=T)
samples <- unique(relatedness2_fvg$INDV1)
m <- matrix(0, nrow=length(samples), ncol=length(samples))
rownames(m) <- colnames(m) <- samples
m[cbind(relatedness2_fvg$INDV1,relatedness2_fvg$INDV2)] <- relatedness2_fvg$RELATEDNESS_PHI
# m <- m + t(m)                                     # symmetrize
#m[upper.tri(m, diag=TRUE)] <- NA                  ## set redundant upper triangular matrix to NA
#write.table(m, "/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/doubletons_analysis/shared_doubletons_all_chroms_counts_matrix.txt", quote=F, sep="\t")   ## 3781 x 3781

related_1 <- relatedness2_fvg[which(relatedness2_fvg$RELATEDNESS_PHI > 0.177 & relatedness2_fvg$RELATEDNESS_PHI != 0.5),]

#this is used to extract uniqe pairs
unique(t(apply(related_1[,1:2],1,sort)))

#load sample file
samples_table <- read.table('/nfs/users/nfs_m/mc14/Work/SANGER/FVG/SEQ_CALLING/DOCS/samples_resume.txt',header=T,sep="\t")
#create village array for FVG
villages <- unique(samples_table$village)
villages <- tolower(as.character(villages))

#read kinship generated with KING software for the all sampleset for FVG from genotypes
king_kinship <- read.table('/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/MERGED/KINSHIP/king.kin0',head=T)

#read kinship generated for each village for FVG
library(pheatmap)
require(Hmisc)

for (i in villages) {

	current_village <- paste('/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/MERGED/AUTO/BY_VILLAGE/FVG_seq_subset_MERGED_sorted_right_sex_auto_',i,'.bed.kin0',sep="")
	gkin_name <- paste(i,"gkin",sep="_")
	sample_name <- paste(i,"samples",sep="_")
	matrix_name <- paste(i,"m",sep="_")

	#read the king output
	assign(gkin_name,read.table(current_village,header=T))
	current_gkin <- get(gkin_name)
	
	#create a relatedness matrix for the heatmap foreach village
	assign(sample_name,unique(c(current_gkin$ID1,current_gkin$ID2)))
	current_samples <- get(sample_name)

	current_matrix <- matrix(0, nrow=length(current_samples), ncol=length(current_samples))
	rownames(current_matrix) <- colnames(current_matrix) <- current_samples
	#we need to assign the correct level to each factor!
	# the problem here is that we have only the unique pairs, no duplicates, so this is different from the file generated from vcf tools
	# erto_m[cbind(factor(erto_gkin$ID1,levels=levels(as.factor(c(erto_gkin$ID1, erto_gkin$ID2)))),factor(erto_gkin$ID2,levels=levels(as.factor(c(erto_gkin$ID1, erto_gkin$ID2)))))] <- erto_gkin$Kinship
	current_matrix[cbind(factor(current_gkin$ID1,levels=sort(unique(c(current_gkin$ID1, current_gkin$ID2)))),factor(current_gkin$ID2,levels=sort(unique(c(current_gkin$ID1, current_gkin$ID2)))))] <- current_gkin$Kinship
	current_matrix <- current_matrix + t(current_matrix)
	diag(current_matrix) <- 0.5

	assign(matrix_name, current_matrix)


	jpeg(paste(i,"geno_relatedness_heatmap.jpg",sep="_"),width=1000, height=1000)
		pheatmap(current_matrix)
	dev.off()

	# now use the data from the WGS
	#In this case we don't need to do the factors trick becouse we have redundant informations in the result file
	current_W_village <- paste('/lustre/scratch113/projects/fvg_seq/REL-2014-01-09/v1/BY_VILLAGE/fvg.vqsr.beagle.impute2.anno.20140109.csq.pop.vqslod.fixed.noX.',capitalize(i),'.vcf.gz.pop_info.relatedness2',sep="")
	
	gkin_W_name <- paste(capitalize(i),"W_gkin",sep="_")
	sample_W_name <- paste(capitalize(i),"W_samples",sep="_")
	matrix_W_name <- paste(capitalize(i),"W_m",sep="_")

	assign(gkin_W_name,read.table(current_W_village,header=T))
	current_W_gkin <- get(gkin_W_name)
	#set the same names as the genotypes
	# current_W_merged1 <- merge(current_W_gkin,samples_table,by.x="INDV2",by.y="called_id")
	# current_W_merged2 <- merge(current_W_merged1,samples_table,by.x="INDV1",by.y="called_id")
	# current_W_gkin <- current_W_merged2
	# current_W_gkin$clinic_id.x <- as.factor(current_W_merged2$clinic_id.x)
	# current_W_gkin$clinic_id.y <- as.factor(current_W_merged2$clinic_id.y)
	# current_W_gkin <- t(apply(current_W_gkin[,1:2],1,sort))

	#create a relatedness matrix for the heatmap foreach village
	assign(sample_W_name,unique(current_W_gkin$INDV1))
	current_W_samples <- get(sample_W_name)

	current_W_matrix <- matrix(0, nrow=length(current_W_samples), ncol=length(current_W_samples))
	rownames(current_W_matrix) <- colnames(current_W_matrix) <- current_W_samples
	current_W_matrix[cbind(current_W_gkin$INDV1,current_W_gkin$INDV2)] <- current_W_gkin$RELATEDNESS_PHI
	# current_W_matrix[cbind(current_W_gkin$clinic_id.x,current_W_gkin$clinic_id.y)] <- current_W_gkin$RELATEDNESS_PHI
	assign(matrix_W_name, current_W_matrix)

	jpeg(paste(i,"W_relatedness_heatmap.jpg",sep="_"),width=1000, height=1000)
    	pheatmap(current_W_matrix)
  	dev.off()

}

#upload unrelated list of samples based on genotypes:
unrelated_3 <- read.table('/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/MERGED/AUTO/kingunrelated.txt',header=F)
colnames(unrelated_3) <- c("FID","IID")

#extract sequence id for unrelated samples
unrelated_sequenced <- samples_table[which(samples_table$clinic_id %in% unrelated_3$IID),]

# > summary(unrelated_sequenced)
#    called_id     clinic_id      seq_pipeline_id     sex         village  
#  591041 :  1   Min.   :582105   591041 :  1     Female:102   Erto   :40  
#  591669 :  1   1st Qu.:591345   591669 :  1     Male  : 59   Illegio:31  
#  591675 :  1   Median :593575   591675 :  1                  Resia  :52  
#  591690 :  1   Mean   :593098   591690 :  1                  Sauris :38  
#  591695 :  1   3rd Qu.:594962   591695 :  1                              
#  593569 :  1   Max.   :605885   593569 :  1                              

#write table of unrelated people with corresponding sequence id
write.table(unrelated_sequenced, file='/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/MERGED/AUTO/FVG_unrelated_sequenced.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = FALSE)

