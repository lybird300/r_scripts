#script for statistics on maf...
rm(list=ls())
require(gplots)

#Set arguments to use from command line
# 
#read commandline args
args <- commandArgs(trailing=TRUE)

filepath <- args[[1]]
# filepath <- "/lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140410_ANNOTATE/MAF/ALL.maf_table.snp.tab"
pop_name <- args[[2]]
# pop_name <- "FVG"
#chr <- args[[2]]

#Upload data:
#all_maf <- read.table(filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="",colClasses=c("character","character","numeric","numeric","numeric","numeric","numeric","numeric"))
#all_maf <- read.table(filepath,header=F,sep=' ', stringsAsFactors=F, comment.char="")
# all_maf <- read.table(filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="")
all_maf <- read.table(filepath,header=T,sep=' ', stringsAsFactors=F, comment.char="")
#FIXME: the best would be to have an headed file...
#colnames(all_maf) <- c('CHROM','POS','ID','AN','AC','DP','AF','MAF')
# colnames(all_maf) <- c('CHROM','POS','ID','AC','AN','AF','MAF')
# colnames(all_maf) <- c('CHROM','POS','REF','ALT','MAF')
#FIXME:different header for different format...maybe we could set up so sort of recognition for columns we want to use...
##################################################
#BEST FIX: provide the header as a parameter!!!###
##################################################
# colnames(all_maf) <- c('rsid','pos','allele_A','allele_B','index','average_maximum_posterior_call','info','cohort_1_AA','cohort_1_AB','all_BB','all_NULL','MAF','missing_data_proportion','cohort_1_hwe','all_A_freq','all_B_freq','minor_ALL')
#colnames(all_maf) <- c('CHROM','ID','POS','exp_freq_a1','info','certainty','type','info_type1','concord_type1','r2_type1','info_type0','concord_type0','r2_type0','MAF')
#colnames(all_maf) <- c('CHROM','ID','POS','exp_freq_a1','info','certainty','type','MAF')

#if you want to extract only SNPs after uploading the table:
all_maf_snps <- all_maf[which(nchar(all_maf$ALT) == nchar(all_maf$REF)),]

#first write a summary (for snps only)
# write.table(summary(all_maf),file="complete_summary.txt",sep="\t",col.names=T,quote=F,row.names=F)
write.table(summary(all_maf_snps),file="complete_snp_summary.txt",sep="\t",col.names=T,quote=F,row.names=F)

#now plot maf density
jpeg(paste(pop_name,"maf_density_plot.jpg",sep="_"),width=800, height=800)
par(lab=c(10,10,12),cex=2)
  plot(density(all_maf_snps$MAF), main="Maf density distribution", xlab="Maf",ylab="Maf density")
dev.off()

#count monomorphic sites
all_mono_snps <- all_maf_snps[which(all_maf_snps$MAF == 0),]
dim(all_mono_snps)
#count recursively until AC = 5
all_count_snps <- NULL

for (i in 1:50) {
  all_current_ac <- all_maf_snps[which(all_maf_snps$AC == i | all_maf_snps$AC == (all_maf_snps$AN - i)),]
  # print(paste("Current allele count:",i,sep=""))
  # print(dim(all_current_ac))
  # assign(paste("all_ac_",i,sep=""),all_current_ac)
  all_count_current <- c(i,length(all_current_ac$CHROM))
  all_count_snps <- rbind(all_count_snps,all_count_current)
}

#set the class to data.frame
all_count_snps <- as.data.frame(all_count_snps)
rownames(all_count_snps) <- NULL
colnames(all_count_snps) <- c("AC","count")

#plot those numbers
jpeg(paste(pop_name,"site_allele_count_plot.jpg",sep="_"),width=800, height=800)
  barplot(all_count_snps$count,names.arg=all_count_snps$AC,col=colors()[72])
dev.off()

#now split all maf in frequencies bins
#this is a function that splits in bins of maf and writes some summaries
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
maf_classes <- c(0,0.005,0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)
# maf_classes <- c(0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
# maf_classes <- c(0,0.02,0.05,0.1,0.2,0.3,0.4,0.5)

all_maf_classes <- split_bins(maf_classes,all_maf_snps[which(all_maf_snps$MAF > 0),],pop_name)
# all_maf_classes <- split_bins(maf_classes,maf_table,pop_name)
gc()
#write a cute output
sink(paste(pop_name,'maf_bin_resume.txt',sep="_"))
print(all_maf_classes)
sink()

jpeg(paste(pop_name,"site_count_plot.jpg",sep="_"),width=800, height=800)
	barplot(as.matrix(all_maf_classes),names.arg=maf_classes[2:length(maf_classes)],col=colors()[72])
dev.off()

q(save='no')
