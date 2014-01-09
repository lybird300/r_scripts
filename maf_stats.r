#script for statistics on maf...
rm(list=ls())
require(gplots)

#Set arguments to use from command line
# 
#read commandline args
args <- commandArgs(trailing=TRUE)

filepath <- args[[1]]
pop_name <- args[[2]]
#chr <- args[[2]]

#Upload data:
#all_maf <- read.table(filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="",colClasses=c("character","character","numeric","numeric","numeric","numeric","numeric","numeric"))
#all_maf <- read.table(filepath,header=F,sep=' ', stringsAsFactors=F, comment.char="")
# all_maf <- read.table(filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="")
all_maf <- read.table(filepath,header=T,sep=' ', stringsAsFactors=F, comment.char="")
#FIXME: the best would be to have an headed file...
#colnames(all_maf) <- c('CHROM','POS','ID','AN','AC','DP','AF','MAF')
# colnames(all_maf) <- c('CHROM','POS','ID','AC','AN','AF','MAF')
#FIXME:different header for different format...maybe we could set up so sort of recognition for columns we want to use...
##################################################
#BEST FIX: provide the header as a parameter!!!###
##################################################
colnames(all_maf) <- c('rsid','pos','allele_A','allele_B','index','average_maximum_posterior_call','info','cohort_1_AA','cohort_1_AB','all_BB','all_NULL','MAF','missing_data_proportion','cohort_1_hwe','all_A_freq','all_B_freq','minor_ALL')
#colnames(all_maf) <- c('CHROM','ID','POS','exp_freq_a1','info','certainty','type','info_type1','concord_type1','r2_type1','info_type0','concord_type0','r2_type0','MAF')
#colnames(all_maf) <- c('CHROM','ID','POS','exp_freq_a1','info','certainty','type','MAF')

#first write a summary
write.table(summary(all_maf),file="complete_summary.txt",sep="\t",col.names=T,quote=F,row.names=F)

#now plot maf density
jpeg("maf_density_plot.jpg",width=1000, height=1000)
par(lab=c(10,10,12),cex=2)
  plot(density(all_maf$MAF), main="Maf density distribution", xlab="Maf",ylab="Maf density")
dev.off()

#now split all maf in frequencies bins
#this is a function that splits in bins of maf and writes some summaries
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
maf_classes <- c(0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
# maf_classes <- c(0,0.02,0.05,0.1,0.2,0.3,0.4,0.5)

all_maf_classes <- split_bins(maf_classes,all_maf,pop_name)
gc()
#write a cute output
sink('maf_bin_resume.txt')
print(all_maf_classes)
sink()

q(save='no')
