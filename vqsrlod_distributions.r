#script to perform some qc on VBI seq data
rm(list=ls())

#read commandline args
args <- commandArgs(trailing=TRUE)

#Extract SNPs common between EUR and VBI (the same process would be applied to TSI)
#Extraction is performed using chr and position match, because some variants don't have an rsID
#
#Set arguments to use from command line
#vbi_file_path <- args[[1]]
#alternative_file_path <- args[[2]]
#i <- args[[3]]

#manual set arguments...manual test on VBI.chr20.snps.vcf
vbi_file_path <- '/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/VBI/VBI.23chr.vqslod.snps.vcf'

#VBI_snps <- read.table(vbi_file_path, sep='\t',header=T,comment.char="",colClasses=c("character","numeric","character","numeric"))
#FIXME:this is fixed for files with header extracted from Klaudia perl script...to have a file with standardized header could be better
#standardize header
#colnames(VBI_snps) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","DP","AC1","AF1","VQSLOD","IMP2","HWE","AC","AN","AF_EUR","AF_MAX")
#FIXME: this could be done using perl/bash
#remove useless columns

#now plot the density distribution for VQSRLOD
#jpeg("VQSRLOD_overall_density.jpg", width=1000, height=1000)
#plot(density(VBI_snps$VQSLOD),main="VQSLOD distribution",xlab="VQSLOD",ylab="VQSLOD density distribution")
#dev.off()
#jpeg("VQSRLOD_overall_hist.jpg", width=1000, height=1000)
#hist(VBI_snps$VQSLOD,main="VQSLOD distribution",xlab="VQSLOD",ylab="VQSLOD density distribution")
#dev.off()
#rm(VBI_snps)
#gc()
#then we need to do it on multiallelic snps
#VBI_multi_basepath <- "/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/QC_OUT/MULTIALLELIC/"
vbi_multi_basepath <- "/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/QC_OUT/MULTIALLELIC/VQSLOD/multiallelic_table.vqslod.csv"
#VBI_multi <- NULL
#for (i in 1:22) {
#	print(i)
#	current_multi_file <- paste(VBI_multi_basepath,"multiallelic_chr",i,"_table.vqslod.csv",sep="")
	current_multi_snps <- read.table(vbi_multi_basepath, sep='\t',header=F,comment.char="", colClasses=c("character","numeric","character","numeric"))
	colnames(current_multi_snps) <- c("CHROM","POS","ID","VQSLOD")
	VBI_multi <- current_multi_snps
#	rm(current_multi_snps)
#	gc()
#}

#plot VQSLOD for multiallelic SNPs
jpeg("VQSRLOD_multiallelic_density.jpg", width=1000, height=1000)
plot(density(VBI_multi$VQSLOD),main="VQSLOD distribution",xlab="VQSLOD",ylab="VQSLOD density distribution")
dev.off()
jpeg("VQSRLOD_multiallelic_hist.jpg", width=1000, height=1000)
hist(VBI_multi$VQSLOD,main="VQSLOD distribution",xlab="VQSLOD",ylab="VQSLOD count")
dev.off()
rm(VBI_multi)
gc()
#then we need to do it on biallelic snps
VBI_bi_basepath <- "/nfs/users/nfs_m/mc14/lustre_home/GENOTIPI/COMPARISON/VBSEQ_QC/QC_OUT/SNPS_FILES/VQSLOD/VBI_table.vqslod.csv"
VBI_bi <- NULL
#for (i in 1:22) {
#	print(i)
#	current_bi_file <- paste(VBI_bi_basepath,"VBI_chr",i,".snps.vqslod.csv",sep="")
#	print("path read!")
#	current_bi_snps <- read.table(vbi_file_path, sep='\t',header=T,row.names = NULL,comment.char="",colClasses=c("numeric","numeric","character","character","character","numeric","character","numeric","numeric","numeric","numeric","character","numeric","character","numeric","numeric","numeric"))
#	current_bi_snps <- read.table(vbi_file_path, sep='\t',header=T,row.names = NULL,comment.char="",colClasses=c("character","numeric","character","numeric"))
#	print("file read!")
	current_bi_snps <- read.table(VBI_bi_basepath, sep='\t',header=F,comment.char="", colClasses=c("character","numeric","character","numeric"))
        colnames(current_bi_snps) <- c("CHROM","POS","ID","VQSLOD")
	VBI_bi <- current_bi_snps
#	print("rbind done!")
#	rm(current_bi_snps)
#	gc()
#}

#plot VQSLOD for multiallelic SNPs
jpeg("VQSRLOD_biallelic_density.jpg", width=1000, height=1000)
plot(density(VBI_bi$VQSLOD),main="VQSLOD distribution",xlab="VQSLOD",ylab="VQSLOD density distribution")
dev.off()
jpeg("VQSRLOD_biallelic_hist.jpg", width=1000, height=1000)
hist(VBI_bi$VQSLOD,main="VQSLOD distribution",xlab="VQSLOD",ylab="VQSLOD count")
dev.off()
rm(VBI_bi)
gc()
q(save="no")
