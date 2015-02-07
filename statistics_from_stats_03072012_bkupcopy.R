#script to estract statisics fron vcf-stats results
#read commandline args
args <- commandArgs(trailing=TRUE)

#Args list:
#args[[1]]:list with individuals and sex (absolute path) (i.e '/nfs/users/nfs_m/mc14/Work/SANGER/DOCS/ValBorbera/sequenced_id_conversion_sex.txt')
#args[[2]]:the path of vcf-stats files to process (absolute path)
#args[[3]]: prefix used to assign names to vcf-stats output files

#Assign args
indiv_file <- args[[1]]
vcfstats_filepath <- args[[2]]
vcfstats_file_prefix <- args[[3]]

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

#INDELS
indels <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'indels',sep="."),sep="/"), header=F, skip=2,col.names=c('count','sample'))
indels$sample <- gsub("samples/","",indels$sample)

indels_summary <- summary(indels)

write.table(indels_summary, file='indels_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("indels_density.jpg",width=1000, height=1000)
plot(density(indels$count),main ="Indels density distribution in VBSEQ")
dev.off()

jpeg("indels_boxplot.jpg", width=1000, height=1000)
indels_boxplot <- boxplot(indels$count,main ="Indels in VBSEQ", xlab="All individuals", ylab="Indels count")
dev.off()

#male/female boxplot
jpeg("indels_boxplot_male_female.jpg", width=1000, height=1000)
indels_boxplot_mf <- boxplot(indels[which(indels$sample %in% male_inds$SEQ_ID),]$count,indels[which(indels$sample %in% female_inds$SEQ_ID),]$count,main ="Indels in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Indels count")
dev.off()

#PRIVATE SNPS
private_s <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'private',sep="."),sep="/"), header=F, skip=1,col.names=c('count','sample'))

private_s_summary <- summary(private_s)

write.table(private_s_summary, file='private_snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("private_snps_density.jpg", width=1000, height=1000)
plot(density(private_s$count),main ="Private SNPs density distribution in VBSEQ")
dev.off()

jpeg("private_snps_boxplot.jpg", width=1000, height=1000)
private_s_boxplot <- boxplot(private_s$count,main ="Private SNPs in VBSEQ", xlab="All individuals", ylab="Private SNPs count")
dev.off()

#male/female boxplot
jpeg("private_snps_boxplot_male_female.jpg", width=1000, height=1000)
private_s_boxplot_mf <- boxplot(private_s[which(private_s$sample %in% male_inds$SEQ_ID),]$count,private_s[which(private_s$sample %in% female_inds$SEQ_ID),]$count,main ="Private SNPs in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Private SNPs count")
dev.off()

#SHARED SNPS
#Number of sites having a non-reference allele in 0,1,2,etc samples

shared_s <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'shared',sep="."),sep="/"), header=F, skip=1,col.names=c('ind_count','freq'))

shared_s_summary <- summary(shared_s)

write.table(shared_s_summary, file='shared_snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("shared_snps_density.jpg", width=1000, height=1000)
plot(density(shared_s$freq),main ="Density distribution of the number of sites having a non-reference allele in VBSEQ")
dev.off()

jpeg("shared_snps_hist.jpg", width=1000, height=1000)
barplot(shared_s$freq,names.arg=shared_s$ind_count,main ="Number of sites having a non-reference allele in VBSEQ", xlab="Individuals", ylab="Number of sites having a non-reference allel (count)")
dev.off()

#SNPS
#Number of positions with SNPs

snps <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'snps',sep="."),sep="/"), header=F, skip=2,col.names=c('count','sample'))
snps$sample <- gsub("samples/","",snps$sample)

snps_summary <- summary(snps)

write.table(snps_summary, file='snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("snps_density.jpg", width=1000, height=1000)
plot(density(snps$count),main ="SNPs density distribution in VBSEQ")
dev.off()

jpeg("snps_boxplot.jpg", width=1000, height=1000)
snps_boxplot <- boxplot(snps$count,main ="SNPs in VBSEQ", xlab="All individuals", ylab="SNPs count")
dev.off()

#male/female boxplot
jpeg("snps_boxplot_male_female.jpg", width=1000, height=1000)
snps_boxplot_mf <- boxplot(snps[which(snps$sample %in% male_inds$SEQ_ID),]$count,snps[which(snps$sample %in% female_inds$SEQ_ID),]$count,main ="SNPs in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="SNPs count")
dev.off()


#SAMPLE TS/TV
#samples-tstv

s_tstv <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'samples-tstv',sep="."),sep="/"), header=F,skip=1,col.names=c('Transitions','Transversions','ts_tv_ratio','sample'))

s_tstv_summary <- summary(s_tstv)

write.table(s_tstv_summary, file='s_tstv_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("s_ts_density.jpg", width=1000, height=1000)
plot(density(s_tstv$Transitions),col='red',main ="Transitions density distribution in VBSEQ")
dev.off()

jpeg("s_tv_density.jpg", width=1000, height=1000)
plot(density(s_tstv$Transversions),col='green',main ="Transversions density distribution in VBSEQ")
dev.off()

jpeg("s_tstv_density.jpg", width=1000, height=1000)
plot(density(s_tstv$ts_tv_ratio),col='blue',main ="Ts/Tv ratio density distribution in VBSEQ")
dev.off()

jpeg("s_ts_boxplot.jpg", width=1000, height=1000)
s_ts_boxplot <- boxplot(s_tstv$Transitions,main ="Transition in VBSEQ", xlab="All individuals", names=c("Transitions"))
dev.off()

jpeg("s_tv_boxplot.jpg", width=1000, height=1000)
s_tv_boxplot <- boxplot(s_tstv$Transversions,main ="Transversion in VBSEQ", xlab="All individuals", names=c("Transversions"))
dev.off()

jpeg("s_tstv_boxplot.jpg", width=1000, height=1000)
s_tstv_boxplot <- boxplot(s_tstv$ts_tv_ratio,main ="Ts/Tv ratio in VBSEQ", xlab="All individuals", names=c("Ts/Tv"))
dev.off()
#to retrieve the outliers:
outliers <- s_tstv[which(s_tstv$Transversions %in%  s_tv_boxplot$out),]$sample

#male/female boxplot
jpeg("s_ts_boxplot_male_female.jpg", width=1000, height=1000)
s_ts_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$Transitions,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$Transitions,main ="Transitions in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="SNPs count")
dev.off()

jpeg("s_tv_boxplot_male_female.jpg", width=1000, height=1000)
s_tv_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$Transversions,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$Transversions,main ="Transversions in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="SNPs count")
dev.off()

jpeg("s_tstv_boxplot_male_female.jpg", width=1000, height=1000)
s_tstv_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$ts_tv_ratio,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$ts_tv_ratio,main ="Ts/Tv ratio in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="SNPs count")
dev.off()

