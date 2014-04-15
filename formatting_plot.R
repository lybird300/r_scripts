args=commandArgs(trailing=TRUE)
#args[[1]] = trait
#args[[2]] = chr
#args[[3]] = result_file
#args[[4]] = maf_file


# MANHATTAN PLOT
source("/nfs/users/nfs_m/mc14/Work/r_scripts/qqman.R")

# PLOT ALL SNPS
# REDUCING # SNPs TO PLOT 
#for(i in 1:22){
#system(paste("awk -F, 'BEGIN{FS=\" \";OFS=\" \"}{if ($10 ==",i," && $16 < 0.01) print $0 }' *new.palinear | cut -f 1,10,11,16 -d ' ' >> to_plot_unfiltered", sep=""))
#system(paste("awk -F, 'BEGIN{FS=\" \";OFS=\" \"}{if ($10 ==",i," && $16 >= 0.01) print $0 }' *new.palinear | sort -t ' ' -k2 | head -n 150000 | cut -f 1,10,11,16 -d ' ' >> to_plot_unfiltered", sep=""))
#}

to_plot_unfiltered = read.csv(args[[3]], sep="\t", header=T)
# to_plot_unfiltered = read.csv("chr6.txt.assoc.logistic.mod", sep="\t", header=T)

to_plot_mafs = read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/GWA/gwa_data_freq.frq.mod",header=T,sep="\t")

# to add maf info
merged_new <- merge(to_plot_unfiltered,to_plot_mafs,by.x="SNP",by.y="SNP",all.x)
merged_new$SNP <- as.character(merged_new$SNP)
merged_new$CHR.y <- NULL
merged_new$A1.y <- NULL
names(merged_new)[2] <- "CHR"
names(merged_new)[4] <- "A1"

dim(to_plot_unfiltered)
dim(merged_new)

#colnames(to_plot_unfiltered)=c("SNP","CHR","BP","P")
#colnames(to_plot_unfiltered)=c("CHR","SNP","BP","P")
#colnames(to_plot_unfiltered)=c("CHR","BP","P")

dim(to_plot_unfiltered[which(to_plot_unfiltered$P < 1 & to_plot_unfiltered$P >= 0.1),])
dim(to_plot_unfiltered[which(to_plot_unfiltered$P < 0.1 & to_plot_unfiltered$P >= 0.01),])
dim(to_plot_unfiltered[which(to_plot_unfiltered$P < 0.01),])

#extract snps to color
low_to_comm <- merged_new[which(merged_new$MAF >= 0.01 & merged_new$MAF < 0.05),]$SNP
#remove snps with pval >= 0.2
low_to_comm <- to_plot_unfiltered[-which(to_plot_unfiltered$P < 1 & to_plot_unfiltered$P >= 0.2),])

jpeg(paste(args[[1]], "_", args[[2]], ".manhattan.unfiltered.jpg", sep=""), width = 750, height = 300, pointsize = 10)
manhattan(merged_new, annotate=low_to_comm, pch=16, main=paste(args[[1]], "_", args[[2]], " unfiltered", sep=""))
# jpeg(paste("VLDL", "_", "all", ".manhattan.unfiltered.jpg", sep=""), width = 1300, height = 600, pointsize = 16)
# manhattan(to_plot_unfiltered, pch=16, main=paste("VLDL", "_", "all", " unfiltered", sep=""))
dev.off()

jpeg(paste(args[[1]], "_", args[[2]], ".qq.unfiltered.jpg", sep=""), width = 750, height = 300, pointsize = 10)
qq(to_plot_unfiltered$P, main=paste(args[[1]], "_", args[[2]], " unfiltered", sep=""))
# jpeg(paste("VLDL", "_", "all", ".qq.unfiltered.jpg", sep=""), width = 1300, height = 600, pointsize = 16)
# qq(to_plot_unfiltered$P, pch=16, main=paste("VLDL", "_", "all", " unfiltered", sep=""))
dev.off()

# PLOT ONLY SNPS WITH RSQ > 0.4
#for(i in 1:22){
#system(paste("awk -F, 'BEGIN{FS=\" \";OFS=\" \"}{if ($10 ==",i," && $16 < 0.01 && $7 >= 0.4) print $0 }' *new.palinear | cut -f 1,10,11,16 -d ' ' >> to_plot_filtered", sep=""))
#system(paste("awk -F, 'BEGIN{FS=\" \";OFS=\" \"}{if ($10 ==",i," && $16 >= 0.01 && $7 >= 0.4) print $0 }' *new.palinear | sort -t ' ' -k2 | head -n 150000 | cut -f 1,10,11,16 -d ' ' >> to_plot_filtered", sep=""))
#}

#to_plot_filtered = read.csv("to_plot_filtered", sep=" ", header=F)
#dim(to_plot_filtered)
#colnames(to_plot_filtered)=c("SNP","CHR","BP","P")

#dim(to_plot_filtered[which(to_plot_filtered$P < 1 & to_plot_filtered$P >= 0.1),])
#dim(to_plot_filtered[which(to_plot_filtered$P < 0.1 & to_plot_filtered$P >= 0.01),])
#dim(to_plot_filtered[which(to_plot_filtered$P < 0.01),])

#jpeg(paste(args[[1]], "_", args[[2]], ".palinear.manhattan.filtered.jpg", sep=""), width = 1300, height = 600, pointsize = 16)
#manhattan(to_plot_filtered, pch=16, main=paste(args[[1]], "_", args[[2]], " filtered", sep=""))
#dev.off()

###########################################################
#require(GWAStoolbox)
#require(GWAtoolbox)
# build script

#write("# Description of input data columns", file="to_GWAtoolbox")
#write("MARKER name", append=T, file="to_GWAtoolbox")
#write("CHR chrom", append=T, file="to_GWAtoolbox")
#write("POSITION position", append=T, file="to_GWAtoolbox")
#write("N n", append=T, file="to_GWAtoolbox")
#write("ALLELE A1 A2", append=T, file="to_GWAtoolbox")
#write("EFFECT beta_SNP_add", append=T, file="to_GWAtoolbox")
#write("STDERR sebeta_SNP_add", append=T, file="to_GWAtoolbox")
#write("PVALUE pgc", append=T, file="to_GWAtoolbox")
#write("FREQLABEL Mean_predictor_allele", append=T, file="to_GWAtoolbox")
#write("IMP_QUALITY Rsq", append=T, file="to_GWAtoolbox")
#
#write("# High quality filters", append=T, file="to_GWAtoolbox")
#write("HQ_SNP 0.01 0.3", append=T, file="to_GWAtoolbox")
#
#write("# Plotting filters", append=T, file="to_GWAtoolbox")
#write("MAF 0.01 0.05", append=T, file="to_GWAtoolbox")
#write("IMP 0.4 0.5", append=T, file="to_GWAtoolbox")
#
#write("# Prefix for output files", append=T, file="to_GWAtoolbox")
#write("PREFIX Result_", append=T, file="to_GWAtoolbox")
#
#write("# Input file with GWA data", append=T, file="to_GWAtoolbox")
#write(paste("PROCESS ", args[[1]], "_", args[[2]], ".palinear", sep=""), append=T, file="to_GWAtoolbox")
#write("SEPARATOR WHITESPACE", append=T, file="to_GWAtoolbox")
#
#gwasqc("to_GWAtoolbox")
##################################################################

q(save="no")


