#Commands for statistics on imputation with 1000GP + VBI

info <- read.table("CHR7/chr7.geno_info",header=F , sep=" ")
colnames(info) <- c("snp_id","rs_id","position","exp_freq_a1","info","certainty","type","info_type0","concord_type0","r2_type0")

info$snp_id<-NULL

summary(info)

#we want a summary for the concordance and a summary for the af
# exp_freq_a1    
#Min.    0.00000
#1st Qu. 0.00000
#Median  0.00100
#Mean    0.08243
#3rd Qu. 0.03600
#Max.    1.00000

#     info     
#Min.    0.0000
#1stQu. 0.0100
#Median  0.3640
#Mean    0.4276
#3rdQu. 0.8520
#Max.    1.0000

#keep only genotyped sites
info_purged_conc <- info[which(info$concord_type0 != -1),]

#summary for genotype data in imputation files
# exp_freq_a1  
#Min.    0.0100
#1st Qu. 0.1750
#Median  0.3440
#Mean    0.3897
#3rd Qu. 0.5760
#Max.    0.9870

#  info_type0     concord_type0       r2_type0    
#Min.   :0.0000   Min.    0.0620   Min.    -1.0000
#1st Qu.:0.9120   1stQu. 0.9450   1stQu.  0.8570
#Median :0.9640   Median  0.9800   Median   0.9490
#Mean   :0.9307   Mean    0.9546   Mean     0.8856
#3rd Qu.:0.9880   3rdQu. 0.9950   3rdQu.  0.9850
#Max.   :1.0000   Max.    1.0000   Max.     1.0000



plot_name <- "DENSITY_INFO.jpg"
jpeg(plot_name, width=1000, height=1000)
  par(lab=c(15,15,17),cex=2)
  plot(density(info$info), main="Distribution of INFO score",xlab="INFO score",ylab="INFO score density")
dev.off()

plot_name <- "DENSITY_r2.jpg"
jpeg(plot_name, width=1000, height=1000)
  par(lab=c(15,15,17),cex=2)
  plot(density(info_purged_conc$r2_type0), main="Distribution of R2 ",xlab="R2",ylab="R2 density")
dev.off()
plot_name <- "DENSITY_concordance.jpg"
jpeg(plot_name, width=1000, height=1000)
  par(lab=c(15,15,17),cex=2)
  plot(density(info_purged_conc$concord_type0), main="Distribution of Genotype concordance ",xlab="Genotype concordance",ylab="Genotype concordance density")
dev.off()

# info_maf <- read.table("CHR7/chr7.geno_info_maf",header=T , sep=" ")
info_vbuk_maf <- read.table("CHR7/chr7.geno_info_maf",header=T , sep=" ")
info_uk1kg_maf <- read.table("CHR7/chr7.geno_info_maf",header=T , sep=" ")
# info_maf_purged_conc <- info[which(info_maf$concord_type0 != -1),]
# info_vbuk_maf_purged_conc <- info_vbuk_maf[which(info_vbuk_maf$concord_type0 != -1),]
info_uk1kg_maf_purged_conc <- info_uk1kg_maf[which(info_uk1kg_maf$concord_type0 != -1),]
source('~/Work/r_scripts/maf_bins_splitter.r')
#maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
# maf_classes <- c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
maf_classes <- c(0,0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)

# split_bins(maf_classes,info_maf,"VBI_1kg_5")
# split_bins(maf_classes,info_maf_purged_conc,"VBI1kgp_panel")
split_bins(maf_classes,info_uk1kg_maf_purged_conc,"UK10K1KP_panel")


#Script and function to plot rsquared correlation splitted in frequency bins
# basepath <- "/lustre/scratch110/sanger/mc14/GENOTIPI/TEST_IMPUTATION"
# output3_info <- read.table(paste(basepath,"OUPUT3/MERGED/SNPS/MACH_FORMAT/0/16.machinfo.mod.new",sep="/"),sep="\t",header=T)
# summary(output3_info)
# output3_rsqr <- read.table(paste(basepath,"OUPUT3/MERGED/SNPS/MACH_FORMAT/0/RSQR/chr16_vbi.doseR2",sep="/"),header=T,sep=" ")
# summary(output3_rsqr)
# output3_merged <- merge(output3_rsqr,output3_info,by.x="SNP",by.y="name",all=FALSE,sort=FALSE)
# summary(output3_merged)

# output3_info <- read.table(paste(basepath,"OUPUT3/MERGED/RSQUARED_NEW/chr16.complete_R2_uniq",sep="/"),sep=" ",header=T)
# output4_info <- read.table(paste(basepath,"OUPUT4/MERGED/RSQUARED_NEW/chr16.complete_R2_uniq",sep="/"),sep=" ",header=T)
# output5_info <- read.table(paste(basepath,"OUPUT5/MERGED/RSQUARED_NEW/chr16.complete_R2_uniq",sep="/"),sep=" ",header=T)

# #define maf classes
# maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
# maf_classes <- c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)

# source('~/Work/r_scripts/maf_bins_splitter.r')
# split_bins(maf_classes,output3_info,"OUT3")
# split_bins(maf_classes,output4_info,"OUT4")
# split_bins(maf_classes,output5_info,"OUT5")

# jpeg("scatter_r2.jpg", width=1000, height=1000)
# plot(output3_info$R2,lwd=2,pch=20,col="red")
# points(output4_info$R2,lwd=2,pch=21,col="green")
# points(output5_info$R2,lwd=2,pch=22,col="deepskyblue3")
# dev.off()


#some addition to plot other stuff
# bin_average <- NULL
# bin_median <- NULL
# for(i in 2:length(maf_classes)){
#   bin_name <- paste("bin",maf_classes[i],sep="_")

#   current_bin <- read.table(paste(basepath,"/SUMMARIES/IMP_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  
#   assign(bin_name, current_bin)
#   bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
#   bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
  
# }

# ##########################

# basepath <- "/lustre/scratch110/sanger/mc14/GENOTIPI/SEQ_IMP_RSQR/RSQUARED"
# chr7_info2 <- read.table(paste(basepath,"chr7.complete_R2_2_uniq",sep="/"),sep=" ",header=T)
# maf_classes <- c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
# split_bins(maf_classes,chr7_info2,"IMP2")

bin_1_r2_average <- NULL
bin_2_r2_average <- NULL
bin_3_r2_average <- NULL
bin_4_r2_average <- NULL
# bin_r2_average <- NULL
# bin_r2_median <- NULL
# bin_r2_imp_average <- NULL
# bin_r2_imp_median <- NULL

for(i in 2:length(maf_classes)){
  # bin_1_name <- paste("bin_1",maf_classes[i],sep="_")
  # bin_2_name <- paste("bin_2",maf_classes[i],sep="_")
  #bin_3_name <- paste("bin_3",maf_classes[i],sep="_")
  bin_4_name <- paste("bin_4",maf_classes[i],sep="_")

  # current_1_bin <- read.table(paste("1KGP_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_1_bin <- current_1_bin[which(current_1_bin$concord_type0 != -1),]
  # current_1_bin$r2_type0 <- as.numeric(as.character(current_1_bin$r2_type0))
  # current_2_bin <- read.table(paste("VBI1kgp_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_2_bin <- current_2_bin[which(current_2_bin$concord_type0 != -1),]
  # current_2_bin$r2_type0 <- as.numeric(as.character(current_2_bin$r2_type0))
  # current_3_bin <- read.table(paste("VBIUK10K_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_3_bin <- current_3_bin[which(current_3_bin$concord_type0 != -1),]
  # current_3_bin$r2_type0 <- as.numeric(as.character(current_3_bin$r2_type0))
  current_4_bin <- read.table(paste("UK10K1KP_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  current_4_bin <- current_4_bin[which(current_4_bin$concord_type0 != -1),]
  current_4_bin$r2_type0 <- as.numeric(as.character(current_4_bin$r2_type0))
  
  # assign(bin_1_name, current_1_bin)
  # assign(bin_2_name, current_2_bin)
  # assign(bin_3_name, current_3_bin)
  assign(bin_4_name, current_4_bin)
  #bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
  #bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
  #bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
  #bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
  # bin_1_r2_average <- c(bin_1_r2_average,mean(current_1_bin$r2_type0,na.rm=TRUE))
  # bin_2_r2_average <- c(bin_2_r2_average,mean(current_2_bin$r2_type0,na.rm=TRUE))
  # bin_3_r2_average <- c(bin_3_r2_average,mean(current_3_bin$r2_type0,na.rm=TRUE))
  bin_4_r2_average <- c(bin_4_r2_average,mean(current_4_bin$r2_type0,na.rm=TRUE))
  # bin_r2_median <- c(bin_r2_median,median(current_bin$r2_type0,na.rm=TRUE))
  # bin_r2_imp_average <- c(bin_r2_imp_average,mean(current_bin$IMP_R2,na.rm=TRUE))
  # bin_r2_imp_median <- c(bin_r2_imp_median,median(current_bin$IMP_R2,na.rm=TRUE))
  
}


jpeg("boxplot_r2.jpg", width=800, height=800)
boxplot(bin_0.01$r2_type0,bin_0.02$r2_type0,bin_0.03$r2_type0,bin_0.04$r2_type0,bin_0.05$r2_type0,bin_0.1$r2_type0,bin_0.2$r2_type0,bin_0.3$r2_type0,bin_0.4$r2_type0,bin_0.5$r2_type0,names=maf_classes[2:11],xlab="Freq bin", ylab="R2",main="R2 score",varwidth=TRUE)
dev.off()

# jpeg("boxplot_imp_r2.jpg", width=800, height=800)
# boxplot(bin_0.01$IMP_R2,bin_0.02$IMP_R2,bin_0.03$IMP_R2,bin_0.04$IMP_R2,bin_0.05$IMP_R2,bin_0.1$IMP_R2,bin_0.2$IMP_R2,bin_0.3$IMP_R2,bin_0.4$IMP_R2,bin_0.5$IMP_R2,names=maf_classes[2:11],xlab="Freq bin", ylab="R2",main="R2 from IMPUTE",varwidth=TRUE)
# dev.off()

#jpeg("scatter_r2.jpg", width=500, height=500)
jpeg("scatter_r2_mean_compare.jpg", width=800, height=800)
par(lab=c(length(maf_classes[2:length(maf_classes)]),5,length(maf_classes[2:length(maf_classes)])))
plot(bin_1_r2_average, type='c',lty=2,lwd=3, main="R2 in different frequence bins", ylab="R2",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
points(bin_1_r2_average,pch=23,cex=3, lwd=1,bg="red")
# plot(bin_r2_average, pch=19, main="R2 in different frequence bins", ylab="R2",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
# lines(bin_r2_imp_average,type='p', pch=19, col='blue')
lines(bin_2_r2_average,type='c',lty=2,lwd=3, col='purple')
points(bin_2_r2_average,pch=23,cex=3, lwd=1,bg='purple')
lines(bin_3_r2_average,type='c',lty=2,lwd=3, col='green')
points(bin_3_r2_average,pch=23,cex=3, lwd=1,bg='green')
lines(bin_4_r2_average,type='c',lty=2,lwd=3, col='blue')
points(bin_4_r2_average,pch=25,cex=3, lwd=1,bg='yellow')
# lines(bin_r2_imp_median,type='p', pch=19, col='purple')
axis(1, at=1:(length(maf_classes)-1), labels=maf_classes[2:length(maf_classes)])
# legend(6, 0.2, c("r2 average","r2 median"), cex=0.8, col=c("red","purple"), pch=19)
legend(6, 0.2, c("1000GP panel","1000GP + VBI panel","VBI+UK10K panel","1000GP+UK10K panel"), cex=1.5, col=c("red","purple","green","blue"), lty=2,pch=c(23,23,23,25),pt.bg=c("red","purple","green","yellow"))
dev.off()


bin_1_info_average <- NULL
bin_2_info_average <- NULL
bin_3_info_average <- NULL
bin_4_info_average <- NULL


for(i in 2:length(maf_classes)){
  bin_1_name <- paste("bin_1_i",maf_classes[i],sep="_")
  bin_2_name <- paste("bin_2_i",maf_classes[i],sep="_")
  bin_3_name <- paste("bin_3_i",maf_classes[i],sep="_")
  bin_4_name <- paste("bin_4_i",maf_classes[i],sep="_")

  current_1_bin <- read.table(paste("1KGP_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  current_1_bin$info <- as.numeric(as.character(current_1_bin$info))
  current_2_bin <- read.table(paste("VBI1kgp_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  current_2_bin$info <- as.numeric(as.character(current_2_bin$info))
  current_3_bin <- read.table(paste("VBIUK10K_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  current_3_bin$info <- as.numeric(as.character(current_3_bin$info))
  current_4_bin <- read.table(paste("UK10K1KP_panel_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  current_4_bin$info <- as.numeric(as.character(current_4_bin$info))

  assign(bin_1_name, current_1_bin)
  assign(bin_2_name, current_2_bin)
  assign(bin_3_name, current_3_bin)
  assign(bin_4_name, current_4_bin)
  bin_1_info_average <- c(bin_1_info_average,mean(current_1_bin$info ,na.rm=TRUE))
  bin_2_info_average <- c(bin_2_info_average,mean(current_2_bin$info,na.rm=TRUE))
  bin_3_info_average <- c(bin_3_info_average,mean(current_3_bin$info,na.rm=TRUE))
  bin_4_info_average <- c(bin_4_info_average,mean(current_4_bin$info,na.rm=TRUE))

}



for(i in 2:length(maf_classes)){
  bin_name <- paste("bin",maf_classes[i],sep="_")

  current_bin <- read.table(paste("VBI_1kg_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_bin <- current_bin[which(current_bin$concord_type0 != -1),]
  current_bin$info <- as.numeric(as.character(current_bin$info))
  
  assign(bin_name, current_bin)
  #bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
  #bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
  #bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
  #bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
  bin_r2_average <- c(bin_r2_average,mean(current_bin$info,na.rm=TRUE))
  bin_r2_median <- c(bin_r2_median,median(current_bin$info,na.rm=TRUE))
  # bin_r2_imp_average <- c(bin_r2_imp_average,mean(current_bin$IMP_R2,na.rm=TRUE))
  # bin_r2_imp_median <- c(bin_r2_imp_median,median(current_bin$IMP_R2,na.rm=TRUE))
  
}


jpeg("boxplot_info.jpg", width=800, height=800)
boxplot(bin_0.01$info,bin_0.02$info,bin_0.03$info,bin_0.04$info,bin_0.05$info,bin_0.1$info,bin_0.2$info,bin_0.3$info,bin_0.4$info,bin_0.5$info,names=maf_classes[2:11],xlab="Freq bin", ylab="INFO",main="INFO score",varwidth=TRUE)
dev.off()

# jpeg("boxplot_imp_r2.jpg", width=800, height=800)
# boxplot(bin_0.01$IMP_R2,bin_0.02$IMP_R2,bin_0.03$IMP_R2,bin_0.04$IMP_R2,bin_0.05$IMP_R2,bin_0.1$IMP_R2,bin_0.2$IMP_R2,bin_0.3$IMP_R2,bin_0.4$IMP_R2,bin_0.5$IMP_R2,names=maf_classes[2:11],xlab="Freq bin", ylab="R2",main="R2 from IMPUTE",varwidth=TRUE)
# dev.off()

jpeg("scatter_info.jpg", width=500, height=500)
par(lab=c(length(maf_classes[2:11]),5,length(maf_classes[2:11])))
plot(bin_r2_average, pch=19, main="INFO in different frequence bins", ylab="INFO",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
# lines(bin_r2_imp_average,type='p', pch=19, col='blue')
lines(bin_r2_median,type='p', pch=19, col='purple')
# lines(bin_r2_imp_median,type='p', pch=19, col='purple')
axis(1, at=1:10, labels=maf_classes[2:11])
legend(6, 0.2, c("INFO score average","INFO score median"), cex=0.8, col=c("red","purple"), pch=19)
dev.off()
