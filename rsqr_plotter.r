#Script and function to plot rsquared correlation splitted in frequency bins

# basepath <- "/lustre/scratch110/sanger/mc14/GENOTIPI/TEST_IMPUTATION"
basepath <- "/lustre/scratch113/teams/soranzo/users/mc14/GENOTIPI/SEQ_IMP_RSQR"

# output3_info <- read.table(paste(basepath,"OUPUT3/MERGED/SNPS/MACH_FORMAT/0/16.machinfo.mod.new",sep="/"),sep="\t",header=T)

# summary(output3_info)

# output3_rsqr <- read.table(paste(basepath,"OUPUT3/MERGED/SNPS/MACH_FORMAT/0/RSQR/chr16_vbi.doseR2",sep="/"),header=T,sep=" ")

# summary(output3_rsqr)

# output3_merged <- merge(output3_rsqr,output3_info,by.x="SNP",by.y="name",all=FALSE,sort=FALSE)

# summary(output3_merged)

# output3_info <- read.table(paste(basepath,"OUPUT3/MERGED/RSQUARED_NEW/chr16.complete_R2_uniq",sep="/"),sep=" ",header=T)
# output4_info <- read.table(paste(basepath,"OUPUT4/MERGED/RSQUARED_NEW/chr16.complete_R2_uniq",sep="/"),sep=" ",header=T)
# output5_info <- read.table(paste(basepath,"OUPUT5/MERGED/RSQUARED_NEW/chr16.complete_R2_uniq",sep="/"),sep=" ",header=T)

panel1 <- read.table(paste(basepath,"VBI_1000GP/chr7.complete_R2_uniq",sep="/"),sep=" ",header=T)
panel2 <- read.table(paste(basepath,"VBI_UK10K/chr7.complete_R2_uniq",sep="/"),sep=" ",header=T)
panel3 <- read.table(paste(basepath,"UK10K_1000GP/chr7.complete_R2_uniq",sep="/"),sep=" ",header=T)
panel4 <- read.table(paste(basepath,"IMP/chr7.complete_R2_uniq",sep="/"),sep=" ",header=T)

#define maf classes
# maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
# maf_classes <- c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
maf_classes <- c(0,0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
source('~/Work/r_scripts/maf_bins_splitter.r')
# split_bins(maf_classes,output3_info,"OUT3")
# split_bins(maf_classes,output4_info,"OUT4")
# split_bins(maf_classes,output5_info,"OUT5")

# jpeg("scatter_r2.jpg", width=1000, height=1000)
# plot(output3_info$R2,lwd=2,pch=20,col="red")
# points(output4_info$R2,lwd=2,pch=21,col="green")
# points(output5_info$R2,lwd=2,pch=22,col="deepskyblue3")
# dev.off()


# #some addition to plot other stuff
# bin_average <- NULL
# bin_median <- NULL
# for(i in 2:length(maf_classes)){
#   bin_name <- paste("bin",maf_classes[i],sep="_")

#   current_bin <- read.table(paste(basepath,"/SUMMARIES/IMP_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  
#   assign(bin_name, current_bin)
#   bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
#   bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
  
# }

##########################

# basepath <- "/lustre/scratch110/sanger/mc14/GENOTIPI/SEQ_IMP_RSQR/RSQUARED"
# chr7_info2 <- read.table(paste(basepath,"chr7.complete_R2_2_uniq",sep="/"),sep=" ",header=T)
# maf_classes <- c(0,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)
# split_bins(maf_classes,chr7_info2,"IMP2")
split_bins(maf_classes,panel1,"VBI_1000GP")
split_bins(maf_classes,panel2,"VBI_UK10K")
split_bins(maf_classes,panel3,"UK10K_1000GP")
split_bins(maf_classes,panel4,"1000GP")

# bin_r2_average <- NULL
# bin_r2_median <- NULL
# bin_r2_imp_average <- NULL
# bin_r2_imp_median <- NULL
# for(i in 2:length(maf_classes)){
#   bin_name <- paste("bin",maf_classes[i],sep="_")

#   current_bin <- read.table(paste(basepath,"/SUMMARIES2/IMP2_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
#   current_bin$R2 <- as.numeric(as.character(current_bin$R2))
  
#   assign(bin_name, current_bin)
#   #bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
#   #bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
#   #bin_average <- c(bin_average,mean(current_bin$R2,na.rm=TRUE))
#   #bin_median <- c(bin_median,median(current_bin$R2,na.rm=TRUE))
#   bin_r2_average <- c(bin_r2_average,mean(current_bin$R2,na.rm=TRUE))
#   bin_r2_median <- c(bin_r2_median,median(current_bin$R2,na.rm=TRUE))
#   bin_r2_imp_average <- c(bin_r2_imp_average,mean(current_bin$IMP_R2,na.rm=TRUE))
#   bin_r2_imp_median <- c(bin_r2_imp_median,median(current_bin$IMP_R2,na.rm=TRUE))
  
# }


# jpeg("boxplot_r2.jpg", width=800, height=800)
# boxplot(bin_0.01$R2,bin_0.02$R2,bin_0.03$R2,bin_0.04$R2,bin_0.05$R2,bin_0.1$R2,bin_0.2$R2,bin_0.3$R2,bin_0.4$R2,bin_0.5$R2,names=maf_classes[2:11],xlab="Freq bin", ylab="R2",main="R2 between sequenced and imputed genotypes on 110 samples",varwidth=TRUE)
# dev.off()

# jpeg("boxplot_imp_r2.jpg", width=800, height=800)
# boxplot(bin_0.01$IMP_R2,bin_0.02$IMP_R2,bin_0.03$IMP_R2,bin_0.04$IMP_R2,bin_0.05$IMP_R2,bin_0.1$IMP_R2,bin_0.2$IMP_R2,bin_0.3$IMP_R2,bin_0.4$IMP_R2,bin_0.5$IMP_R2,names=maf_classes[2:11],xlab="Freq bin", ylab="R2",main="R2 from IMPUTE",varwidth=TRUE)
# dev.off()

# jpeg("scatter_r2.jpg", width=500, height=500)
# par(lab=c(length(maf_classes[2:11]),5,length(maf_classes[2:11])))
# plot(bin_r2_average, pch=19, main="R2 in different frequence bins", ylab="R2",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
# lines(bin_r2_imp_average,type='p', pch=19, col='blue')
# lines(bin_r2_median,type='p', pch=19, col='green')
# lines(bin_r2_imp_median,type='p', pch=19, col='purple')
# axis(1, at=1:10, labels=maf_classes[2:11])
# legend(6, 0.2, c("IMP vs WGS r2 average","IMP vs GWAS r2 average","IMP vs WGS r2 median","IMP vs GWAS r2 median"), cex=0.8, col=c("red","blue","green","purple"), pch=19)
# dev.off()

bin_1_r2_average <- NULL
bin_2_r2_average <- NULL
bin_3_r2_average <- NULL
bin_4_r2_average <- NULL

for(i in 2:length(maf_classes)){
  bin_1_name <- paste("bin_1",maf_classes[i],sep="_")
  bin_2_name <- paste("bin_2",maf_classes[i],sep="_")
  bin_3_name <- paste("bin_3",maf_classes[i],sep="_")
  bin_4_name <- paste("bin_4",maf_classes[i],sep="_")

  current_1_bin <- read.table(paste("VBI_1000GP_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_1_bin <- current_1_bin[which(current_1_bin$concord_type0 != -1),]
  current_1_bin$R2 <- as.numeric(as.character(current_1_bin$R2))
  current_2_bin <- read.table(paste("VBI_UK10K_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_2_bin <- current_2_bin[which(current_2_bin$concord_type0 != -1),]
  current_2_bin$R2 <- as.numeric(as.character(current_2_bin$R2))
  current_3_bin <- read.table(paste("UK10K_1000GP_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_3_bin <- current_3_bin[which(current_3_bin$concord_type0 != -1),]
  current_3_bin$R2 <- as.numeric(as.character(current_3_bin$R2))
  current_4_bin <- read.table(paste("1000GP_maf_lte_",maf_classes[i],"_table.txt",sep=""),sep="\t",header=T)
  # current_4_bin <- current_4_bin[which(current_4_bin$concord_type0 != -1),]
  current_4_bin$R2 <- as.numeric(as.character(current_4_bin$R2))
  
  assign(bin_1_name, current_1_bin)
  assign(bin_2_name, current_2_bin)
  assign(bin_3_name, current_3_bin)
  assign(bin_4_name, current_4_bin)
  bin_1_r2_average <- c(bin_1_r2_average,mean(current_1_bin$R2,na.rm=TRUE))
  bin_2_r2_average <- c(bin_2_r2_average,mean(current_2_bin$R2,na.rm=TRUE))
  bin_3_r2_average <- c(bin_3_r2_average,mean(current_3_bin$R2,na.rm=TRUE))
  bin_4_r2_average <- c(bin_4_r2_average,mean(current_4_bin$R2,na.rm=TRUE))
  
}


#OBSERVED R2 calculation plot
jpeg("scatter_OBSERVED_r2_mean_compare.jpg", width=800, height=800)
par(lab=c(length(maf_classes[2:length(maf_classes)]),5,length(maf_classes[2:length(maf_classes)])))
plot(bin_4_r2_average, type='c',lty=2,lwd=3, main="R2 in different frequence bins", ylab="R2",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
points(bin_4_r2_average,pch=23,cex=3, lwd=1,bg="red")
# plot(bin_r2_average, pch=19, main="R2 in different frequence bins", ylab="R2",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
# lines(bin_r2_imp_average,type='p', pch=19, col='blue')
lines(bin_1_r2_average,type='c',lty=2,lwd=3, col='purple')
points(bin_1_r2_average,pch=23,cex=3, lwd=1,bg='purple')
lines(bin_2_r2_average,type='c',lty=2,lwd=3, col='green')
points(bin_2_r2_average,pch=23,cex=3, lwd=1,bg='green')
lines(bin_3_r2_average,type='c',lty=2,lwd=3, col='blue')
points(bin_3_r2_average,pch=25,cex=3, lwd=1,bg='yellow')
# lines(bin_r2_imp_median,type='p', pch=19, col='purple')
axis(1, at=1:(length(maf_classes)-1), labels=maf_classes[2:length(maf_classes)])
# legend(6, 0.2, c("r2 average","r2 median"), cex=0.8, col=c("red","purple"), pch=19)
legend(6, 0.2, c("1000GP panel","1000GP + VBI panel","VBI+UK10K panel","1000GP+UK10K panel"), cex=1.5, col=c("red","purple","green","blue"), lty=2,pch=c(23,23,23,25),pt.bg=c("red","purple","green","yellow"))
dev.off()


jpeg("scatter_R2ONE_mean_compare.jpg", width=800, height=800)
par(lab=c(length(maf_classes[2:length(maf_classes)]),5,length(maf_classes[2:length(maf_classes)])))
plot(bin_3_r2_average, type='c',lty=2,lwd=3, main="R2 in different frequence bins", ylab="R2",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
points(bin_3_r2_average,pch=23,cex=3, lwd=1,bg="red")
# plot(bin_r2_average, pch=19, main="R2 in different frequence bins", ylab="R2",xlab="Frequence bins",xaxt="n",ylim=c(0,1),col='red')
# lines(bin_r2_imp_average,type='p', pch=19, col='blue')
# lines(bin_r2_imp_median,type='p', pch=19, col='purple')
axis(1, at=1:(length(maf_classes)-1), labels=maf_classes[2:length(maf_classes)])
# legend(6, 0.2, c("r2 average","r2 median"), cex=0.8, col=c("red","purple"), pch=19)
legend(6, 0.2, c("1000GP panel","1000GP + VBI panel","VBI+UK10K panel","1000GP+UK10K panel"), cex=1.5, col=c("red","purple","green","blue"), lty=2,pch=c(23,23,23,25),pt.bg=c("red","purple","green","yellow"))
dev.off()
