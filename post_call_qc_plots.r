#Script to generate plots for QC

#1)Plots of HET vs NRD
#1a) By sample

rm(list=ls())


pop <- "VBI"
pop <- "FVG"

if (pop == "VBI"){
  pop_col <- "blue"
  #read data by sample
  # pop_NRD_file <- "../../DISCORDANCE/QC_OUT/all_sample_concordance_discordance.txt"
  pop_NRD_file <- "all_sample_concordance_discordance.txt"
  pop_HET_file <-"wgs_all_chr_updated_no_triallelic_het_rate.het"

}else{
  pop_col <- "red"
  #read data by sample
  pop_HET_file <-"wgs_all_chr_cleaned_no_mismatch_het_rate.het"
  # pop_NRD_file <- "../../20140522_WGS_GWAS_DISC/QC_OUT/all_sample_concordance_discordance.txt"
  pop_NRD_file <- "all_sample_concordance_discordance.txt"
}

pop_NRD <- read.table(pop_NRD_file,header=T)
pop_HET <- read.table(pop_HET_file,header=T)

#merge the data in a single file
pop_HET_NRD <- merge(pop_HET,pop_NRD, by.x="IID",by.y="SAMPLE_ID",all=T)
pop_HET_NRD$HET_RATE <- (pop_HET_NRD$HET_RATE)*100
pop_HET_NRD$NRD <- (pop_HET_NRD$NRD)*100


#plot HET vs NRD in a fancy way with sd lines for both measures
jpeg(paste(pop,"HET_NRD_by_sample.jpg",sep="_"),width=800, height=800)
par(cex=2,lwd=3)
plot(pop_HET_NRD$HET_RATE,pop_HET_NRD$NRD,main=paste("HET vs NRD by sample in ",pop,sep=""),ylab="discordance (%)",xlab="het rate (%)",col=pop_col)
#stats for NRD
abline(h=mean(pop_HET_NRD$NRD)+2*sd(pop_HET_NRD$NRD),col="red",lty="dotdash")
abline(h=mean(pop_HET_NRD$NRD)+3*sd(pop_HET_NRD$NRD),col="darkred",lty="dotdash")
# abline(h=mean(pop_HET_NRD$NRD),col="red",lty="dotdash")
# abline(h=mean(pop_HET_NRD$NRD)-sd(pop_HET_NRD$NRD),col="red",lty="dotdash")

#stats for HET
abline(v=mean(pop_HET_NRD$HET)+2*sd(pop_HET_NRD$HET),col="green",lty="dotdash")
abline(v=mean(pop_HET_NRD$HET)+3*sd(pop_HET_NRD$HET),col="darkgreen",lty="dotdash")
# abline(v=mean(pop_HET_NRD$HET),col="green",lty="dotdash")
# abline(v=mean(pop_HET_NRD$HET)-sd(pop_HET_NRD$HET),col="green",lty="dotdash")
dev.off()

#1b) by site: the HET file from plink is based on rsIDs, while the NRD is based on position. We need a map file to have complete information
pop_map_file <- "/home/max/Work/valborbera/NEW_SEQ/POST_REF_QC/WGS/wgs_all_chr.bim"
pop_site_HET_file <- "/home/max/Work/valborbera/NEW_SEQ/POST_REF_QC/WGS/wgs_all_chr_hardy.hwe"
pop_site_NRD_file <- "/home/max/Work/valborbera/NEW_SEQ/POST_REF_QC/WGS/DISCORDANCE/QC_OUT/all_sites_all_chr_concordance_discordance_table.txt"

pop_site_HET <- read.table(pop_site_HET_file,header=T,strip.white=T) 
pop_site_NRD <- read.table(pop_site_NRD_file,header=T,strip.white=T)
pop_map <- read.table(pop_map_file,header=F,strip.white=T)
pop_map$V3 <- NULL
pop_map$V5 <- NULL
pop_map$V6 <- NULL

colnames(pop_map) <- c("CHR","SNP","POS")

#merge HET datset with map info and order back with chr and pos
pop_site_HET_pos <- merge(pop_site_HET,pop_map,by.x="SNP",by.y="SNP",all=T)
pop_site_HET_pos <- (pop_site_HET_pos[order(pop_site_HET_pos$CHR.x,pop_site_HET_pos$POS),])

#we need to work by chromosome in this, and plot by chr too
for(chr in 1:22 ){
  current_chr_HET <- pop_site_HET_pos[which(pop_site_HET_pos$CHR.x == chr),]
  current_chr_NRD <- pop_site_NRD[which(pop_site_NRD$CHR == chr),]

  current_chr_merge <- merge(current_chr_HET,current_chr_NRD,by.x="POS",by.y="POS",all=T)

  summary(current_chr_merge)

  jpeg(paste(pop,"HET_NRD_by_site_chr",chr,".jpg",sep="_"),width=800, height=800,pointsize = 10)
  plot(current_chr_merge$NRD,current_chr_merge$O.HET.,main=paste("HET vs NRD by site in",pop,"chr",chr,sep=" "),xlab="NRD",ylab="HET RATE")
  #stats for NRD
  abline(v=mean(current_chr_merge$NRD)+sd(current_chr_merge$NRD),col="red",lty="dotdash")
  abline(v=mean(current_chr_merge$NRD),col="red",lty="dotdash")
  abline(v=mean(current_chr_merge$NRD)-sd(current_chr_merge$NRD),col="red",lty="dotdash")

  #stats for HET
  abline(h=mean(current_chr_merge$O.HET.)+sd(current_chr_merge$O.HET.),col="green",lty="dotdash")
  abline(h=mean(current_chr_merge$O.HET.),col="green",lty="dotdash")
  abline(h=mean(current_chr_merge$O.HET.)-sd(current_chr_merge$O.HET.),col="green",lty="dotdash")
  dev.off()  

  #also plot HET by POS
  jpeg(paste(pop,"HET_POS_chr",chr,".jpg",sep="_"),width=800, height=800,pointsize = 10)
  plot(current_chr_merge$POS,current_chr_merge$O.HET.,main=paste("HET vs POS by site in",pop,"chr",chr,sep=" "),xlab="POS",ylab="HET RATE")
  #stats for HET
  abline(h=mean(current_chr_merge$O.HET.)+sd(current_chr_merge$O.HET.),col="green",lty="dotdash")
  abline(h=mean(current_chr_merge$O.HET.),col="green",lty="dotdash")
  abline(h=mean(current_chr_merge$O.HET.)-sd(current_chr_merge$O.HET.),col="green",lty="dotdash")
  dev.off()  
}

#Pie chart plot for different sites category:
library(plotrix)
# library(ggplot2)

pop <- "FVG"
pop <- "VBI"

FVG_pie <- as.data.frame(cbind(11568618,1516565,1214932,73880))
VBI_pie <- as.data.frame(cbind(14141052,3266988,1202096,99424))

colnames(FVG_pie) <- colnames(VBI_pie) <- c("SNPs","SiS_SNPs","INDELs","SiS_INDELs")

# assign(paste(pop,"pie",sep="_"), pop_maf_resume)

