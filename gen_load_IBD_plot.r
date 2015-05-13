#r script to plot for enza
rm(list=ls())

######################################################################################################
# Plot 7: IBD
###############################################################################
#Same as roh, but we consider samples's pairs:
# input_format <- "PLINK"
chr <- "10"
# pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
pops <- c("CEU","TSI","VBI","CARL","Erto","Illegio","Resia","Sauris")
LOD <- 5
# setwd("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/ROH/20140811/BEAGLE/")
# base_folder <- getwd() #only if we're lazy
input_format <- "BEAGLE"
xmax <- NULL
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    # pop_ibd_file <- paste(pop,"roh.length.ibd",sep=".")
    # pop_ibd_file <- paste(pop,"roh.length.4.ibd",sep=".")
    pop_ibd_file <- paste("CHR",chr,"/",pop,".roh.length.",LOD,".ibd",sep="")
    pop_ibd_table <- read.table(pop_ibd_file,header=FALSE,sep=" ")
    colnames(pop_ibd_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","PAIR","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_ibd_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_ibd_table <- read.table(pop_ibd_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_ibd_table <- pop_ibd_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_ibd_table$IID2 <- as.character(pop_ibd_table$IID2)
  pop_ibd_table$IID1 <- as.character(pop_ibd_table$IID1)
  pop_ibd_table$PAIR <- as.character(pop_ibd_table$PAIR)
  current_pop_table_name <- paste(pop,"ibd",sep="_")
  assign(current_pop_table_name,pop_ibd_table)

  tot_ibd <- tapply(pop_ibd_table$LENGTH,pop_ibd_table$PAIR,sum)
  tot_ibd <- as.data.frame(tot_ibd, row.names=NULL)
  tot_ibd$ID <- rownames(tot_ibd)
  colnames(tot_ibd) <- c("IBD_tot","ID")
  tot_ibd$ID <- as.character(tot_ibd$ID)
  tot_ibd$IBD_tot <- as.numeric(as.character(tot_ibd$IBD_tot))
  tot_ibd$IBD_tot <- tot_ibd$IBD_tot/1000000

  assign(paste(pop,"tot_ibd",sep="_"),tot_ibd)
  
  assign(paste("M_ibd",pop,sep="_"),ecdf(tot_ibd$IBD_tot))
  xmax <- c(xmax,summary(ecdf(tot_ibd$IBD_tot))[6])

}

q_CEU_ibd <- quantile(CEU_tot_ibd$IBD_tot,c(.95))
q_TSI_ibd <- quantile(TSI_tot_ibd$IBD_tot,c(.95))
q_FVG_ibd <- quantile(FVG_tot_ibd$IBD_tot,c(.95))
q_VBI_ibd <- quantile(VBI_tot_ibd$IBD_tot,c(.95))
q_CARL_ibd <- quantile(CARL_tot_ibd$IBD_tot,c(.95))
q_Sauris_ibd <- quantile(Sauris_tot_ibd$IBD_tot,c(.95))
q_Erto_ibd <- quantile(Erto_tot_ibd$IBD_tot,c(.95))
q_Illegio_ibd <- quantile(Illegio_tot_ibd$IBD_tot,c(.95))
q_Resia_ibd <- quantile(Resia_tot_ibd$IBD_tot,c(.95))



#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.49898 
# TSI: 95% -> 1.242603 
# FVG: 95% -> 5.354979 
# VBI: 95% -> 2.017785 

#W50 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.473412
# TSI: 95% -> 1.32598
# FVG: 95% -> 5.255203
# VBI: 95% -> 1.994098

#W24 (Beagle RefinedIBD, filtered at LOD = 4)
# q_CEU_ibd : 95% -> 0.2779214 
# q_TSI_ibd : 95% -> 0.3357718 
# q_FVG_ibd : 95% -> 1.440876 
# q_VBI_ibd : 95% -> 0.6909859 
# q_CARL_ibd : 95% -> 0.5372137 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# q_CEU_ibd : 95% ->5.505728 
# q_TSI_ibd : 95% ->4.427433 
# q_FVG_ibd : 95% ->16.18415 
# q_VBI_ibd : 95% ->5.645495 
# q_CARL_ibd : 95% ->7.41046 
# q_Sauris_ibd : 95% ->35.33875 
# q_Erto_ibd : 95% ->25.70715 
# q_Illegio_ibd : 95% ->32.76607 
# q_Resia_ibd : 95% ->35.3291 


# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"/8_ibd.jpg",sep=""),width=800, height=800)
#   par(lwd=3,cex=1.8)
#   plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
#   lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col="red", verticals=TRUE, pch=46)
#   lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col="green", verticals=TRUE, pch=46)
#   lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col="blue", verticals=TRUE, pch=46)
#   lines(M_ibd_CARL,CARL_tot_ibd$IBD_tot,col="yellow", verticals=TRUE, pch=46)
#   legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI","CARL"),col=c("black","red","green","blue","yellow"),ncol=5)
# dev.off()


jpeg(paste(base_folder,"/8_ibd_all_5POP_lod5_",chr,".jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=1.5)
  plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col="green", verticals=TRUE, pch=46)
  lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col="blue", verticals=TRUE, pch=46)
  lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col="red", verticals=TRUE, pch=46)
  lines(M_ibd_CARL,CARL_tot_ibd$IBD_tot,col="yellow", verticals=TRUE, pch=46)
  lines(M_ibd_Erto,Erto_tot_ibd$IBD_tot,col="darkred", verticals=TRUE, pch=46)
  lines(M_ibd_Illegio,Illegio_tot_ibd$IBD_tot,col="indianred", verticals=TRUE, pch=46)
  lines(M_ibd_Resia,Resia_tot_ibd$IBD_tot,col="mediumvioletred", verticals=TRUE, pch=46)
  lines(M_ibd_Sauris,Sauris_tot_ibd$IBD_tot,col="orangered", verticals=TRUE, pch=46)
  legend("bottomright",pch =c(rep(19,length(pops))),legend=pops,col=c("black","green","blue","red","yellow","darkred","indianred","mediumvioletred","orangered"),ncol=length(pops))
dev.off()

#plot the ROH cumulative dstribution (for length ?)

######################################################################################################
# Plot 7_A_IBD GENOMEWIDE for each POPULATION
###############################################################################
#Same as roh, but we consider samples's pairs:
# input_format <- "PLINK"
# pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
rm(list=ls())
require(ggplot2)
require(reshape2)

pops <- c("CEU","TSI","CARL","VBI","Erto","Illegio","Resia","Sauris")
pops_c <- c("CEU","TSI","CAR","VBI","FVG-E","FVG-I","FVG-R","FVG-S")
LOD <- 5
# setwd("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/ROH/20140811/BEAGLE/")
base_folder <- getwd() #only if we're lazy
input_format <- "BEAGLE"
xmax <- NULL
for (pop in pops) {

  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    pop_ibd_file <- paste(pop,"/",pop,".WG.roh.length.",LOD,".ibd",sep="")
    pop_ibd_table <- read.table(pop_ibd_file,header=FALSE,sep=" ")
    colnames(pop_ibd_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","PAIR","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_ibd_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_ibd_table <- read.table(pop_ibd_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_ibd_table <- pop_ibd_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_ibd_table$IID2 <- as.character(pop_ibd_table$IID2)
  pop_ibd_table$IID1 <- as.character(pop_ibd_table$IID1)
  pop_ibd_table$PAIR <- as.character(pop_ibd_table$PAIR)
  current_pop_table_name <- paste(pop,"ibd",sep="_")
  assign(current_pop_table_name,pop_ibd_table)

  tot_ibd <- tapply(as.numeric(pop_ibd_table$LENGTH),pop_ibd_table$PAIR,sum)
  tot_ibd <- as.data.frame(tot_ibd, row.names=NULL)
  tot_ibd$ID <- rownames(tot_ibd)
  colnames(tot_ibd) <- c("IBD_tot","ID")
  tot_ibd$ID <- as.character(tot_ibd$ID)
  tot_ibd$IBD_tot <- as.numeric(as.character(tot_ibd$IBD_tot))
  tot_ibd$IBD_tot <- tot_ibd$IBD_tot/1000000

  assign(paste(pop,"tot_ibd",sep="_"),tot_ibd)
  
  assign(paste("M_ibd",pop,sep="_"),ecdf(tot_ibd$IBD_tot))
  xmax <- c(xmax,summary(ecdf(tot_ibd$IBD_tot))[6])

}

# q_FVG_ibd <- quantile(FVG_tot_ibd$IBD_tot,c(.95))
q_CEU_ibd <- quantile(CEU_tot_ibd$IBD_tot,c(.95))
q_TSI_ibd <- quantile(TSI_tot_ibd$IBD_tot,c(.95))
q_VBI_ibd <- quantile(VBI_tot_ibd$IBD_tot,c(.95))
q_CARL_ibd <- quantile(CARL_tot_ibd$IBD_tot,c(.95))
q_Sauris_ibd <- quantile(Sauris_tot_ibd$IBD_tot,c(.95))
q_Erto_ibd <- quantile(Erto_tot_ibd$IBD_tot,c(.95))
q_Illegio_ibd <- quantile(Illegio_tot_ibd$IBD_tot,c(.95))
q_Resia_ibd <- quantile(Resia_tot_ibd$IBD_tot,c(.95))


# summary_FVG_ibd <- summary(FVG_tot_ibd$IBD_tot)
summary_CEU_ibd <- summary(CEU_tot_ibd$IBD_tot)
summary_TSI_ibd <- summary(TSI_tot_ibd$IBD_tot)
summary_VBI_ibd <- summary(VBI_tot_ibd$IBD_tot)
summary_CARL_ibd <- summary(CARL_tot_ibd$IBD_tot)
summary_Sauris_ibd <- summary(Sauris_tot_ibd$IBD_tot)
summary_Erto_ibd <- summary(Erto_tot_ibd$IBD_tot)
summary_Illegio_ibd <- summary(Illegio_tot_ibd$IBD_tot)
summary_Resia_ibd <- summary(Resia_tot_ibd$IBD_tot)

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.49898 
# TSI: 95% -> 1.242603 
# FVG: 95% -> 5.354979 
# VBI: 95% -> 2.017785 

#W50 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.473412
# TSI: 95% -> 1.32598
# FVG: 95% -> 5.255203
# VBI: 95% -> 1.994098

#W24 (Beagle RefinedIBD, filtered at LOD = 4)
# q_CEU_ibd : 95% -> 0.2779214 
# q_TSI_ibd : 95% -> 0.3357718 
# q_FVG_ibd : 95% -> 1.440876 
# q_VBI_ibd : 95% -> 0.6909859 
# q_CARL_ibd : 95% -> 0.5372137 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# q_CEU_ibd : 95% ->5.505728 
# q_TSI_ibd : 95% ->4.427433 
# q_FVG_ibd : 95% ->16.18415 
# q_VBI_ibd : 95% ->5.645495 
# q_CARL_ibd : 95% ->7.41046 
# q_Sauris_ibd : 95% ->35.33875 
# q_Erto_ibd : 95% ->25.70715 
# q_Illegio_ibd : 95% ->32.76607 
# q_Resia_ibd : 95% ->35.3291 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)GENOME WIDE 1/8/2014
# q_CEU_ibd: 95% -> 80.43721 
# q_TSI_ibd: 95% -> 61.59844 
# q_FVG_ibd: 95% -> 246.3632 
# q_VBI_ibd: 95% -> 61.37638 
# q_CARL_ibd: 95% -> 91.72759 
# q_Sauris_ibd: 95% -> 443.422 
# q_Erto_ibd: 95% -> 363.6143 
# q_Illegio_ibd: 95% -> 449.9667 
# q_Resia_ibd: 95% -> 412.6811 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)GENOME WIDE 5/8/2014 ONLY 46 SAMPLES
#/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/ROH/20140804/BEAGLE
# q_CEU_ibd: 95% -> 136.0912 
# q_TSI_ibd: 95% -> 122.4737 
# q_FVG_ibd: 95% -> 328.0547 
# q_VBI_ibd: 95% -> 162.723 
# q_CARL_ibd: 95% -> 157.1352 
# q_Sauris_ibd: 95% -> 461.4193 
# q_Erto_ibd: 95% -> 337.5515 
# q_Illegio_ibd: 95% -> 444.8172 
# q_Resia_ibd: 95% -> 429.6575 

# IBD           
# POP Min.  1st_Qu. Median  Mean  3rd_Qu. Max.
# summary_CEU_ibd 87.09 111.7 118.8 119.2 125.4 312.7
# summary_TSI_ibd 85.98 102.1 107.4 108.2 113.3 144.4
# summary_FVG_ibd 36.17 117.7 127.6 158.5 147.2 933.8
# summary_VBI_ibd 61.24 125.2 133.8 134.8 143.9 519.1
# summary_CARL_ibd  50.79 83.86 94.51 104 109.5 612.5
# summary_Sauris_ibd  84.33 212.1 258 275.3 315 1242
# summary_Erto_ibd  95.73 171.1 226.7 229.3 269 1151
# summary_Illegio_ibd 49.36 176.4 256.8 261 320 1101
# summary_Resia_ibd 89.86 247.4 292.5 301.8 336.6 2521

source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
pop_colors <- col_pop(pops_c)
base_folder <- getwd()



# jpeg(paste(base_folder,"/7_A_ibd_all_5POP_lod5_WG_1500.jpg",sep=""),width=1000, height=1000)
jpeg(paste(base_folder,"/7_A_ibd_all_5POP_lod5_WG_1500_NO_FVG.jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=2)
  # plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  # lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "FVG"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "CEU"),1],main="",xlab="IBD  genome per pair of individuals (Mb)",ylab="Cumulative frequency", xlim=c(0,1000), verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "TSI"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "VBI"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_CARL,CARL_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "CARL"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Erto,Erto_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Erto"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Illegio,Illegio_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Illegio"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Resia,Resia_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Resia"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Sauris,Sauris_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Sauris"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  abline(h=0.95,col='grey',lty='dashed')
  # leg_txt <- c(pop_colors[which(pop_colors$pop == "CEU"),2],pop_colors[which(pop_colors$pop == "TSI"),2],pop_colors[which(pop_colors$pop == "VBI"),2],pop_colors[which(pop_colors$pop == "FVG"),2],pop_colors[which(pop_colors$pop == "CARL"),2],pop_colors[which(pop_colors$pop == "Erto"),2],pop_colors[which(pop_colors$pop == "Illegio"),2],pop_colors[which(pop_colors$pop == "Resia"),2],pop_colors[which(pop_colors$pop == "Sauris"),2])
  leg_txt <- c(pop_colors[which(pop_colors$pop == "CEU"),2],pop_colors[which(pop_colors$pop == "TSI"),2],pop_colors[which(pop_colors$pop == "VBI"),2],pop_colors[which(pop_colors$pop == "CARL"),2],pop_colors[which(pop_colors$pop == "Erto"),2],pop_colors[which(pop_colors$pop == "Illegio"),2],pop_colors[which(pop_colors$pop == "Resia"),2],pop_colors[which(pop_colors$pop == "Sauris"),2])
  # bkg <- c(pop_colors[which(pop_colors$pop == "CEU"),1],pop_colors[which(pop_colors$pop == "TSI"),1],pop_colors[which(pop_colors$pop == "VBI"),1],pop_colors[which(pop_colors$pop == "FVG"),1],pop_colors[which(pop_colors$pop == "CARL"),1],pop_colors[which(pop_colors$pop == "Erto"),1],pop_colors[which(pop_colors$pop == "Illegio"),1],pop_colors[which(pop_colors$pop == "Resia"),1],pop_colors[which(pop_colors$pop == "Sauris"),1])
  bkg <- c(pop_colors[which(pop_colors$pop == "CEU"),1],pop_colors[which(pop_colors$pop == "TSI"),1],pop_colors[which(pop_colors$pop == "VBI"),1],pop_colors[which(pop_colors$pop == "CARL"),1],pop_colors[which(pop_colors$pop == "Erto"),1],pop_colors[which(pop_colors$pop == "Illegio"),1],pop_colors[which(pop_colors$pop == "Resia"),1],pop_colors[which(pop_colors$pop == "Sauris"),1])
  legend("bottomright",pch =c(rep(22,length(pops))),legend=leg_txt, pt.lwd=2,pt.cex=2,pt.bg=bkg,col=c(rep('black',length(pops))),ncol=4,bty="n")
dev.off()

#############################################################################################
#Replot everything with a ggplot function
#Create a data frame with Ibd values
ibd_all_cum <- data.frame(IBD_tot=c(CEU_tot_ibd$IBD_tot,
  TSI_tot_ibd$IBD_tot,
  CARL_tot_ibd$IBD_tot,
  VBI_tot_ibd$IBD_tot,
  Erto_tot_ibd$IBD_tot,
  Illegio_tot_ibd$IBD_tot,
  Resia_tot_ibd$IBD_tot,
  Sauris_tot_ibd$IBD_tot),
 pop=rep(pops_c,rep(dim(CEU_tot_ibd)[1],8)),
 ecdf=c(M_ibd_CEU(CEU_tot_ibd$IBD_tot),
  M_ibd_TSI(TSI_tot_ibd$IBD_tot),
  M_ibd_CARL(CARL_tot_ibd$IBD_tot),
  M_ibd_VBI(VBI_tot_ibd$IBD_tot),
  M_ibd_Erto(Erto_tot_ibd$IBD_tot),
  M_ibd_Illegio(Illegio_tot_ibd$IBD_tot),
  M_ibd_Resia(Resia_tot_ibd$IBD_tot),
  M_ibd_Sauris(Sauris_tot_ibd$IBD_tot)))

#manually calculate ecdf (no ecdf function used)
# ibd_all_cum <- ibd_all_cum[order(ibd_all_cum$pop),]
ibd_all_cum$pop <- factor(ibd_all_cum$pop,level=names(pop_colors))
# ibd_all_cum$ecdf <- ave(ibd_all_cum$IBD_tot, ibd_all_cum$pop, FUN=function(IBD_tot) seq_along(IBD_tot)/length(IBD_tot))

#plot 

pl <- ggplot(ibd_all_cum)
pl <- pl + aes(x = IBD_tot, y = ecdf, colour=pop)
pl <- pl + scale_color_manual("Cohorts", values=pop_colors)
pl <- pl + geom_line(size=1.5)
pl <- pl + geom_hline(aes(yintercept=0.95), linetype=2,colour="Lightgrey",size=1.2)
pl <- pl + xlab("IBD genome per pair of individuals (Mb)") + ylab("Cumulative frequency")
pl <- pl + scale_x_continuous(limits=c(0,1250))
pl <- pl + scale_y_continuous(limits=c(0,1))
pl <- pl + theme_bw()
pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))

ggsave(filename=paste(base_folder,"/7_A_ibd_all_5POP_lod5_WG_1200_NO_FVG_ggplot.jpeg",sep=""),width=8, height=8,dpi=400,plot=pl)
# ggsave(filename=paste(base_folder,"/test_IBD_2.jpeg",sep=""),width=8, height=8,dpi=400,plot=pl)

# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/7_roh.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"PLOTS/7_roh_density.jpg",sep=""),width=800, height=800)
  plot(density(VBI_tot_roh$ROH_tot),main="",xlab="Total ROH homozigosity", col="blue")
  lines(density(FVG_tot_roh$ROH_tot),col="red")
  lines(density(TSI_tot_roh$ROH_tot),col="green")
  lines(density(CEU_tot_roh$ROH_tot))
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()

