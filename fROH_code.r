###########################################################################################
#Plot 6_A_ROH GENOMEWIDE BY POPULATION
###########################################################################################

rm(list=ls())
require(ggplot2)
require(reshape2)

source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
pops <- c("CEU","TSI","CARL","VBI","Erto","Illegio","Resia","Sauris")
pops_c <- c("CEU","TSI","CAR","VBI","FVG-E","FVG-I","FVG-R","FVG-S")
input_format <- "BEAGLE"
xmax <- NULL
LOD <- 5
setwd("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/ROH/20140811/BEAGLE/")
base_folder <- getwd() #only if we're lazy
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    pop_roh_file <- paste(pop,"/",pop,".WG.roh.length.",LOD,".hbd",sep="")
    pop_roh_table <- read.table(pop_roh_file,header=FALSE,sep=" ")
    colnames(pop_roh_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_roh_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_roh_table <- read.table(pop_roh_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_roh_table <- pop_roh_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_roh_table$IID2 <- as.character(pop_roh_table$IID2)
  pop_roh_table$IID1 <- as.character(pop_roh_table$IID1)
  current_pop_table_name <- paste(pop,"roh",sep="_")
  assign(current_pop_table_name,pop_roh_table)

  tot_roh <- tapply(pop_roh_table$LENGTH,pop_roh_table$IID1,sum)
  tot_roh <- as.data.frame(tot_roh, row.names=NULL)
  tot_roh$ID <- rownames(tot_roh)
  colnames(tot_roh) <- c("ROH_tot","ID")
  tot_roh$ID <- as.character(tot_roh$ID)
  tot_roh$ROH_tot <- as.numeric(as.character(tot_roh$ROH_tot))
  tot_roh$ROH_tot <- tot_roh$ROH_tot/1000000

  assign(paste(pop,"tot_roh",sep="_"),tot_roh)
  
  assign(paste("M",pop,sep="_"),ecdf(tot_roh$ROH_tot))
  xmax <- c(xmax,summary(ecdf(tot_roh$ROH_tot))[6])

}

# q_FVG <- quantile(FVG_tot_roh$ROH_tot,c(.95))
q_CEU <- quantile(CEU_tot_roh$ROH_tot,c(.95))
q_TSI <- quantile(TSI_tot_roh$ROH_tot,c(.95))
q_VBI <- quantile(VBI_tot_roh$ROH_tot,c(.95))
q_CARL <- quantile(CARL_tot_roh$ROH_tot,c(.95))
q_Sauris <- quantile(Sauris_tot_roh$ROH_tot,c(.95))
q_Erto <- quantile(Erto_tot_roh$ROH_tot,c(.95))
q_Illegio <- quantile(Illegio_tot_roh$ROH_tot,c(.95))
q_Resia <- quantile(Resia_tot_roh$ROH_tot,c(.95))

# summary_FVG <- summary(FVG_tot_roh$ROH_tot)
summary_CEU <- summary(CEU_tot_roh$ROH_tot)
summary_TSI <- summary(TSI_tot_roh$ROH_tot)
summary_VBI <- summary(VBI_tot_roh$ROH_tot)
summary_CARL <- summary(CARL_tot_roh$ROH_tot)
summary_Sauris <- summary(Sauris_tot_roh$ROH_tot)
summary_Erto <- summary(Erto_tot_roh$ROH_tot)
summary_Illegio <- summary(Illegio_tot_roh$ROH_tot)
summary_Resia <- summary(Resia_tot_roh$ROH_tot)

CEU_tot_roh$pop <- "CEU"
TSI_tot_roh$pop <- "TSI"
VBI_tot_roh$pop <- "VBI"
CARL_tot_roh$pop <- "CARL"
Sauris_tot_roh$pop <- "FVG-S"
Erto_tot_roh$pop <- "FVG-E"
Illegio_tot_roh$pop <- "FVG-I"
Resia_tot_roh$pop <- "FVG-R"

write.table(CEU_tot_roh,file=paste("CEU_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(TSI_tot_roh,file=paste("TSI_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(VBI_tot_roh,file=paste("VBI_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(CARL_tot_roh,file=paste("CARL_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(Sauris_tot_roh,file=paste("Sauris_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(Erto_tot_roh,file=paste("Erto_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(Illegio_tot_roh,file=paste("Illegio_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(Resia_tot_roh,file=paste("Resia_tot_roh.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.335094 
# TSI: 95% -> 1.018848 
# FVG: 95% -> 4.769859 
# VBI: 95% -> 1.972573 

# W50 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU : 95% -> 1.253247
# TSI : 95% -> 1.147975
# FVG : 95% -> 4.582939
# VBI : 95% -> 2.059055

# W50K O7K (IBDseq, filtered at LOD = 5)
# q_CEU : 95% -> 0.4587313 
# q_TSI : 95% -> 0.9422786 
# q_FVG : 95% -> 7.134104 
# q_VBI : 95% -> 6.582102

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# q_CEU: 95% -> 1.269719
# q_TSI: 95% -> 1.138477
# q_FVG: 95% -> 4.863319
# q_VBI: 95% -> 2.069055
# q_CARL: 95% -> 2.658259

#W24 (Beagle RefinedIBD, filtered at LOD = 4)
# q_CEU: 95% -> 0.191057
# q_TSI: 95% -> 0.2552902
# q_FVG: 95% -> 1.918711
# q_VBI: 95% -> 0.7033285
# q_CARL: 95% -> 0.4669344

#W24 (Beagle RefinedIBD, filtered at LOD = 4) 01/08/2014
# q_CEU: 95% -> 4.118258 
# q_TSI: 95% -> 3.991706 
# q_FVG: 95% -> 15.47934 
# q_VBI: 95% -> 5.637363 
# q_CARL: 95% -> 7.155192 
# q_Sauris: 95% -> 13.92794 
# q_Erto: 95% -> 10.56745 
# q_Illegio: 95% -> 17.53563 
# q_Resia: 95% -> 23.23969

#W24 (Beagle RefinedIBD, filtered at LOD = 5) 01/08/2014
# q_CEU: 95% -> 3.407679 
# q_TSI: 95% -> 3.136137 
# q_FVG: 95% -> 15.98222 
# q_VBI: 95% -> 5.134627 
# q_CARL: 95% -> 6.092222 
# q_Sauris: 95% -> 13.14537 
# q_Erto: 95% -> 10.32173 
# q_Illegio: 95% -> 16.37155 
# q_Resia: 95% -> 22.91865 

#W24 (Beagle RefinedIBD, filtered at LOD = 5) 01/08/2014 GENOMEWIDE per pop
# q_CEU: 95% -> 53.80378 
# q_TSI: 95% -> 46.95559 
# q_FVG: 95% -> 143.3815 
# q_VBI: 95% -> 74.93171 
# q_CARL: 95% -> 79.10135 
# q_Sauris: 95% -> 186.39 
# q_Erto: 95% -> 174.3446 
# q_Illegio: 95% -> 210.5669 
# q_Resia: 95% -> 188.7471 

#W24 (Beagle RefinedIBD, filtered at LOD = 5) ROH 05/08/2014 GENOMEWIDE per pop (subsampled at 46 individuals)
# q_CEU: 95% -> 99.67381 
# q_TSI: 95% -> 94.80848 
# q_FVG: 95% -> 235.589 
# q_VBI: 95% -> 189.19 
# q_CARL: 95% -> 118.7458 
# q_Sauris: 95% -> 211.1938 
# q_Erto: 95% -> 163.8415 
# q_Illegio: 95% -> 209.7505 
# q_Resia: 95% -> 206.2326 

# ROH           
# POP Min.  1st_Qu. Median  Mean  3rd_Qu. Max.
# summary_CEU 62.32 77.84 83.75 85.06 89.52 169.2
# summary_TSI 62.41 74.63 79.08 80.66 84.92 120.5
# summary_FVG 24.97 125.9 150.3 152.8 174.1 294.2
# summary_VBI 14.74 85.72 94.04 110.5 114.1 404.5
# summary_CARL  30.14 71.37 83.46 83.92 93.15 126
# summary_Sauris  77.59 101 139.9 136.3 159.7 231.8
# summary_Erto  66.56 101.1 118.4 121.4 141.9 212.2
# summary_Illegio 17.23 100.9 121.3 126.1 153.8 253
# summary_Resia 19.13 127.6 147.1 147.4 172.5 229.7



pop_colors <- col_pop(pops_c)
base_folder <- getwd()



#################################################
#
roh_all_cum <- data.frame(
  ROH_tot=c(CEU_tot_roh$ROH_tot,
    TSI_tot_roh$ROH_tot,
    CARL_tot_roh$ROH_tot,
    VBI_tot_roh$ROH_tot,
    Erto_tot_roh$ROH_tot,
    Illegio_tot_roh$ROH_tot,
    Resia_tot_roh$ROH_tot,
    Sauris_tot_roh$ROH_tot),
 pop=rep(pops_c,rep(46,8)),
 ecdf=c(M_CEU(CEU_tot_roh$ROH_tot),
  M_TSI(TSI_tot_roh$ROH_tot),
  M_CARL(CARL_tot_roh$ROH_tot),
  M_VBI(VBI_tot_roh$ROH_tot),
  M_Erto(Erto_tot_roh$ROH_tot),
  M_Illegio(Illegio_tot_roh$ROH_tot),
  M_Resia(Resia_tot_roh$ROH_tot),
  M_Sauris(Sauris_tot_roh$ROH_tot)))

#manually calculate ecdf (no ecdf function used)
# ibd_all_cum <- ibd_all_cum[order(ibd_all_cum$pop),]
roh_all_cum$pop <- factor(roh_all_cum$pop,level=names(pop_colors))
# ibd_all_cum$ecdf <- ave(ibd_all_cum$IBD_tot, ibd_all_cum$pop, FUN=function(IBD_tot) seq_along(IBD_tot)/length(IBD_tot))

#plot 
pl <- ggplot(roh_all_cum)
pl <- pl + aes(x = ROH_tot, y = ecdf, colour=pop)
pl <- pl + scale_color_manual("Cohorts", values=pop_colors)
pl <- pl + geom_line(size=1.5)
pl <- pl + geom_hline(aes(yintercept=0.95), linetype=2,colour="Lightgrey",size=1.2)
pl <- pl + xlab("Total length of ROH per individual (Mb)") + ylab("Cumulative frequency")
pl <- pl + scale_x_continuous(limits=c(0,260))
pl <- pl + scale_y_continuous(limits=c(0,1))
pl <- pl + theme_bw()
pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
  
# ggsave(filename=paste(base_folder,"/test_IBD_2.jpeg",sep=""),width=8, height=8,dpi=400,plot=pl)
ggsave(filename=paste(base_folder,"/6_ROH_5POP_lod5_no_FVG.WG.ggplot.jpeg",sep=""),width=8, height=8,dpi=400,plot=pl)

#############################BOXPLOT 
require(ggplot2)
require(reshape2)
pops_c <- c("CEU","TSI","CARL","VBI","FVG-E","FVG-I","FVG-R","FVG-S")
all_tot_roh <- rbind(CEU_tot_roh,TSI_tot_roh,VBI_tot_roh,CARL_tot_roh,Sauris_tot_roh,Erto_tot_roh,Illegio_tot_roh,Resia_tot_roh)
all_tot_roh$pop <- as.factor(all_tot_roh$pop)
all_tot_roh$pop <- factor(all_tot_roh$pop, levels = pops_c)
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
pop_colors <- col_pop(pops_c)

pl <- ggplot(all_tot_roh)
pl <- pl + geom_boxplot()
pl <- pl + aes(x = factor(pop), y = ROH_tot, fill=pop)
pl <- pl + ggtitle("Total Runs of homozygosity per individual")
pl <- pl + ylab("Length (Mb)")
pl <- pl + xlab("")
# pl <- pl + guides(fill=guide_legend(title="Cohorts"))
pl <- pl + scale_fill_manual("", values=pop_colors)
pl <- pl + theme_bw()
# pl <- pl + facet_grid(cat~cons, scales="free")
pl <- pl + theme(axis.text.x=element_text(size = rel(1.8)))
pl <- pl + theme(axis.text.y=element_text(size = rel(1.8)))
pl <- pl + theme(axis.title= element_text(size=rel(1.8)))
pl <- pl + theme(legend.text= element_text(size = rel(1.8)), legend.title = element_text(size = rel(1.8)))
pl <- pl + theme(legend.position="none")

ggsave(filename=paste(getwd(),"/figure3cRev_11062015.jpeg",sep=""),width=8, height=7,dpi=400,plot=pl)


#######################ààFROM ENZA ROH BARS
require(ggplot2)
require(reshape2)

myd=read.table("all.roh.class", header=T) 
summary(myd) 
mys=subset(myd, pop!="FVG") 


# ggplot (mys, aes(x=ROH_length, fill=pop) )     

# ggplot (mys,  aes(x=ROH_length, fill=pop ) )
#  + geom_bar(stat = "bin", position ="dodge" )+ 
#  facet_grid(type ~ ., scales = "free")+
#   theme(axis.text= element_text(size=rel (1.5)))  + theme_bw() +   
#   scale_fill_manual(values = c("#F54E4E", "#13256F", "#10DFBC",  "#48867F", "#16C224", "#00631F", "#619FE0", "#CA8A1A"))

pl <- ggplot(all_pop_all_MAF_table_syn_reshaped)
pl <- pl + geom_bar(stat="identity",width=0.81, position = position_dodge(width=0.8),colour="black")
pl <- pl + aes(x = factor(breaks), y = value, fill=variable)
pl <- pl + xlab("MAF")
# pl <- pl + ylab("Proportion of sites (%)")
pl <- pl + ylab("Site count")
pl <- pl + scale_fill_manual("", values=all_cols)
pl <- pl + facet_wrap( ~ cat, ncol=1)
pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
pl <- pl + theme_bw()

pl <- pl + theme(axis.text.x=element_text(size = rel(1.2),angle=90, vjust=0.5))
pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
ggsave(filename=paste(base_folder,"/ALL_MAF_pop_SYNONYMOUS.jpeg",sep=""),width=18, height=7,dpi=300,plot=pl)


pl <- ggplot(mys)
pl <- pl + geom_bar(stat = "bin", position = "dodge",colour="black")
pl <- pl + aes(x = ROH_length, fill=pop)
# pl <- pl + ggtitle("Total Runs of homozygosity per individual")
pl <- pl + ylab("Count")
pl <- pl + xlab("ROH length")
# pl <- pl + guides(fill=guide_legend(title="Cohorts"))
pl <- pl + scale_fill_manual(values = c("#F54E4E", "#13256F", "#10DFBC",  "#48867F", "#16C224", "#00631F", "#619FE0", "#CA8A1A"))
pl <- pl + facet_grid(type ~ ., scales = "free")
pl <- pl + theme_bw()
pl <- pl + theme(axis.text.x=element_text(size = rel(1.6)))
pl <- pl + theme(axis.text.y=element_text(size = rel(1.6)))
pl <- pl + theme(axis.title= element_text(size=rel(1.6)))
pl <- pl + theme(legend.text= element_text(size = rel(1.6)), legend.title = element_text(size = rel(1.6)))
ggsave(filename=paste(getwd(),"/figureS7.jpeg",sep=""),width=14, height=7,dpi=300,plot=pl)

# pl <- pl + theme(legend.position="none")

