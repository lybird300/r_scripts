#r script to plot for enza
rm(list=ls())
###########################################################################################
#Plot 6 ROH
chr <- "10"
pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")

# input_format <- "PLINK"
input_format <- "BEAGLE"
xmax <- NULL
# LOD <- 4
LOD <- 5
# base_folder <- getwd() #only if we're lazy
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    # pop_roh_file <- paste(pop,"roh.length.hbd",sep=".")
    # pop_roh_file <- paste(pop,"roh.length.4.hbd",sep=".")
    pop_roh_file <- paste("CHR",chr,"/",pop,".roh.length.",LOD,".hbd",sep="")
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
  # jpeg(paste("6_roh_",pop".jpg",sep=""),width=800, height=800)
  # plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity")
  # lines(M_FVG,FVG_tot_roh$ROH_tot,col="red")
  # lines(M_TSI,TSI_tot_roh$ROH_tot,col="green")
  # lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue")
  # legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
  # dev.off()

}

q_CEU <- quantile(CEU_tot_roh$ROH_tot,c(.95))
q_TSI <- quantile(TSI_tot_roh$ROH_tot,c(.95))
q_FVG <- quantile(FVG_tot_roh$ROH_tot,c(.95))
q_VBI <- quantile(VBI_tot_roh$ROH_tot,c(.95))
q_CARL <- quantile(CARL_tot_roh$ROH_tot,c(.95))
q_Sauris <- quantile(Sauris_tot_roh$ROH_tot,c(.95))
q_Erto <- quantile(Erto_tot_roh$ROH_tot,c(.95))
q_Illegio <- quantile(Illegio_tot_roh$ROH_tot,c(.95))
q_Resia <- quantile(Resia_tot_roh$ROH_tot,c(.95))

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

# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all_IBDseq.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all.jpg",sep=""),width=800, height=800)
base_folder <- getwd()

jpeg(paste(base_folder,"/6_roh_all_5POP",chr,".jpg",sep=""),width=1000, height=1000)
jpeg(paste(base_folder,"/6_roh_all_5POP_lod5_",chr,".jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=1.5)
  plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  lines(M_TSI,TSI_tot_roh$ROH_tot,col="green", verticals=TRUE, pch=46)
  lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue", verticals=TRUE, pch=46)
  lines(M_FVG,FVG_tot_roh$ROH_tot,col="red", verticals=TRUE, pch=46)
  lines(M_CARL,CARL_tot_roh$ROH_tot,col="yellow", verticals=TRUE, pch=46)
  lines(M_Erto,Erto_tot_roh$ROH_tot,col="darkred", verticals=TRUE, pch=46)
  lines(M_Illegio,Illegio_tot_roh$ROH_tot,col="indianred", verticals=TRUE, pch=46)
  lines(M_Resia,Resia_tot_roh$ROH_tot,col="mediumvioletred", verticals=TRUE, pch=46)
  lines(M_Sauris,Sauris_tot_roh$ROH_tot,col="orangered", verticals=TRUE, pch=46)
  legend("bottomright",pch =c(rep(19,length(pops))),legend=pops,col=c("black","green","blue","red","yellow","darkred","indianred","mediumvioletred","orangered"),ncol=length(pops))
dev.off()

# jpeg(paste(base_folder,"PLOTS/6_roh_iso_IBDseq.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_iso.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"/6_roh_iso_3POP.jpg",sep=""),width=800, height=800)
  par(lwd=4,cex=2)
  plot(M_FVG,FVG_tot_roh$ROH_tot,verticals=TRUE, pch=46,main="",xlab="Total ROH homozigosity (Mb)",col="red")
  lines(M_VBI,VBI_tot_roh$ROH_tot,verticals=TRUE, pch=46,col="blue")
  lines(M_CARL,CARL_tot_roh$ROH_tot,verticals=TRUE, pch=46,col="yellow")
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("FVG","VBI","CARL"),col=c("red","blue","yellow"),ncol=3)
dev.off()

# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/7_roh.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"PLOTS/7_roh_density.jpg",sep=""),width=800, height=800)
  plot(density(VBI_tot_roh$ROH_tot),main="",xlab="Total ROH homozigosity", col="blue")
  lines(density(FVG_tot_roh$ROH_tot),col="red")
  lines(density(TSI_tot_roh$ROH_tot),col="green")
  lines(density(CEU_tot_roh$ROH_tot))
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()

###########################################################################################
#Plot 6_A_ROH GENOMEWIDE BY POPULATION
###########################################################################################

pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")

# input_format <- "PLINK"
input_format <- "BEAGLE"
xmax <- NULL
# LOD <- 4
LOD <- 5
# base_folder <- getwd() #only if we're lazy
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    # pop_roh_file <- paste(pop,"roh.length.hbd",sep=".")
    # pop_roh_file <- paste(pop,"roh.length.4.hbd",sep=".")
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
  # jpeg(paste("6_roh_",pop".jpg",sep=""),width=800, height=800)
  # plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity")
  # lines(M_FVG,FVG_tot_roh$ROH_tot,col="red")
  # lines(M_TSI,TSI_tot_roh$ROH_tot,col="green")
  # lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue")
  # legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
  # dev.off()

}

q_CEU <- quantile(CEU_tot_roh$ROH_tot,c(.95))
q_TSI <- quantile(TSI_tot_roh$ROH_tot,c(.95))
q_FVG <- quantile(FVG_tot_roh$ROH_tot,c(.95))
q_VBI <- quantile(VBI_tot_roh$ROH_tot,c(.95))
q_CARL <- quantile(CARL_tot_roh$ROH_tot,c(.95))
q_Sauris <- quantile(Sauris_tot_roh$ROH_tot,c(.95))
q_Erto <- quantile(Erto_tot_roh$ROH_tot,c(.95))
q_Illegio <- quantile(Illegio_tot_roh$ROH_tot,c(.95))
q_Resia <- quantile(Resia_tot_roh$ROH_tot,c(.95))

summary_CEU <- summary(CEU_tot_roh$ROH_tot)
summary_TSI <- summary(TSI_tot_roh$ROH_tot)
summary_FVG <- summary(FVG_tot_roh$ROH_tot)
summary_VBI <- summary(VBI_tot_roh$ROH_tot)
summary_CARL <- summary(CARL_tot_roh$ROH_tot)
summary_Sauris <- summary(Sauris_tot_roh$ROH_tot)
summary_Erto <- summary(Erto_tot_roh$ROH_tot)
summary_Illegio <- summary(Illegio_tot_roh$ROH_tot)
summary_Resia <- summary(Resia_tot_roh$ROH_tot)

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



# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all_IBDseq.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all.jpg",sep=""),width=800, height=800)
pop_colors <- col_pop(pops)
base_folder <- getwd()

# jpeg(paste(base_folder,"/6_roh_all_5POP",chr,".jpg",sep=""),width=1000, height=1000)
jpeg(paste(base_folder,"/6_A_roh_all_5POP_lod5.WG.jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=1.5)
  # plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  plot(M_CEU,CEU_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "CEU"),1],main="",xlab="Total ROH homozigosity (Mb)", xlim=c(0,450), verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_TSI,TSI_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "TSI"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_VBI,VBI_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "VBI"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_FVG,FVG_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "FVG"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_CARL,CARL_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "CARL"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Erto,Erto_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Erto"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Illegio,Illegio_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Illegio"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Resia,Resia_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Resia"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Sauris,Sauris_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Sauris"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  abline(h=0.95,col='grey',lty='dashed')
  #define parameters for legend
  
  leg_txt <- c(pop_colors[which(pop_colors$pop == "CEU"),2],pop_colors[which(pop_colors$pop == "TSI"),2],pop_colors[which(pop_colors$pop == "VBI"),2],pop_colors[which(pop_colors$pop == "FVG"),2],pop_colors[which(pop_colors$pop == "CARL"),2],pop_colors[which(pop_colors$pop == "Erto"),2],pop_colors[which(pop_colors$pop == "Illegio"),2],pop_colors[which(pop_colors$pop == "Resia"),2],pop_colors[which(pop_colors$pop == "Sauris"),2])
  bkg <- c(pop_colors[which(pop_colors$pop == "CEU"),1],pop_colors[which(pop_colors$pop == "TSI"),1],pop_colors[which(pop_colors$pop == "VBI"),1],pop_colors[which(pop_colors$pop == "FVG"),1],pop_colors[which(pop_colors$pop == "CARL"),1],pop_colors[which(pop_colors$pop == "Erto"),1],pop_colors[which(pop_colors$pop == "Illegio"),1],pop_colors[which(pop_colors$pop == "Resia"),1],pop_colors[which(pop_colors$pop == "Sauris"),1])

  legend("bottomright",pch =c(rep(22,length(pops))),legend=leg_txt, pt.lwd=2,pt.cex=2,pt.bg=bkg,col=c(rep('black',length(pops))),ncol=4,bty="n")
dev.off()
