#Script to plot stuff for quality check of INGI ref panels
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/assign_bins.r")
require(ggplot2)
library(plyr)

# define populations and reference panes
# pops <- c("FVG")
# pops <- c("FVG","INCIPE2","CARL")
# pops <- c("INCIPE2")
# pops <- c("CARL")
#last mod 4/4/2016
pops <- c("CARL","FVG","INCIPE2","VBI")
# pops <- c("FVG")
# ingi_panels <- c("CARL.shapeit","FVG.shapeit","VBI.shapeit")
# gen_pop_ref_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit", "EUR.shapeit")
# all_panels <- c(ingi_panels,gen_pop_ref_panels)
#for each population, we upload the interesting columns of the info file:
# base_folder <- "/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/IMPUTED"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/31032016_IMPUTATION"
#rsId, position,expected_af, info, type and concordances 

# modes <- c("3BIN","11BIN")
# for (mode in modes) {
# print(mode)

for (pop in pops){
    # pop <- "INCIPE2"
    # pop <- "CARL"
    # pop <- "FVG"
    print(pop)
    current_pop_all_panels_all_chr <- NULL

    if (pop == "CARL"){
        # selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","CARL.shapeit")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","CARL.shapeit")
    } else if (pop == "VBI"){
        # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","VBI.shapeit")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","VBI.shapeit")
    }else if (pop == "FVG"){
        # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","FVG.shapeit")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","FVG.shapeit")
    }else if (pop == "INCIPE2"){
        # selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","1000GP_Phase3.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit")
    }

    for(panel in selected_panels){
        # panel <- "CARL.shapeit"
        print(panel)
        for (chr in 21:21){
            # chr <- 1
            pop_folder <- paste(base_folder,pop,panel,sep="/")
            current_pop_current_panel_current_chr_info_name <- paste(pop_folder,"/",chr,"/chr",chr,".gen_info_partial_t2.gz",sep="")
            # current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,nrows=100000)
            # CHROM RS_ID POS EXP_FREQ_A1 INFO TYPE INFO_TYPE0 CONCORD_TYPE0 r2_TYPE0 COHORT PANEL MAF BIN
            current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,sep=" ",stringsAsFactors=F, comment.char="",colClasses=c("integer","character","integer",rep("numeric",2),"character","NULL","NULL","NULL",rep("character",2),rep("numeric",2)))
            current_pop_all_panels_all_chr <- rbind(current_pop_all_panels_all_chr,current_pop_current_panel_current_chr_info)
            # assign(paste("complete_",pop,"info",sep=""),current_pop_all_panels_all_chr)

        }
    }
    
    current_pop_all_panels_all_chr$BIN3 <- 0
    # maf_bins <- c(0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50)
    # maf_bins <- c(0,0.01,0.02,0.05,0.5)
    maf_bins <- c(0,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)
    mode <- paste(length(maf_bins),"BIN",sep="")
    for (i in 1:(length(maf_bins))){
        if (i == 1){
            print(i)
            try(current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$MAF <= maf_bins[i]),]$BIN3 <- i-1,silent=T)
        } else {
            print(i)
            try(current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$MAF <= maf_bins[i] & current_pop_all_panels_all_chr$MAF > maf_bins[i-1]),]$BIN3 <- i-1,silent=T)

        }
    }

    # remove all monomorphic data plus all the data with maf higher than the last bin
    current_pop_all_panels_all_chr_nomono <- current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$BIN3 != 0),]

    cdata <- ddply(current_pop_all_panels_all_chr_nomono, c("BIN3","PANEL"), summarise,
                   N    = length(INFO),
                   mean = mean(INFO),
                   median = median(INFO),
                   sd   = sd(INFO),
                   min   = min(INFO),
                   max   = max(INFO),
                   p25   = quantile(INFO,c(0.25)),
                   p75   = quantile(INFO,c(0.75)),
                   se   = sd / sqrt(N),
                   variance= var(INFO)
    )

    cdata_nomono <- cdata[which(cdata$BIN3 != 0),]
    cdata_nomono$BIN <- 0

    for (i in 1:(length(maf_bins))){
            print(i)
            try(cdata_nomono[which(cdata_nomono$BIN3 == i-1),]$BIN <- maf_bins[i],silent=T)
    }

    cdata_nomono$PANNELLO <- cdata_nomono$PANEL
    cdata_nomono$PANNELLO <- factor(cdata_nomono$PANNELLO,selected_panels)
    cdata_nomono$BIN <- factor(cdata_nomono$BIN,maf_bins)
   #create folder for each population
   out_folder <- paste(base_folder,"/PLOTS/",pop,"_",mode,sep="")
   dir.create(out_folder, recursive=T)

    # save(cdata_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_nomono.RData",sep=""))
    # save(cdata,file=paste(out_folder,"/",pop,"_",mode,"_cdata.RData",sep=""))
    write.table(cdata,file=paste(out_folder,"/",pop,"_",mode,"_cdata.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
    write.table(cdata_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_nomono.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
    # # load(paste(pop,"_cdata_nomono.RData",sep=""))

    #plot without error bars of the info MEDIAN
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_nomono)
    pl <- pl + aes(x = BIN3, y = median, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2,colour="black")
    # pl <- pl + geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), width=0.2)
    pl <- pl + scale_x_continuous(breaks=cdata_nomono$BIN3,labels=cdata_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Median IMPUTE info score")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    # ggsave(filename=paste(getwd(),"/complete_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
    ggsave(filename=paste(out_folder,"/complete_13042016_",pop,"_",mode,"_info_nobar.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    #plot without error bars of the info AVERAGE
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_nomono)
    pl <- pl + aes(x = BIN3, y = mean, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2,colour="black")
    # pl <- pl + geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), width=0.2)
    pl <- pl + scale_x_continuous(breaks=cdata_nomono$BIN3,labels=cdata_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Average IMPUTE info score")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    # pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = "top")
    # ggsave(filename=paste(getwd(),"/complete_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
    ggsave(filename=paste(out_folder,"/complete_13042016_",pop,"_",mode,"_info_nobar_mean.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)


    #plot without error bars of the info VARIANCE
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_nomono)
    pl <- pl + aes(x = BIN3, y = variance, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2,colour="black")
    # pl <- pl + geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), width=0.2)
    pl <- pl + scale_x_continuous(breaks=cdata_nomono$BIN3,labels=cdata_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Variance of the IMPUTE info score")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = "top")
    # ggsave(filename=paste(getwd(),"/complete_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
    ggsave(filename=paste(out_folder,"/complete_13042016_",pop,"_",mode,"_info_variance.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    # cdata_nomono$LAB <- paste("BIN:",cdata_nomono$BIN,"(N:",cdata_nomono$N,")",sep="")

    #plot also in boxplot form with median data
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_nomono)
    pl <- pl + aes(x = PANNELLO, y = mean)
    # pl <- pl + geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2,colour="black")
    # pl <- pl + geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), width=0.2)
    # pl <- pl + scale_x_continuous(breaks=cdata_nomono$BIN2)
    pl <- pl + geom_boxplot(aes(lower=p25,middle=median,upper=p75,ymax=max,ymin=min,fill=PANNELLO),stat = "identity")
    pl <- pl + facet_wrap( ~ BIN,scales = "free_y")
    # pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("Panels") + ylab("Average IMPUTE info score")
    pl <- pl + theme_bw()
    # pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.x=element_blank())
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    # pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
    # ggsave(filename=paste(getwd(),"/complete_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
    ggsave(filename=paste(out_folder,"/complete_13042016_",pop,"_",mode,"_info_boxplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

}
 
##############################################
# Plot R2 values by maf an cocnord
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/assign_bins.r")
require(ggplot2)
library(plyr)

# define populations and reference panes
#last mod 4/4/2016
# pops <- c("CARL","FVG","INCIPE2","VBI")
pops <- c("CARL","INCIPE2","VBI")
# pops <- c("FVG")
# ingi_panels <- c("CARL.shapeit","FVG.shapeit","VBI.shapeit")
# gen_pop_ref_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit", "EUR.shapeit")
# all_panels <- c(ingi_panels,gen_pop_ref_panels)
#for each population, we upload the interesting columns of the info file:
# base_folder <- "/lustre/scratch114/teams/soranzo/users/mc14/fromscratch113/INGI/05272015_MERGED_REF_PANEL/IMPUTED"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/31032016_IMPUTATION"
#rsId, position,expected_af, info, type and concordances 
# maf_bins <- c(0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50)
# modes <- c("3BIN","11BIN")
# for (mode in modes) {
# print(mode)
# mode <- "3BIN"

for (pop in pops){
    # pop <- "INCIPE2"
    # pop <- "CARL"
    # pop <- "FVG"
    print(pop)
    current_pop_all_panels_all_chr <- NULL

    if (pop == "CARL"){
        # selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","CARL.shapeit")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","CARL.shapeit")
    } else if (pop == "VBI"){
        # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","VBI.shapeit")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","VBI.shapeit")
    }else if (pop == "FVG"){
        # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","FVG.shapeit")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","FVG.shapeit")
    }else if (pop == "INCIPE2"){
        # selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","1000GP_Phase3.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref")
        selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit")
    }

    for(panel in selected_panels){
        # panel <- "1000Gph1.shapeit"
        print(panel)
        for (chr in 21:21){
            # chr <- 1
            pop_folder <- paste(base_folder,pop,panel,sep="/")
            current_pop_current_panel_current_chr_info_name <- paste(pop_folder,"/",chr,"/chr",chr,".gen_info_partial_t2.gz",sep="")
            # current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,nrows=100000)
            # CHROM RS_ID POS EXP_FREQ_A1 INFO TYPE INFO_TYPE0 CONCORD_TYPE0 r2_TYPE0 COHORT PANEL MAF BIN
            current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,sep=" ",stringsAsFactors=F, comment.char="",colClasses=c("integer","character","integer",rep("numeric",2),"character",rep("numeric",3),rep("character",2),rep("numeric",2)))
            current_pop_all_panels_all_chr <- rbind(current_pop_all_panels_all_chr,current_pop_current_panel_current_chr_info)
            # assign(paste("complete_",pop,"info",sep=""),current_pop_all_panels_all_chr)

        }
    }
    
    #select only type 2 snps
    current_pop_all_panels_all_chr <- current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$TYPE == "2"),]

    current_pop_all_panels_all_chr$BIN3 <- 0
    # maf_bins <- c(0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50)
    # maf_bins <- c(0,0.01,0.02,0.05,0.5)
    maf_bins <- c(0,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)
    mode <- paste(length(maf_bins),"BIN",sep="")
    for (i in 1:(length(maf_bins))){
        if (i == 1){
            print(i)
            try(current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$MAF <= maf_bins[i]),]$BIN3 <- i-1,silent=T)
        } else {
            print(i)
            try(current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$MAF <= maf_bins[i] & current_pop_all_panels_all_chr$MAF > maf_bins[i-1]),]$BIN3 <- i-1,silent=T)

        }
    }

    current_pop_all_panels_all_chr_nomono <- current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$BIN3 != 0),]
    current_pop_all_panels_all_chr_nomono <- current_pop_all_panels_all_chr_nomono[which(current_pop_all_panels_all_chr_nomono$r2_TYPE0 != -1),]

    cdata_r2 <- ddply(current_pop_all_panels_all_chr_nomono, c("BIN3","PANEL"), summarise,
                   N    = length(r2_TYPE0),
                   mean = mean(r2_TYPE0),
                   sd   = sd(r2_TYPE0),
                   se   = sd / sqrt(N),
                   median = median(r2_TYPE0),
                   min   = min(r2_TYPE0),
                   max   = max(r2_TYPE0),
                   p25   = quantile(r2_TYPE0,c(0.25)),
                   p75   = quantile(r2_TYPE0,c(0.75)),
                   variance= var(r2_TYPE0)
    )
    cdata_r2_nomono <- cdata_r2[which(cdata_r2$BIN3 != 0),]
    # cdata_r2_nomono <- cdata_r2_nomono[which(!is.na(cdata_r2_nomono$se)),]
    cdata_conc <- ddply(current_pop_all_panels_all_chr_nomono, c("BIN3","PANEL"), summarise,
                   N    = length(CONCORD_TYPE0),
                   mean = mean(CONCORD_TYPE0),
                   sd   = sd(CONCORD_TYPE0),
                   se   = sd / sqrt(N),
                   median = median(CONCORD_TYPE0),
                   min   = min(CONCORD_TYPE0),
                   max   = max(CONCORD_TYPE0),
                   p25   = quantile(CONCORD_TYPE0,c(0.25)),
                   p75   = quantile(CONCORD_TYPE0,c(0.75)),
                   variance= var(CONCORD_TYPE0)
    )
    cdata_conc_nomono <- cdata_conc[which(cdata_conc$BIN3 != 0),]
    # cdata_conc_nomono <- cdata_conc_nomono[which(!is.na(cdata_conc_nomono$se)),]
    #plot the same cohort stratifying by panel
    # The errorbars overlapped, so use position_dodge to move them horizontally
    cdata_r2_nomono$BIN <- 0

    for (i in 1:(length(maf_bins))){
            print(i)
            try(cdata_r2_nomono[which(cdata_r2_nomono$BIN3 == i-1),]$BIN <- maf_bins[i],silent=T)
    }

    cdata_r2_nomono$BIN <- factor(cdata_r2_nomono$BIN,maf_bins)
    cdata_r2_nomono$PANNELLO <- cdata_r2_nomono$PANEL
    cdata_r2_nomono$PANNELLO <- factor(cdata_r2_nomono$PANNELLO,selected_panels)

    cdata_conc_nomono$BIN <- 0
    for (i in 1:(length(maf_bins))){
            print(i)
            try(cdata_conc_nomono[which(cdata_conc_nomono$BIN3 == i-1),]$BIN <- maf_bins[i],silent=T)
    }

    cdata_conc_nomono$BIN <- factor(cdata_conc_nomono$BIN,maf_bins)

    cdata_conc_nomono$PANNELLO <- cdata_conc_nomono$PANEL
    cdata_conc_nomono$PANNELLO <- factor(cdata_conc_nomono$PANNELLO,selected_panels)
   
   #create folder for each population
   out_folder <- paste(base_folder,"/PLOTS/",pop,"_",mode,sep="")
   dir.create(out_folder, recursive=T)

    write.table(cdata_r2,file=paste(out_folder,"/",pop,"_",mode,"_cdata_r2.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
    write.table(cdata_conc,file=paste(out_folder,"/",pop,"_",mode,"_cdata_conc.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
    write.table(cdata_r2_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_r2_nomono.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
    write.table(cdata_conc_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_conc_nomono.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
    
    # plot mean r2
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_r2_nomono)
    pl <- pl + aes(x = BIN3, y = mean, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_r2_nomono$BIN3,labels=cdata_r2_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Average IMPUTE r2")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    ggsave(filename=paste(out_folder,"/complete_13042016_r2_",pop,"_",mode,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    # plot variance
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_r2_nomono)
    pl <- pl + aes(x = BIN3, y = variance, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_r2_nomono$BIN3,labels=cdata_r2_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Variance of IMPUTE r2")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.2, 0.2))
    ggsave(filename=paste(out_folder,"/complete_13042016_r2_",pop,"_",mode,"_variance.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    # plot median
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_r2_nomono)
    pl <- pl + aes(x = BIN3, y = median, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_r2_nomono$BIN3,labels=cdata_r2_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Median of IMPUTE r2")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    ggsave(filename=paste(out_folder,"/complete_13042016_r2_",pop,"_",mode,"_median.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    #plot also in boxplot form with median data
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_r2_nomono)
    pl <- pl + aes(x = PANNELLO, y = mean)
    pl <- pl + geom_boxplot(aes(lower=p25,middle=median,upper=p75,ymax=max,ymin=min,fill=PANNELLO),stat = "identity")
    pl <- pl + facet_wrap( ~ BIN,scales = "free_y")
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("Panels") + ylab("Average IMPUTE r2 ")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_blank())
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    ggsave(filename=paste(out_folder,"/complete_13042016_r2_",pop,"_",mode,"_info_boxplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    #plot mean data for concordance
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_conc_nomono)
    pl <- pl + aes(x = BIN3, y = mean, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_conc_nomono$BIN3,labels=cdata_conc_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Average CONCORDANCE")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    ggsave(filename=paste(out_folder,"/complete_13042016_conc_",pop,"_",mode,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    #plot variance for concordance
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_conc_nomono)
    pl <- pl + aes(x = BIN3, y = variance, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_conc_nomono$BIN3,labels=cdata_conc_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Variance of CONCORDANCE")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = "top")
    ggsave(filename=paste(out_folder,"/complete_13042016_conc_",pop,"_",mode,"_variance.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    #plot median for concordance
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_conc_nomono)
    pl <- pl + aes(x = BIN3, y = median, colour=PANNELLO,group=PANNELLO)
    # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_conc_nomono$BIN3,labels=cdata_conc_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
    # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Median of CONCORDANCE")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.2, 0.2))
    ggsave(filename=paste(out_folder,"/complete_13042016_conc_",pop,"_",mode,"_median.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

    #plot also in boxplot form with median data
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_r2_nomono)
    pl <- pl + aes(x = PANNELLO, y = mean)
    pl <- pl + geom_boxplot(aes(lower=p25,middle=median,upper=p75,ymax=max,ymin=min,fill=PANNELLO),stat = "identity")
    pl <- pl + facet_wrap( ~ BIN,scales = "free_y")
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("Panels") + ylab(" IMPUTE CONCORDANCE ")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_blank())
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    ggsave(filename=paste(out_folder,"/complete_13042016_conc_",pop,"_",mode,"_info_boxplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
   
}

####################################################################
