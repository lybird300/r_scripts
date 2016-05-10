##############################################
# Plot R2 values by maf an cocnord
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/assign_bins.r")
require(ggplot2)
library(plyr)
args=commandArgs(trailing=TRUE)
maf_bins <- args
# define populations and reference panes
#last mod 4/4/2016
# pops <- c("CARL","FVG","INCIPE2","VBI")
pops <- c("CARL","INCIPE2","VBI")
# pops <- c("FVG")
# define panel set for plotting purpouses
panel_set1 <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","EUR.shapeit")
panel_set2 <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","TGP3_ALL.shapeit","EUR.shapeit")
# panel_set3 <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","EUR.shapeit")
panel_set3 <- c("CARL_FVG_VBI.shapeit","uk10k1kg.ref","EUR.shapeit")
all_set <- c("panel_set1","panel_set2","panel_set3")

current_date <- format(Sys.time(),"%d_%m_%Y_%H%M%S")
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
    for (pan_set in all_set){
        current_pan_set <- get(pan_set)
        current_pop_all_panels_all_chr <- NULL

        if (pop == "CARL"){
            # selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","CARL.shapeit")
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","CARL.shapeit")
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","CARL.shapeit")
            selected_panels <- c(current_pan_set,"CARL.shapeit")
            } else if (pop == "VBI"){
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","VBI.shapeit")
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","VBI.shapeit")
            selected_panels <- c(current_pan_set,"VBI.shapeit")
            }else if (pop == "FVG"){
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","FVG.shapeit")
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit","FVG.shapeit")
            selected_panels <- c(current_pan_set,"FVG.shapeit")
            }else if (pop == "INCIPE2"){
            # selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","1000GP_Phase3.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref")
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit")
            # selected_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit","EUR.shapeit")
            selected_panels <- c(current_pan_set)
        }

        for(panel in selected_panels){
            # panel <- "1000Gph1.shapeit"
            print(panel)
            for (chr in 2:2){
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
        # maf_bins <- c(0,0.01,0.02,0.05,0.1,0.2,0.5)
        # maf_bins <- c(0,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)
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
        out_folder <- paste(base_folder,"/PLOTS/",current_date,"/",pan_set,"/",pop,"_",mode,sep="")
        dir.create(out_folder, recursive=T)

        write.table(cdata_r2,file=paste(out_folder,"/",pop,"_",mode,"_cdata_r2.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
        write.table(cdata_conc,file=paste(out_folder,"/",pop,"_",mode,"_cdata_conc.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
        write.table(cdata_r2_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_r2_nomono.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
        write.table(cdata_conc_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_conc_nomono.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

        #######################################################
        #plot without error bars of the r2 on MEDIAN MEAN and VARIANCE
        plot_range <- c("median","mean","variance")

        for (variable in plot_range){
            print(variable)
            # plot mean r2
            pd <- position_dodge(0.001)
            pl <- ggplot(cdata_r2_nomono)
            pl <- pl + aes(x = BIN3, y = cdata_r2_nomono[,variable], colour=PANNELLO,group=PANNELLO)
            # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
            pl <- pl + scale_x_continuous(breaks=cdata_r2_nomono$BIN3,labels=cdata_r2_nomono$BIN)
            pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
            # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
            pl <- pl + scale_color_discrete(name=paste(pop," on:"))
            pl <- pl + xlab("MAF bins") + ylab(paste(toupper(variable),"IMPUTE r2",sep=" "))
            pl <- pl + theme_bw()
            pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
            pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
            pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
            pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
            ggsave(filename=paste(out_folder,"/complete_13042016_r2_",pop,"_",mode,"_info_",toupper(variable),".jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

        } 

        #plot also in boxplot form with median data
        pd <- position_dodge(0.001)
        pl <- ggplot(cdata_r2_nomono)
        pl <- pl + aes(x = PANNELLO, y = mean)
        pl <- pl + geom_boxplot(aes(lower=p25,middle=median,upper=p75,ymax=max,ymin=min,fill=PANNELLO),stat = "identity")
        pl <- pl + facet_wrap( ~ BIN,scales = "free_y")
        pl <- pl + scale_color_discrete(name=paste(pop," on:"))
        pl <- pl + xlab("Panels") + ylab("IMPUTE r2 ")
        pl <- pl + theme_bw()
        pl <- pl + theme(axis.text.x=element_blank())
        pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
        pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
        ggsave(filename=paste(out_folder,"/complete_13042016_r2_",pop,"_",mode,"_info_boxplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
        
        #######################################################
        #plot without error bars of the concordance on MEDIAN MEAN and VARIANCE
        
        plot_range <- c("median","mean","variance")

        for (variable in plot_range){
            print(variable)
            #plot mean data for concordance
            pd <- position_dodge(0.001)
            pl <- ggplot(cdata_conc_nomono)
            pl <- pl + aes(x = BIN3, y = cdata_conc_nomono[,variable], colour=PANNELLO,group=PANNELLO)
            # pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
            pl <- pl + scale_x_continuous(breaks=cdata_conc_nomono$BIN3,labels=cdata_conc_nomono$BIN)
            pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
            # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
            pl <- pl + scale_color_discrete(name=paste(pop," on:"))
            pl <- pl + xlab("MAF bins") + ylab(paste(toupper(variable),"IMPUTE Concordance",sep=" "))
            pl <- pl + theme_bw()
            pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
            pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
            pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
            pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
            ggsave(filename=paste(out_folder,"/complete_13042016_conc_",pop,"_",mode,"_info_",toupper(variable),".jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
        }
           
        #plot also in boxplot form with median data
        pd <- position_dodge(0.001)
        pl <- ggplot(cdata_r2_nomono)
        pl <- pl + aes(x = PANNELLO, y = mean)
        pl <- pl + geom_boxplot(aes(lower=p25,middle=median,upper=p75,ymax=max,ymin=min,fill=PANNELLO),stat = "identity")
        pl <- pl + facet_wrap( ~ BIN,scales = "free_y")
        pl <- pl + scale_color_discrete(name=paste(pop," on:"))
        pl <- pl + xlab("Panels") + ylab("IMPUTE CONCORDANCE ")
        pl <- pl + theme_bw()
        pl <- pl + theme(axis.text.x=element_blank())
        pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
        pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
        ggsave(filename=paste(out_folder,"/complete_13042016_conc_",pop,"_",mode,"_info_boxplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
       
    }
}

############