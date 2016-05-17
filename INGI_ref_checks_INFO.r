#Script to plot stuff for quality check of INGI ref panels
# from INFO SCORE

rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/imputation_test.r")
require(ggplot2)
library(plyr)
args=commandArgs(trailing=TRUE)
maf_bins <- args
# maf_bins <- c(0,0.01,0.02,0.05,0.1,0.2,0.5)
# commandArgs <- function() maf_bins

# define populations and reference panes
#last mod 4/4/2016
pops <- c("CARL","FVG","INCIPE2","VBI")
# pops <- c("CARL")
panel_set1 <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","EUR.shapeit")
panel_set2 <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","TGP3_ALL.shapeit","EUR.shapeit")
panel_set3 <- c("CARL_FVG_VBI.shapeit","uk10k1kg.ref","EUR.shapeit")
panel_set4 <- c("CARL_FVG_VBI.shapeit","uk10k1kg.ref","CARL_FVG_VBI_TGP3_ALL.shapeit","TGP3_ALL.shapeit")
# panel_set3 <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_UK10K_TGP3_ALL.shapeit","uk10k1kg.ref","EUR.shapeit")
all_set <- c("panel_set1","panel_set2","panel_set3","panel_set4")
# gen_pop_ref_panels <- c("CARL_FVG_VBI.shapeit","CARL_FVG_VBI_TSI.shapeit","CARL_FVG_VBI_TGP3_ALL.shapeit","uk10k1kg.ref","TGP3_ALL.shapeit", "EUR.shapeit")
#for each population, we upload the interesting columns of the info file:
current_date <- format(Sys.time(),"%d_%m_%Y_%H%M%S")
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/31032016_IMPUTATION"

for (pop in pops){
    # pop <- "CARL"
    #first we need to read ALL panels, based on the population
    if (pop == "CARL"){
        selected_panels <- c(panel_set1,panel_set2,panel_set3,panel_set4,"CARL.shapeit")
    } else if (pop == "VBI"){
        selected_panels <- c(panel_set1,panel_set2,panel_set3,panel_set4,"VBI.shapeit")
    }else if (pop == "FVG"){
        selected_panels <- c(panel_set1,panel_set2,panel_set3,panel_set4,"FVG.shapeit")
    }else if (pop == "INCIPE2"){
        selected_panels <- c(panel_set1,panel_set2,panel_set3,panel_set4)
    }

    all_panels <- unique(selected_panels)
    current_pop_all_panels_all_chr_complete <- NULL
    #read all panels and save the data ONCE
    for(panel in all_panels){
        # panel <- "CARL.shapeit"
        print(panel)
        for (chr in c(2,21)){
            # chr <- 21
            pop_folder <- paste(base_folder,pop,panel,sep="/")
            current_pop_current_panel_current_chr_info_name <- paste(pop_folder,"/",chr,"/chr",chr,".gen_info_partial_t2.gz",sep="")
            # current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,nrows=100000)
            # CHROM RS_ID POS EXP_FREQ_A1 INFO TYPE INFO_TYPE0 CONCORD_TYPE0 r2_TYPE0 COHORT PANEL MAF BIN
            current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,sep=" ",stringsAsFactors=F, comment.char="",colClasses=c("integer","character","integer",rep("numeric",2),"character","NULL","NULL","NULL",rep("character",2),rep("numeric",2)))
            current_pop_all_panels_all_chr_complete <- rbind(current_pop_all_panels_all_chr_complete,current_pop_current_panel_current_chr_info)
            # assign(paste("complete_",pop,"info",sep=""),current_pop_all_panels_all_chr)
        }
    }

    #now, for each panel set we'll extract the relevant data from the WHOLE current_pop_all_panels_all_chr
    print(pop)
    for (pan_set in all_set){
        # pan_set <- (all_set[1])
        current_pan_set <- get(pan_set)


        if (pop == "CARL"){
            selected_panels <- c(current_pan_set,"CARL.shapeit")
        } else if (pop == "VBI"){
            selected_panels <- c(current_pan_set,"VBI.shapeit")
        }else if (pop == "FVG"){
            selected_panels <- c(current_pan_set,"FVG.shapeit")
        }else if (pop == "INCIPE2"){
            selected_panels <- c(current_pan_set)
        }

        #subset the dataframe according to the selected panel set
        current_pop_all_panels_all_chr <- current_pop_all_panels_all_chr_complete[current_pop_all_panels_all_chr_complete$PANEL %in% current_pan_set,]
        
        current_pop_all_panels_all_chr <- current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$INFO >= 0),]
        current_pop_all_panels_all_chr$BIN3 <- 0
        maf_bins <- c(0,0.02,0.05,0.1,0.2,0.5)
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

        # remove all monomorphic data plus all the data with maf higher than the last bin
        current_pop_all_panels_all_chr_nomono <- current_pop_all_panels_all_chr[which(current_pop_all_panels_all_chr$BIN3 != 0),]

        ##before converting data for plotting and summarizing, I need to test differences in parameter distributions between panels
        # I can do a ks-test to asses differences between all panels
        # imputation_test(selected_panels)
        # we'll loop trough all the alternatives
        all_tests <- NULL
        for (alt in c("two.sided", "less", "greater")){
        print(alt)
        current_alt <- imputation_wilcox_test(current_pop_all_panels_all_chr_nomono,c("BIN3","PANEL"),"INFO",alt)        
        all_tests <- rbind(all_tests,current_alt)
        }


        ###################################################################################
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
        min_freq <- maf_bins[2]
        out_folder <- paste(base_folder,"/PLOTS/",current_date,"_INFO_",min_freq,"/",pan_set,"/",pop,"_",mode,sep="")
        dir.create(out_folder, recursive=T)

        # save(cdata_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_nomono.RData",sep=""))
        # save(cdata,file=paste(out_folder,"/",pop,"_",mode,"_cdata.RData",sep=""))
        write.table(cdata,file=paste(out_folder,"/",pop,"_",mode,"_cdata.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
        write.table(cdata_nomono,file=paste(out_folder,"/",pop,"_",mode,"_cdata_nomono.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
        write.table(all_tests,file=paste(out_folder,"/",pop,"_",mode,"_wilcox_nomono.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=T)
        # # load(paste(pop,"_cdata_nomono.RData",sep=""))

        #######################################################
        #plot without error bars of the info on MEDIAN MEAN and VARIANCE
        plot_range <- c("median","mean","variance")

        for (variable in plot_range){
            # variable <- plot_range[1]
            print(variable)
            
            pd <- position_dodge(0.001)
            pl <- ggplot(cdata_nomono)
            pl <- pl + aes(x = BIN3, y = cdata_nomono[,variable], colour=PANNELLO,group=PANNELLO)
            # pl <- pl + geom_errorbar(aes(ymax=mean+se, ymin=mean-se), width=0.2,colour="black")
            # pl <- pl + geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), width=0.2)
            pl <- pl + scale_x_continuous(breaks=cdata_nomono$BIN3,labels=cdata_nomono$BIN)
            pl <- pl + geom_line(position=pd,aes(group=PANNELLO),size=1.8) + geom_point(position=pd,size=3, shape=21, fill="white")
            pl <- pl + scale_color_discrete(name=paste(pop," on:"))
            pl <- pl + xlab("MAF bins") + ylab(paste(toupper(variable),"IMPUTE info score",sep=" "))
            pl <- pl + theme_bw()
            pl <- pl + theme(axis.text.x=element_text(size = rel(1.2), angle=45,hjust=1))
            pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
            pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
            pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
            # ggsave(filename=paste(getwd(),"/complete_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
            ggsave(filename=paste(out_folder,"/complete_13042016_",pop,"_",mode,"_info_",toupper(variable),".jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

        }

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
}