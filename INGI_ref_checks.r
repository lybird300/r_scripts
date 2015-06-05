#Script to plot stuff for quality check of INGI ref panels
rm(list=ls())
# source("/home/max/Work/script/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/assign_bins.r")
require(plotrix)
require(ggplot2)
library(plyr)

# define populations and reference panes
# pops <- c("CARL")
# pops <- c("FVG","INCIPE2","VBI")
pops <- c("FVG","INCIPE2","CARL")
# pops <- c("CARL","FVG","INCIPE2","VBI")
ingi_panels <- c("INGI.shapeit","CARL.shapeit","FVG.shapeit","VBI.shapeit")
gen_pop_ref_panels <- c("1000Gph1.shapeit","1000GP_Phase3.shapeit","INGI_1000GPh3.shapeit", "uk10k1kg.ref")
all_panels <- c(ingi_panels,gen_pop_ref_panels)
#for each population, we upload the interesting columns of the info file:
base_folder <- "/lustre/scratch113/projects/carl_seq/05272015_MERGED_REF_PANEL/IMPUTED"
#rsId, position,expected_af, info, type and concordances 
# maf_bins <- c(0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50)

for (pop in pops){
    # pop <- "CARL"
    print(pop)
    # current_pop_all_panels_all_chr <- NULL

    # if (pop == "CARL"){
    #     selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","CARL.shapeit")
    # } else if (pop == "VBI"){
    #     selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","VBI.shapeit")
    # }else if (pop == "FVG"){
    #     selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","FVG.shapeit")
    # }else if (pop == "INCIPE2"){
    #     selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref")
    # }

    # for(panel in selected_panels){
    #     # panel <- "1000Gph1.shapeit"
    #     print(panel)
    #     for (chr in 1:1){
    #         # chr <- 1
    #         pop_folder <- paste(base_folder,pop,panel,sep="/")
    #         current_pop_current_panel_current_chr_info_name <- paste(pop_folder,"/chr",chr,".gen_info_partial.gz",sep="")
    #         # current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,nrows=100000)
    #         current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,stringsAsFactors=F, comment.char="",colClasses=c("integer","character","integer",rep("numeric",2),"integer","NULL",rep("numeric",2),rep("character",2),rep("numeric",2)))
    #         current_pop_all_panels_all_chr <- rbind(current_pop_all_panels_all_chr,current_pop_current_panel_current_chr_info)

    #     }
    # }
    # # assign(paste("complete_",pop,"info",sep=""),current_pop_all_panels_all_chr)
    # cdata <- ddply(current_pop_all_panels_all_chr, c("BIN","PANEL"), summarise,
    #                N    = length(INFO),
    #                mean = mean(INFO),
    #                sd   = sd(INFO),
    #                se   = sd / sqrt(N)
    # )
    # cdata_nomono <- cdata[which(cdata$BIN != 0),]
    # cdata_nomono$PANNELLO <- cdata_nomono$PANEL
    # cdata_nomono$PANNELLO <- factor(cdata_nomono$PANNELLO,selected_panels)
    # save(cdata_nomono,file=paste(pop,"_cdata_nomono.RData",sep=""))
    load(paste(pop,"_cdata_nomono.RData",sep=""))
    #plot the same cohort stratifying by panel
    # The errorbars overlapped, so use position_dodge to move them horizontally

    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_nomono)
    pl <- pl + aes(x = BIN, y = mean, colour=PANNELLO,group=PANNELLO)
    pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005)
    pl <- pl + scale_x_continuous(breaks=cdata_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO)) + geom_point(position=pd,size=3, shape=21, fill="white")
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Average IMPUTE info score")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(0.8), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(0.8)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    ggsave(filename=paste(getwd(),"/complete_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
}

##############################################
# Plot R2 values by maf an cocnord
rm(list=ls())
# source("/home/max/Work/script/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/assign_bins.r")
require(plotrix)
require(ggplot2)
library(plyr)

# define populations and reference panes
pops <- c("CARL","FVG","INCIPE2","VBI")
# pops <- c("CARL")
ingi_panels <- c("INGI.shapeit","CARL.shapeit","FVG.shapeit","VBI.shapeit")
gen_pop_ref_panels <- c("1000Gph1.shapeit","1000GP_Phase3.shapeit","INGI_1000GPh3.shapeit", "uk10k1kg.ref")
all_panels <- c(ingi_panels,gen_pop_ref_panels)
#for each population, we upload the interesting columns of the info file:
base_folder <- "/lustre/scratch113/projects/carl_seq/05272015_MERGED_REF_PANEL/IMPUTED"

for (pop in pops){
    # pop <- "CARL"
    current_pop_all_panels_all_chr <- NULL

    if (pop == "CARL"){
        selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","CARL.shapeit")
    } else if (pop == "VBI"){
        selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","VBI.shapeit")
    }else if (pop == "FVG"){
        selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref","FVG.shapeit")
    }else if (pop == "INCIPE2"){
        selected_panels <- c("INGI.shapeit","1000Gph1.shapeit","INGI_1000GPh3.shapeit","uk10k1kg.ref")
    }

    for(panel in selected_panels){
        # panel <- "1000Gph1.shapeit"
        for (chr in 1:1){
            # chr <- 1
            pop_folder <- paste(base_folder,pop,panel,sep="/")
            current_pop_current_panel_current_chr_info_name <- paste(pop_folder,"/chr",chr,".gen_info_partial_t2.gz",sep="")
            # current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T,nrows=100000)
            current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T)
            current_pop_all_panels_all_chr <- rbind(current_pop_all_panels_all_chr,current_pop_current_panel_current_chr_info)

        }
    }
    # assign(paste("complete_",pop,"info",sep=""),current_pop_all_panels_all_chr)
    

    cdata_r2 <- ddply(current_pop_all_panels_all_chr, c("BIN","PANEL"), summarise,
                   N    = length(r2_TYPE0),
                   mean = mean(r2_TYPE0),
                   sd   = sd(r2_TYPE0),
                   se   = sd / sqrt(N)
    )
    cdata_r2_nomono <- cdata_r2[which(cdata_r2$BIN != 0),]
    cdata_r2_nomono <- cdata_r2_nomono[which(!is.na(cdata_r2_nomono$se)),]

    cdata_conc <- ddply(current_pop_all_panels_all_chr, c("BIN","PANEL"), summarise,
                   N    = length(CONCORD_TYPE0),
                   mean = mean(CONCORD_TYPE0),
                   sd   = sd(CONCORD_TYPE0),
                   se   = sd / sqrt(N)
    )
    cdata_conc_nomono <- cdata_conc[which(cdata_conc$BIN != 0),]
    cdata_conc_nomono <- cdata_conc_nomono[which(!is.na(cdata_conc_nomono$se)),]
    #plot the same cohort stratifying by panel
    # The errorbars overlapped, so use position_dodge to move them horizontally
    
    cdata_r2_nomono$PANNELLO <- cdata_r2_nomono$PANEL
    cdata_r2_nomono$PANNELLO <- factor(cdata_r2_nomono$PANNELLO,selected_panels)
    
    cdata_conc_nomono$PANNELLO <- cdata_conc_nomono$PANEL
    cdata_conc_nomono$PANNELLO <- factor(cdata_conc_nomono$PANNELLO,selected_panels)


    # pd <- position_dodge(0.001)
    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_r2_nomono)
    pl <- pl + aes(x = BIN, y = mean, colour=PANNELLO,group=PANNELLO)
    pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_conc_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO)) + geom_point(position=pd,size=3, shape=21, fill="white")
    # pl <- pl + scale_fill_manual("Ref. Panel", values=factor(PANNELLO),guide=FALSE)
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Average R2")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(0.8), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(0.8)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    ggsave(filename=paste(getwd(),"/complete_r2_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
    # ggsave(filename=paste(getwd(),"/complete_r2_",pop,"_infotest.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)


    pd <- position_dodge(0.001)
    pl <- ggplot(cdata_conc_nomono)
    pl <- pl + aes(x = BIN, y = mean, colour=PANNELLO,group=PANNELLO)
    pl <- pl + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.005,colour="black")
    pl <- pl + scale_x_continuous(breaks=cdata_conc_nomono$BIN)
    pl <- pl + geom_line(position=pd,aes(group=PANNELLO)) + geom_point(position=pd,size=3, shape=21, fill="white")
    pl <- pl + scale_color_discrete(name=paste(pop," on:"))
    pl <- pl + xlab("MAF bins") + ylab("Concordance")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(0.8), angle=45,hjust=1))
    pl <- pl + theme(axis.text.y=element_text(size = rel(0.8)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)),legend.position = c(0.8, 0.2))
    ggsave(filename=paste(getwd(),"/complete_conc_",pop,"_info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
}

