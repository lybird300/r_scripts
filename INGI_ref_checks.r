#Script to plot stuff for quality check of INGI ref panels
rm(list=ls())
# source("/home/max/Work/script/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/assign_bins.r")
require(plotrix)
require(ggplot2)

# define populations and reference panes
pops <- c("CARL","FVG","INCIPE2","VBI")
ingi_panels <- c("INGI.shapeit","CARL.shapeit","FVG.shapeit","VBI.shapeit")
gen_pop_ref_panels <- c("1000Gph1.shapeit","1000GP_Phase3.shapeit","INGI_1000GPh3.shapeit", "uk10k1kg.ref")
all_panels <- c(ingi_panels,gen_pop_ref_panels)
#for each population, we upload the interesting columns of the info file:
base_folder <- "/lustre/scratch113/projects/carl_seq/05272015_MERGED_REF_PANEL/IMPUTED"
#rsId, position,expected_af, info, type and concordances 
# maf_bins <- c(0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50)

for (pop in pops){
    # pop <- "CARL"
    current_pop_all_panels_all_chr <- NULL
    for(panel in all_panels){
        # panel <- "1000Gph1.shapeit"
        for (chr in 1:22){
            # chr <- 1
            pop_folder <- paste(base_folder,pop,panel,sep="/")
            current_pop_current_panel_current_chr_info_name <- paste(pop_folder,"/chr",chr,".gen_info_partial.gz",sep="")
            current_pop_current_panel_current_chr_info <-  read.table(current_pop_current_panel_current_chr_info_name,header=T)
            #add a panel and cohort column
            current_pop_current_panel_current_chr_info$panel <- panel
            current_pop_current_panel_current_chr_info$cohort <- pop
            current_pop_current_panel_current_chr_info$maf_bin <- "NA"
            # current_pop_current_panel_current_chr_info$maf_bin <- apply(current_pop_current_panel_current_chr_info,2,assign_bins(current_pop_current_panel_current_chr_info,"MAF",maf_bins))

            for (j in 1:(length(current_pop_current_panel_current_chr_info$MAF))){
                for (i in 1:length(maf_bins)){
                    if ( i == 1 ){
                        if ( current_pop_current_panel_current_chr_info$MAF[j] <= maf_bins[i] ){
                            current_pop_current_panel_current_chr_info$maf_bin[j] <- maf_bins[i]
                        }
                    }else{
                        if ( current_pop_current_panel_current_chr_info$MAF[j]  <= maf_bins[i] & current_pop_current_panel_current_chr_info$MAF[j]  > maf_bins[i-1] ){
                            current_pop_current_panel_current_chr_info$maf_bin[j]  <- maf_bins[i]
                        }
                    }
                }
            }

            current_pop_all_panels_all_chr <- rbind(current_pop_all_panels_all_chr,current_pop_current_panel_current_chr_info)

        }
    }
    assign(paste("complete_",pop,"info",sep=""),current_pop_all_panels_all_chr)
    
    #plot the same cohort stratifying by panel
    pl <- ggplot(current_pop_all_panels_all_chr)
    pl <- pl + aes(x = bin, y = INFO, colour=panel)
    # pl <- pl + scale_color_manual("Cohorts", values=pop_colors)
    pl <- pl + geom_line(size=1.5)
    # pl <- pl + geom_hline(aes(yintercept=0.95), linetype=2,colour="Lightgrey",size=1.2)
    pl <- pl + xlab("MAF bins") + ylab("IMPUTE info score")
    # pl <- pl + scale_x_continuous(limits=c(0,260))
    # pl <- pl + scale_y_continuous(limits=c(0,1))
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
    ggsave(filename=paste(getwd(),"/complete_",pop,"info.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
}

