######################################################################################################
# Plot 10: Ne from IBD
###############################################################################
# We calculated the IBD sharing segments using GERMLINE , than the sharing density for each chromosome.
# Now we'd like to: visualize this sharing density, select regions to remove because of too much sharing
# 
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")

pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
base_folder <- getwd()

ylab <- "Sharing density"
xlab <- "Chrom segment"

for (chr in seq(1,22)){
  current_shared_path <- paste(base_folder,"/CHR",chr,sep="")
  main <- paste("IBD sharing density across chr",chr)
  current_chr_all_pop <- NULL
  current_chr_all_pop_names <- NULL
  current_chr_all_pop_sd <- NULL

  for (pop in pops){
    current_chr_current_pop_path <- paste(current_shared_path,"/",pop,".",chr,".non_missing.match.shareDens",sep="")
    current_chr_current_pop <- read.table(current_chr_current_pop_path,header=F)
    current_chr_current_pop$cohort <- pop
    colnames(current_chr_current_pop) <- c("bin","share","cohort")
    print(pop)
    # current_chr_all_pop <- cbind(current_chr_all_pop,current_chr_current_pop[,2])
    current_chr_all_pop <- rbind(current_chr_all_pop,current_chr_current_pop)
    current_sd <- sd(current_chr_current_pop$share)
    current_mean <- mean(current_chr_current_pop$share)
    current_chr_current_pop_sd <- data.frame(share_sd=current_sd, mean_sd=(current_mean+current_sd),
      sd2=(current_mean+(2*current_sd)),
      sd5=(current_mean+(5*current_sd)),
      cohort=pop)
    current_chr_all_pop_sd <- rbind(current_chr_all_pop_sd, current_chr_current_pop_sd)
  }
  
  # current_chr_all_pop <- as.data.frame(cbind(current_chr_current_pop[,1],current_chr_all_pop))
  # colnames(current_chr_all_pop) <- c("bin",current_chr_all_pop_names)
  pl <- ggplot(current_chr_all_pop, aes(x = bin, y = share,col=cohort)) + 
    geom_point(size = 3) +
    # geom_line(size = 2) +
    # geom_boxplot(aes(fill = factor(cohort))) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    facet_wrap(~cohort) +
    + geom_hline(aes(yintercept = mean_sd), current_chr_all_pop_sd)
    + geom_hline(aes(yintercept = sd2), current_chr_all_pop_sd)
    + geom_hline(aes(yintercept = sd5), current_chr_all_pop_sd)
    # element_text()
    # + theme_bw() 
    # + theme(axis.title=element_text(face="bold.italic", 
    # size="12", color="brown"), legend.position="top")

  # jpeg(paste(current_shared_path,"/",chr,"_lines_dens.jpg",sep=""),width=1800, height=800)
  # jpeg(paste(current_shared_path,"/",chr,"_boxplot_dens.jpg",sep=""),width=1800, height=800)
  jpeg(paste(current_shared_path,"/",chr,"_point_dens.jpg",sep=""),width=1800, height=800)
    # par(oma=c(3,3,3,3),cex=1.4)
    pl
  dev.off()
}

# hist(current_chr_current_pop,
# col=all_cols$color,
# freq=FALSE,
# breaks=20,
# ylim=c(0, 180),
# legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
#we need to choose which analysis's windows size results we want to use
# win_size <- "w_20k"

    require(plotrix)
    all_cols <- col_pop(pops)

    #make some test plots....
    jpeg(paste(base_folder,"/10_ibd_cluster_all_size_",win_size,"_",den_size,"_plotrix.jpg",sep=""),width=1800, height=800)
      par(oma=c(3,3,3,3),cex=1.4)
      multhist(all_pop_cluster_size,
       col=all_cols$color,
       # freq=FALSE,
       ylab="Cluster Number",
       xlab="Sample number in cluster",
       # breaks=20,
       # ylim=c(0, 180),
       main="Sample number in cluster in all populations")
      legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
    dev.off()

    jpeg(paste(base_folder,"/10_ibd_cluster_all_length_",win_size,"_",den_size,"_plotrix.jpg",sep=""),width=1800, height=800)
      par(oma=c(3,3,3,3),cex=1.4)
      multhist(all_pop_cluster_length,
       col=all_cols$color,
       # freq=FALSE,
       ylab="Cluster Number",
       xlab="Cluster length (in Mb)",
       # breaks=12,
       # ylim=c(0, 80),
       main="Cluster length in all populations")
      legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
    dev.off()

    # cluster ibd shared 200K window LOD 5
    # q_CEU_w_200k: 95% -> 0.5123605 
    # q_TSI_w_200k: 95% -> 0.404625 
    # q_FVG_w_200k: 95% -> 1.020236 
    # q_VBI_w_200k: 95% -> 0.7859278 
    # q_CARL_w_200k: 95% -> 0.4817628 
    # q_Sauris_w_200k: 95% -> 0.7032557 
    # q_Erto_w_200k: 95% -> 0.7394982 
    # q_Illegio_w_200k: 95% -> 0.701432 
    # q_Resia_w_200k: 95% -> 1.086016 

    # cluster ibd shared 20K window LOD 5
    # q_CEU_w_20k: 95% -> 0.3514158 
    # q_TSI_w_20k: 95% -> 0.3765026 
    # q_FVG_w_20k: 95% -> 0.8347775 
    # q_VBI_w_20k: 95% -> 0.580388 
    # q_CARL_w_20k: 95% -> 0.3499446 
    # q_Sauris_w_20k: 95% -> 0.6053076 
    # q_Erto_w_20k: 95% -> 0.5024974 
    # q_Illegio_w_20k: 95% -> 0.4840495 
    # q_Resia_w_20k: 95% -> 0.84094 

    pop_colors <- col_pop(pops)
    base_folder <- getwd()

    # jpeg(paste(base_folder,"/6_roh_all_5POP",chr,".jpg",sep=""),width=1000, height=1000)
    # jpeg(paste(base_folder,"/10_ibd_cluster_all_length_",win_size,"_ECDF.jpg",sep=""),width=1000, height=1000)
    jpeg(paste(base_folder,"/10_ibd_cluster_all_length_",win_size,"_",den_size,"_lines.jpg",sep=""),width=1000, height=1000)
      par(lwd=4,cex=1.5)
      plot(M_CEU,CEU_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "CEU"),1],main="",xlab="Cluster size (Mb)",ylab="% of Samples sharing cluster with the same length", xlim=c(0,max(xmax)), yaxs='i',verticals=TRUE, pch=46)
      # plot(CEU_tot_cluster$LENGTH/1000000,CEU_tot_cluster$SIZE,col=pop_colors[which(pop_colors$pop == "CEU"),1],main="",xlab="Sample number per cluster",ylab="Cluster size in Mb")
      lines(M_TSI,TSI_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "TSI"),1],yaxs='i',verticals=TRUE, pch=46)
      lines(M_VBI,VBI_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "VBI"),1],yaxs='i',verticals=TRUE,pch=46)
      lines(M_FVG,FVG_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "FVG"),1],yaxs='i',verticals=TRUE,pch=46)
      lines(M_CARL,CARL_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "CARL"),1],yaxs='i',verticals=TRUE,pch=46)
      lines(M_Erto,Erto_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "Erto"),1],yaxs='i',verticals=TRUE,pch=46)
      lines(M_Illegio,Illegio_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "Illegio"),1],yaxs='i',verticals=TRUE,pch=46)
      lines(M_Resia,Resia_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "Resia"),1],yaxs='i',verticals=TRUE,pch=46)
      lines(M_Sauris,Sauris_tot_cluster$LENGTH/1000000,col=pop_colors[which(pop_colors$pop == "Sauris"),1],yaxs='i',verticals=TRUE,pch=46)
      abline(h=0.95,col='grey',lty='dashed')
      #define parameters for legend
      leg_txt <- c(pop_colors[which(pop_colors$pop == "CEU"),2],pop_colors[which(pop_colors$pop == "TSI"),2],pop_colors[which(pop_colors$pop == "VBI"),2],pop_colors[which(pop_colors$pop == "FVG"),2],pop_colors[which(pop_colors$pop == "CARL"),2],pop_colors[which(pop_colors$pop == "Erto"),2],pop_colors[which(pop_colors$pop == "Illegio"),2],pop_colors[which(pop_colors$pop == "Resia"),2],pop_colors[which(pop_colors$pop == "Sauris"),2])
      bkg <- c(pop_colors[which(pop_colors$pop == "CEU"),1],pop_colors[which(pop_colors$pop == "TSI"),1],pop_colors[which(pop_colors$pop == "VBI"),1],pop_colors[which(pop_colors$pop == "FVG"),1],pop_colors[which(pop_colors$pop == "CARL"),1],pop_colors[which(pop_colors$pop == "Erto"),1],pop_colors[which(pop_colors$pop == "Illegio"),1],pop_colors[which(pop_colors$pop == "Resia"),1],pop_colors[which(pop_colors$pop == "Sauris"),1])
      legend("bottomright",pch =c(rep(22,length(pops))),legend=leg_txt, pt.lwd=2,pt.cex=2,pt.bg=bkg,col=c(rep('black',length(pops))),ncol=4,bty="n")
    dev.off()
  
}
}

