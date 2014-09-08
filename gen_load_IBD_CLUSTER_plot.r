######################################################################################################
# Plot 10: IBD clustering
###############################################################################
#Same as roh, but we consider samples's pairs:
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
# input_format <- "PLINK"
# chr <- "22"
pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
base_folder <- getwd()

#we need to choose which analysis's windows size results we want to use
win_sizes <- c("w_20k","w_100k","w_200k","w_400k")
dens_sizes <- c("0.4","0.5","0.6","0.7","0.8")

# win_size <- "w_20k"
for (win_size in win_sizes){
  for (den_size in dens_sizes){
    #inizialize variables
    xmax <- NULL
    all_pop_cluster <- NULL
    all_pop_cluster_size <- NULL
    all_pop_cluster_length <- NULL
    all_pop_cluster_names <- NULL
    for (pop in pops) {
      current_pop_cluster <- NULL
      for (chr in 1:22){
        current_pop_current_chr_ibd_cluster_file <- paste(base_folder,"/",win_size,"/",den_size,"/",pop,"/",pop,".chr",chr,".ibd.out.clst.uniq.size",sep="")
        current_pop_current_chr_cluster <- try(read.table(current_pop_current_chr_ibd_cluster_file,header=F))
        if(!inherits(current_pop_current_chr_cluster,'try-error')){
          colnames(current_pop_current_chr_cluster) <- c("CLST","CHROM","START","END","SIZE")
          current_pop_current_chr_cluster$LENGTH <- current_pop_current_chr_cluster$END - current_pop_current_chr_cluster$START
          current_pop_cluster <- rbind(current_pop_cluster,current_pop_current_chr_cluster)
        }
      }
      # all_pop_cluster <- cbind(all_pop_cluster,current_pop_cluster$SIZE,current_pop_cluster$LENGTH)
      all_pop_cluster_size <- append(all_pop_cluster_size,list(current_pop_cluster$SIZE))
      all_pop_cluster_length <- append(all_pop_cluster_length,list(current_pop_cluster$LENGTH/1000000))

      assign(paste(pop,"tot_cluster",sep="_"),current_pop_cluster)  
      assign(paste("M",pop,sep="_"),ecdf(current_pop_cluster$LENGTH/1000000))
      xmax <- c(xmax,summary(ecdf(current_pop_cluster$LENGTH/1000000))[6])
      # all_pop_cluster_names <- c(all_pop_cluster_names,paste(pop,"_SIZE",sep=""),paste(pop,"_LENGTH",sep=""))
    }

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

    #plot the same way used for roh and ibd with ecdf
    assign(paste("q_CEU",win_size,den_size,sep="_"),quantile(CEU_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_TSI",win_size,den_size,sep="_"),quantile(TSI_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_FVG",win_size,den_size,sep="_"),quantile(FVG_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_VBI",win_size,den_size,sep="_"),quantile(VBI_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_CARL",win_size,den_size,sep="_"),quantile(CARL_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_Sauris",win_size,den_size,sep="_"),quantile(Sauris_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_Erto",win_size,den_size,sep="_"),quantile(Erto_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_Illegio",win_size,den_size,sep="_"),quantile(Illegio_tot_cluster$LENGTH/1000000,c(.95)))
    assign(paste("q_Resia",win_size,den_size,sep="_"),quantile(Resia_tot_cluster$LENGTH/1000000,c(.95)))

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

