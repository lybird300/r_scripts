#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"

########################################
#Plot 9: some check stats on kinship data
i_pops <- c("FVG","VBI")

for (i_pop in i_pops){
  if(i_pop == "FVG"){
    new_kinship <- read.table("/lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/PLINK/KINSHIP/FVG_20140517.all.no_1st_deg.vcf.ibs0",header=T)
  }else if (i_pop == "VBI"){
    new_kinship <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/PLINK/KINSHIP/VBI_20140518.all.no_1st_deg.vcf.ibs0",header=T)
  }

  #check how many 1st degree samples we have still
  head(new_kinship[order(new_kinship$Kinship,decreasing=T),])

  # summaries
  summary(new_kinship)

  #boxplot for kinship
  jpeg(paste(base_folder,"PLOTS/9_kinship_",i_pop,".jpg",sep=""),width=800, height=800)
    par(lwd=2)
    boxplot(new_kinship$Kinship)
    # plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
    # lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col="red", verticals=TRUE, pch=46)
    # lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col="green", verticals=TRUE, pch=46)
    # lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col="blue", verticals=TRUE, pch=46)
    # legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
  dev.off()

  #plot also a heatmap

  
}

# merged_daf[which(merged_daf$VBI > 0 & (merged_daf$TSI > 0 | merged_daf$CEU > 0)),]
# merged_daf[which(merged_daf$FVG > 0 & (merged_daf$TSI > 0 | merged_daf$CEU > 0)),]

#########################################################################
# UTILITIES

# 1) Kinship heatmap

kinship_file <- "/lustre/scratch113/projects/carl_seq/20140410/20140410_VQSR2.5_reapply_138_excl/20140710_RELEASE/PLINK/KINSHIP/CARL_20140710.all.ibs0"
kinship_file <- "/lustre/scratch113/projects/carl_seq/20140410/20140410_VQSR2.5_reapply_138_excl/20140710_RELEASE/PLINK/KINSHIP/CARL_20140710.all.no1stdeg.ibs0"
current_gkin <- read.table(kinship_file,header=T)

current_gkin$FID1 <- as.character(current_gkin$FID1)
current_gkin$ID1 <- as.character(current_gkin$ID1)
current_gkin$FID2 <- as.character(current_gkin$FID2)
current_gkin$ID2 <- as.character(current_gkin$ID2)

#create a relatedness matrix for the heatmap foreach village
current_samples <- unique(c(current_gkin$ID1,current_gkin$ID2))

current_matrix <- matrix(0, nrow=length(current_samples), ncol=length(current_samples))

rownames(current_matrix) <- colnames(current_matrix) <- current_samples

#we need to assign the correct level to each factor!
# the problem here is that we have only the unique pairs, no duplicates, so this is different from the file generated from vcf tools
# erto_m[cbind(factor(erto_gkin$ID1,levels=levels(as.factor(c(erto_gkin$ID1, erto_gkin$ID2)))),factor(erto_gkin$ID2,levels=levels(as.factor(c(erto_gkin$ID1, erto_gkin$ID2)))))] <- erto_gkin$Kinship
current_matrix[cbind(factor(current_gkin$ID1,levels=sort(unique(c(current_gkin$ID1, current_gkin$ID2)))),factor(current_gkin$ID2,levels=sort(unique(c(current_gkin$ID1, current_gkin$ID2)))))] <- current_gkin$Kinship
current_matrix <- current_matrix + t(current_matrix)
#first check that the diagonal is all zeros, than assign 0.5 to reflect the kinship with himself
diag(current_matrix) <- 0.5

plot_path <- getwd() 
plot_name <- "CARL_geno_relatedness_heatmap_no1stdeg.jpg"
# jpeg(paste("/lustre/scratch113/teams/soranzo/users/mc14/INGI_VB/ARRAY/300700_overlap/SEQ_OVERLAP/KINSHIP/VBI_geno_relatedness_heatmap.jpg",sep="_"),width=10000, height=10000,)
library(pheatmap)

jpeg(paste(plot_path,plot_name,sep="/"),width=1000, height=1000)
  pheatmap(current_matrix)
dev.off()

