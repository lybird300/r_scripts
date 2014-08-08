######################################################################################################
# Plot 10: IBD clustering
###############################################################################
#Same as roh, but we consider samples's pairs:
# input_format <- "PLINK"
chr <- "22"
pops <- c("CEU","TSI","VBI","FVG")

input_format <- "BEAGLE"
xmax <- NULL
for (pop in pops) {
  pop_ibd_cluster_file <- paste(pop,"clst.uniq.size",sep=".")
  current_cluster <- read.table(pop_ibd_cluster_file,header=F)
  colnames(current_cluster) <- c("CLST","START","END","SIZE")
  current_cluster$LENGTH <- current_cluster$END - current_cluster$START
  assign(paste(pop,"ibd_cluster",sep="_"),current_cluster)

}

ymax <- c(max((density(FVG_ibd_cluster$SIZE))$y),max((density(CEU_ibd_cluster$SIZE))$y),max((density(VBI_ibd_cluster$SIZE))$y),max((density(TSI_ibd_cluster$SIZE))$y))

jpeg(paste(base_folder,"PLOTS/10_ibd_cluster.jpg",sep=""),width=800, height=800)
  par(lwd=3,cex=1.2)
  plot(density(CEU_ibd_cluster$SIZE),main="",xlab="IBD cluster size", ylim=c(0,max(ymax)), pch=46)
  lines(density(FVG_ibd_cluster$SIZE),col="red", pch=46)
  lines(density(TSI_ibd_cluster$SIZE),col="green", pch=46)
  lines(density(VBI_ibd_cluster$SIZE),col="blue", pch=46)
  legend("topright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()

jpeg(paste(base_folder,"PLOTS/10_ibd_cluster_CEU.jpg",sep=""),width=800, height=800)
hist(CEU_ibd_cluster$SIZE,col="black",main="CEU")
dev.off()

jpeg(paste(base_folder,"PLOTS/10_ibd_cluster_FVG.jpg",sep=""),width=800, height=800)
hist(FVG_ibd_cluster$SIZE,col="red",main="FVG")
dev.off()

jpeg(paste(base_folder,"PLOTS/10_ibd_cluster_TSI.jpg",sep=""),width=800, height=800)
hist(TSI_ibd_cluster$SIZE,col="green",main="TSI")
dev.off()

jpeg(paste(base_folder,"PLOTS/10_ibd_cluster_VBI.jpg",sep=""),width=800, height=800)
hist(VBI_ibd_cluster$SIZE,col="blue",main="VBI")
dev.off()

