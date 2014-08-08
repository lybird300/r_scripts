#r script to plot for enza
rm(list=ls())
# source("/home/max/Work/script/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"

########################################################################
#PLOT 1A: same as plot 1 but based on MAF: for the DAF categories, we plot the maf spectrum
#upload data for each population:
# chr <- "22"
#the order of all pops is due to the merge order
pops <- c("CEU","TSI","VBI","FVG","CARL")
#wrapped in this for-loop, I'll write a table for each chr, so we can easily import the thing to plot, after...
all_pop_MAF <- NULL
for (pop in pops) {
  print(pop)
  all_chr_MAF <- NULL
  for (chr in 1:22) {
  # for (chr in 20:22) {
  print(chr)
  base_folder <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/CHR",chr,sep="")

    # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
    pop_table_file <- paste(base_folder,"/",pop,".chr",chr,".tab.gz",sep="")
    # for files with not mac=1
    # pop_table_file <- paste(pop,".chr",chr,".not_fixed.not_MAC1.tab.gz",sep="")

    pop_table <- read.table(pop_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
    # pop_table <- read.table(pop_table_file,header=F,skip=1,nrows=100000,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
    colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
    pop_table$MAF <- as.numeric(as.character(pop_table$MAF))

    #remove nas
    #we want to keep the NA's because those are the novel variants!
    # pop_table <- pop_table[which(!is.na(pop_table$DAF)),]
    
    #remove fixed variants (the ones with DAF == 0 or == 1) ---> those are also the ones with MAF 0, basically!
    pop_table <- pop_table[which(pop_table$MAF != 0),]
    pop_table <- pop_table[which(pop_table$MAF != 1),]

    #remove MAC = 1
    # pop_table <- pop_table[which(pop_table$MAC > 1),]

    # current_hist <- hist(pop_table$MAF,plot=F,breaks=20)

    # #use relative frequency sites per bin count/total variants number
    # current_pop_MAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pop,length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS)),chr_length=(length(pop_table$POS))))
    
    all_chr_MAF <- rbind(all_chr_MAF,pop_table)

  }
  all_pop_MAF <- append(all_pop_MAF,list(all_chr_MAF$MAF))
}

pop_col <- col_pop(pops)

all_pop_hist <- multhist(all_pop_MAF,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_MAF_table <- as.data.frame(cbind((all_pop_hist[[1]]$mids),t(all_pop_hist[[2]])))
colnames(all_pop_MAF_table) <- c("breaks",pops)

write.table(all_pop_MAF_table,file=paste("all_pop_MAF_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
  
require(plotrix)
base_folder <- getwd()
jpeg(paste(base_folder,"/1_all_pop_MAF_plotrix.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  multhist(all_pop_MAF,
   col=pop_col$color,
   freq=FALSE,
   ylab="Relative Frequency (N sites/Tot sites in freq bin)(%)",
   xlab="MAF",
   breaks=20,
   ylim=c(0, 40),
   main="MAF in all populations")
  legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
dev.off()

################################################################################################################
#WE NEED TO UPLOAD DATA ALSO FOR NOVEL SITES from here:
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
base_novel_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/NOVEL"
pops_ingi_novel <- c("VBI_n","FVG_n","CARL_n")
pops_ingi <- c("VBI","FVG","CAR")

#we work by chr and by population
#upload data for each population about novel sites
all_pop_novel_MAF <- NULL
for (pop in pops_ingi) {
  print(pop)
  current_pop_novel <- NULL
  for (chr in 1:22) {
    print(chr)
    #read the current file for this population and this chromosome
    current_chr_novel <- read.table(paste(base_novel_folder,"/",pop,".",chr,".n.maf",sep=""), sep="\t",header=F)
    colnames(current_chr_novel) <- c("CHROM","POZ","POS","TYPE","CEU","TSI","VBI","FVG","CAR")
    current_pop_col <- grep(pop,colnames(current_chr_novel))
    current_chr_novel <- current_chr_novel[which(current_chr_novel[,current_pop_col] != 0),]
    
    current_pop_novel <-rbind(current_pop_novel,current_chr_novel[,c(1,2,3,4,current_pop_col)])
  }
  
  current_pop_col <- grep(pop,colnames(current_pop_novel))
  all_pop_novel_MAF <- append(all_pop_novel_MAF,list(current_pop_novel[,current_pop_col]))

}
pop_col <- col_pop(pops_ingi_novel)
base_folder <- getwd()

novel_hist <- multhist(all_pop_novel_MAF,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_novel_MAF_table <- as.data.frame(cbind((novel_hist[[1]]$mids),t(novel_hist[[2]])))
colnames(all_pop_novel_MAF_table) <- c("breaks",pops_ingi_novel)

write.table(all_pop_novel_MAF_table,file=paste("all_ingi_pop_MAF_novel_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

require(plotrix)
jpeg(paste(base_folder,"/1_all_INGI_novel_MAF_plotrix.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  multhist(all_pop_novel_MAF,
   col=pop_col$color,
   freq=FALSE,
   ylab="Relative Frequency (N sites/Tot sites in freq bin)(%)",
   xlab="MAF",
   breaks=20,
   ylim=c(0, 40),
   main="MAF in all populations")
  legend("topright",pch =c(rep(22,length(pop_col[,1]))),pt.lwd=2,pt.cex=2,pt.bg=pop_col[,1],col=c(rep('black',length(pops_ingi_novel))),legend=pop_col[,2], ncol=2,bty="n")
dev.off()

################################################################################################################
#now upload the private and shared plots and add them to the previous plot
# merged_daf <- read.table(paste0(base_folder,"INPUT_FILES/INGI_chr",chr,".merged_daf.tab.gz"),skip=1,header=F)
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
#we need to use the files with MAF info, but created after filtering by DAF: we need to put all chr together
pops_ingi_class <- c("VBI_p","FVG_p","CARL_p","VBI_s","FVG_s","CARL_s")
pops_ingi <- c("VBI","FVG","CARL")
var_class <- c("private","shared")

for (type in var_class) {
  print(type)
  all_pop_current_type <- NULL

  for (pop in pops_ingi) {
    print(pop)
    current_pop_current_type <- NULL
    for (chr in 1:22) {
      print(chr)
      #read the current file for this population and this chromosome
      base_type_folder <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/CHR",chr,sep="")
      current_filename <- paste(base_type_folder,"/INGI_chr",chr,".merged_maf.tab.gz.",pop,".",type,".tab.gz",sep="")
      current_chr_current_type <- read.table(current_filename, sep="\t",header=F)
      colnames(current_chr_current_type) <- c("CHROM","POZ","POS","TYPE","CEU","TSI","VBI","FVG","CARL")
      current_pop_col <- grep(pop,colnames(current_chr_current_type))
      current_chr_current_type <- current_chr_current_type[which(current_chr_current_type[,current_pop_col] != 0),]
      
      current_pop_current_type <-rbind(current_pop_current_type,current_chr_current_type[,c(1,2,3,4,current_pop_col)])
    }
    
    current_pop_col <- grep(pop,colnames(current_pop_current_type))
    all_pop_current_type <- append(all_pop_current_type,list(current_pop_current_type[,current_pop_col]))

  }
  all_pop_current_type_table_name <- paste("all_pop_MAF_",type,sep="")
  assign(all_pop_current_type_table_name,all_pop_current_type)
  
}

#we now have lists for private and shared sites
all_pop_MAF_private_shared <- append(all_pop_MAF_private,all_pop_MAF_shared)

#now we should be able to proceed as before...and plot
pop_col <- col_pop(pops_ingi_class)
base_folder <- getwd()

all_pop_MAF_private_shared_hist <- multhist(all_pop_MAF_private_shared,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_MAF_private_shared_table <- as.data.frame(cbind((all_pop_MAF_private_shared_hist[[1]]$mids),t(all_pop_MAF_private_shared_hist[[2]])))
colnames(all_pop_MAF_private_shared_table) <- c("breaks",pops_ingi_class)

write.table(all_pop_MAF_private_shared_table,file=paste("all_pop_MAF_private_shared_table.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

require(plotrix)
jpeg(paste(base_folder,"/1_all_INGI_private_shared_MAF_plotrix.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  multhist(all_pop_MAF_private_shared,
   col=pop_col$color,
   density=pop_col$density,
   freq=FALSE,
   ylab="Relative Frequency (N sites/Tot sites in freq bin)(%)",
   xlab="MAF",
   breaks=20,
   ylim=c(0, 50),
   main="MAF in all populations")
  legend("topright",pch =c(rep(22,length(pop_col[,1]))),pt.lwd=2,pt.cex=2,pt.bg=pop_col[,1],col=c(rep('black',length(pops_ingi_novel))),legend=pop_col[,2], ncol=2,bty="n")
dev.off()

#now we need to put all together to plot the complete spectum for all classes!

all_pop_all_MAF <- append(append(all_pop_MAF,all_pop_novel_MAF),all_pop_MAF_private_shared)

#we need pops in the same order as the list
all_pops <- c(pops,pops_ingi_novel,pops_ingi_class)