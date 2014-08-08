#define function for population colors in plots
col_pop <- function(pops){
  all_cols <- NULL

  for(i in 1:length(pops)){
  if(pops[i] == "CEU"){
    cur_col <- "#E8D0A9"
  }
  if(pops[i] == "TSI"){
    cur_col <- "#B7AFA3"
  }
  if(pops[i] == "FVG"){
    cur_col <- "#5A5E61"
  }
  if(pops[i] == "FVG_p"){
    cur_col <- "#5A5E61"
  }
  if(pops[i] == "FVG_s"){
    cur_col <- "#5A5E61"
  }
  if(pops[i] == "Erto"){
    cur_col <- "#96D4CA"
  }
  if(pops[i] == "Illegio"){
    cur_col <- "#95C0BB"
  }
  if(pops[i] == "Resia"){
    cur_col <- "#85C0E7"
  }
  if(pops[i] == "Sauris"){
    cur_col <- "#66939E"
  }
  if(pops[i] == "VBI"){
    cur_col <- "#DF5E5E"
  }
  if(pops[i] == "VBI_p"){
    cur_col <- "#DF5E5E"
  }
  if(pops[i] == "VBI_s"){
    cur_col <- "#DF5E5E"
  }
   if(pops[i] == "CARL"){
    cur_col <- "#F7A6A6"
  }
  if(pops[i] == "CARL_p"){
    cur_col <- "#F7A6A6"
  }
  if(pops[i] == "CARL_s"){
    cur_col <- "#F7A6A6"
  }
  pop_col <- cbind(cur_col,pops[i])
  all_cols <- rbind(all_cols,pop_col)
}
all_cols <- as.data.frame(all_cols)
colnames(all_cols) <- c("color","pop")
all_cols$color <- as.character(all_cols$color)
all_cols$pop <- as.character(all_cols$pop)


return (all_cols)
}

#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"
########################################################################
#plot 1 DAF spectrun for each population
#upload data for each population:

chr <- "22"
pops <- c("CEU","TSI","VBI","FVG","CARL")

all_pop_DAF <- NULL

for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  # pop_table_file <- paste(pop,".chr",chr,".tab.gz",sep="")
  # for files with not mac=1
  pop_table_file <- paste(pop,".chr",chr,".not_fixed.not_MAC1.tab.gz",sep="")

  pop_table <- read.table(pop_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="")
  colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  pop_table$DAF <- as.numeric(as.character(pop_table$DAF))

  #remove nas
  pop_table <- pop_table[which(!is.na(pop_table$DAF)),]
  #remove fixed variants (the ones with DAF == 0 or == 1)
  pop_table <- pop_table[which(pop_table$DAF != 0),]
  pop_table <- pop_table[which(pop_table$DAF != 1),]

  #remove MAC = 1
  pop_table <- pop_table[which(pop_table$MAC > 1),]


  current_hist <- hist(pop_table$DAF,plot=F)

  # current_pop_DAF <- cbind(as.data.frame(current_hist$breaks),as.data.frame(current_hist$counts),as.data.frame(rep(pop,length(current_hist$counts))))
  # current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=(current_hist$counts/length(pop_table$POS)),pop=rep(pop,length(current_hist$counts))))
  #use relative frequency sites per bin count/total variants number
  current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pop,length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS))))
  
  all_pop_DAF <- rbind(all_pop_DAF,current_pop_DAF)

}

all_pop_DAF$breaks <- as.numeric(as.character(all_pop_DAF$breaks))
all_pop_DAF$counts <- as.numeric(as.character(all_pop_DAF$counts))
all_pop_DAF$rel_count <- as.numeric(as.character(all_pop_DAF$rel_count))
all_pop_DAF$pop <- as.character(all_pop_DAF$pop)

all_pop_DAF_table <- all_pop_DAF[which(all_pop_DAF$pop == "CEU"),1]
all_pop_DAF_table_names <- c("breaks")

for(pop_i in pops){
  actual_pop <- all_pop_DAF[which(all_pop_DAF$pop == pop_i),]
  all_pop_DAF_table_names <- c(all_pop_DAF_table_names,paste(pop_i,"count",sep="_"),paste(pop_i,"rel_count",sep="_"))

  all_pop_DAF_table <- cbind(all_pop_DAF_table,cbind(actual_pop$counts,actual_pop$rel_count*100))
}
all_pop_DAF_table <- as.data.frame(all_pop_DAF_table)
colnames(all_pop_DAF_table) <- all_pop_DAF_table_names

# write.table(all_pop_DAF_table,file=paste("All_pop_DAF_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_DAF_table,file=paste("All_pop_DAF_count_no_MAC1_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#absolute count
barplot3 <- as.matrix(t(all_pop_DAF_table[,c(2,4,6,8)]))

#relative site count
barplot4 <- as.matrix(t(all_pop_DAF_table[,c(3,5,7,9)]))


all_cols <- col_pop(pops)

# jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF.jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_no_MAC1.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot3,
    beside=T,
    main="DAF in all populations",
    names.arg=all_pop_DAF_table$breaks,
    xlab="",
    ylab="",col=all_cols[,1])
    legend("topright",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=length(all_cols[,1]))
    mtext(1, text = "Annotations", line = 4,cex=1.4)
    mtext(2, text = "Site count", line = 4,cex=1.4)
dev.off()

# jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_rf.jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_rf_no_MAC1.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot4,
    beside=T,
    main="DAF in all populations",
    names.arg=all_pop_DAF_table$breaks,
    xlab="",
    ylab="",col=all_cols[,1])
    legend("topright",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=length(all_cols[,1]))
    mtext(1, text = "Annotations", line = 4,cex=1.4)
    mtext(2, text = "Relative Frequency (N sites/Tot sites with DAF information )(%)", line = 4,cex=1.4)
dev.off()


#now upload the private and shared plots and add them to the previous plot
# merged_daf <- read.table(paste0(base_folder,"INPUT_FILES/INGI_chr",chr,".merged_daf.tab.gz"),skip=1,header=F)
# colnames(merged_daf) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","FVG")

# merged_daf$CEU <- as.numeric(as.character(merged_daf$CEU))
# merged_daf$TSI <- as.numeric(as.character(merged_daf$TSI))
# merged_daf$VBI <- as.numeric(as.character(merged_daf$VBI))
# merged_daf$FVG <- as.numeric(as.character(merged_daf$FVG))

pop_VBI_private <- read.table(paste(base_folder,"INPUT_FILES/VBI_private_chr",chr,".merged_daf.tab.gz",sep=""),header=F)
pop_FVG_private <- read.table(paste(base_folder,"INPUT_FILES/FVG_private_chr",chr,".merged_daf.tab.gz",sep=""),header=F)
pop_VBI_shared <- read.table(paste(base_folder,"INPUT_FILES/VBI_shared_chr",chr,".merged_daf.tab.gz",sep=""),header = F)
pop_FVG_shared <- read.table(paste(base_folder,"INPUT_FILES/FVG_shared_chr",chr,".merged_daf.tab.gz",sep=""),header = F)

# for mac=1 exclusion
pop_VBI_private <- read.table(paste(base_folder,"INPUT_FILES/VBI_private_chr",chr,".merged_daf.not_MAC1.tab.gz",sep=""),header=F)
pop_FVG_private <- read.table(paste(base_folder,"INPUT_FILES/FVG_private_chr",chr,".merged_daf.not_MAC1.tab.gz",sep=""),header=F)
pop_VBI_shared <- read.table(paste(base_folder,"INPUT_FILES/VBI_shared_chr",chr,".merged_daf.not_MAC1.tab.gz",sep=""),header = F)
pop_FVG_shared <- read.table(paste(base_folder,"INPUT_FILES/FVG_shared_chr",chr,".merged_daf.not_MAC1.tab.gz",sep=""),header = F)

colnames(pop_VBI_private) <- c("CHR","POZ","POS","VT","CEU","TSI","DAF","FVG")
pop_VBI_private$CEU <- as.numeric(as.character(pop_VBI_private$CEU))
pop_VBI_private$TSI <- as.numeric(as.character(pop_VBI_private$TSI))
pop_VBI_private$DAF <- as.numeric(as.character(pop_VBI_private$DAF))
pop_VBI_private$FVG <- as.numeric(as.character(pop_VBI_private$FVG))

colnames(pop_FVG_private) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","DAF")
pop_FVG_private$CEU <- as.numeric(as.character(pop_FVG_private$CEU))
pop_FVG_private$TSI <- as.numeric(as.character(pop_FVG_private$TSI))
pop_FVG_private$VBI <- as.numeric(as.character(pop_FVG_private$VBI))
pop_FVG_private$DAF <- as.numeric(as.character(pop_FVG_private$DAF))

colnames(pop_VBI_shared) <- c("CHR","POZ","POS","VT","CEU","TSI","DAF","FVG")
pop_VBI_shared$CEU <- as.numeric(as.character(pop_VBI_shared$CEU))
pop_VBI_shared$TSI <- as.numeric(as.character(pop_VBI_shared$TSI))
pop_VBI_shared$DAF <- as.numeric(as.character(pop_VBI_shared$DAF))
pop_VBI_shared$FVG <- as.numeric(as.character(pop_VBI_shared$FVG))

colnames(pop_FVG_shared) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","DAF")
pop_FVG_shared$CEU <- as.numeric(as.character(pop_FVG_shared$CEU))
pop_FVG_shared$TSI <- as.numeric(as.character(pop_FVG_shared$TSI))
pop_FVG_shared$VBI <- as.numeric(as.character(pop_FVG_shared$VBI))
pop_FVG_shared$DAF <- as.numeric(as.character(pop_FVG_shared$DAF))

summary(pop_VBI_private)
summary(pop_FVG_private)
summary(pop_VBI_shared)
summary(pop_FVG_shared)

dim(pop_VBI_private)
dim(pop_FVG_private)
dim(pop_VBI_shared)
dim(pop_FVG_shared)

# after removal of MAC=1
# VBI_private : 21170
# FVG_private : 22631
# VBI_shared : 109341
# FVG_shared : 104707


#this is useless now...
# pop_VBI_private <- pop_VBI_private[which(pop_VBI_private$DAF != 1),]
# pop_FVG_private <- pop_FVG_private[which(pop_FVG_private$DAF != 1),]
# pop_VBI_shared <- pop_VBI_shared[which(pop_VBI_shared$DAF != 1),]
# pop_FVG_shared <- pop_FVG_shared[which(pop_FVG_shared$DAF != 1),]

pops_ingi <- c("VBI_p","FVG_p","VBI_s","FVG_s")
pops_ingi_files <- c("pop_VBI_private","pop_FVG_private","pop_VBI_shared","pop_FVG_shared")

for (k in 1:length(pops_ingi_files)){
  current_pop <- get(pops_ingi_files[k])

  current_hist <- hist(current_pop$DAF,plot=F)

  #added relative frequency sites per bin count/total variants number for the current category! shared or private!!
  current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pops_ingi[k],length(current_hist$counts)),rel_count=(current_hist$counts/length(current_pop$POS))))
  current_pop_DAF$breaks <- as.numeric(as.character(current_pop_DAF$breaks))
  current_pop_DAF$counts <- as.numeric(as.character(current_pop_DAF$counts))
  current_pop_DAF$rel_count <- as.numeric(as.character(current_pop_DAF$rel_count))
  current_pop_DAF$pop <- as.character(current_pop_DAF$pop)

  #add again to all_pop_DAF to have a complete plot
  all_pop_DAF <- rbind(all_pop_DAF,current_pop_DAF)
  
}

all_pop_DAF_table_sp <- all_pop_DAF[which(all_pop_DAF$pop == "CEU"),1]
all_pop_DAF_table_sp_names <- c("breaks")

pops2 <- sort(c(pops,pops_ingi))

for(pop_i2 in pops2){
  actual_pop2 <- all_pop_DAF[which(all_pop_DAF$pop == pop_i2),]
  all_pop_DAF_table_sp_names <- c(all_pop_DAF_table_sp_names,paste(pop_i2,"count",sep="_"),paste(pop_i2,"rel_count",sep="_"))

  all_pop_DAF_table_sp <- cbind(all_pop_DAF_table_sp,cbind(actual_pop2$counts,actual_pop2$rel_count*100))
}

all_pop_DAF_table_sp <- as.data.frame(all_pop_DAF_table_sp)
colnames(all_pop_DAF_table_sp) <- all_pop_DAF_table_sp_names

write.table(all_pop_DAF_table_sp,file=paste("All_pop_DAF_sp_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
# without MAC=1
write.table(all_pop_DAF_table_sp,file=paste("All_pop_DAF_sp_count_no_MAC1_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#absolute count
barplot5 <- as.matrix(t(all_pop_DAF_table_sp[,c(2,4,6,8,10,12,14,16)]))

#relative site count
barplot6 <- as.matrix(t(all_pop_DAF_table_sp[,c(3,5,7,9,11,13,15,17)]))

#relative count but only for private sites
barplot7 <- as.matrix(t(all_pop_DAF_table_sp[,c(3,5,9,11,13,15)]))

all_cols <- col_pop(pops)
# jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_sp.jpg",sep=""),width=1800, height=800)
# no MAC=1
jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_sp_no_MAC1.jpg",sep=""),width=1800, height=800)
    par(oma=c(3,3,3,3),cex=1.4)
    barplot(barplot5,
      beside=T,
      main="DAF in all populations",
      names.arg=all_pop_DAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[,1]
      )
      mtext(1, text = "DAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Site count", line = 4,cex=1.4)
      legend(x="top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=8)
  dev.off()

  # jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_sp_rf.jpg",sep=""),width=1800, height=800)
  jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_sp_rf_no_MAC1.jpg",sep=""),width=1800, height=800)
      par(oma=c(3,3,3,3),cex=1.4)
    barplot(barplot6,
      beside=T,
      main="DAF in all populations",
      names.arg=all_pop_DAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[,1],ylim=c(0,60))
      mtext(1, text = "DAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Relative Frequency (N sites/Tot sites per category)(%)", line = 4,cex=1.4)
      legend("top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=8)
  dev.off()

  #plot only private, not shared
  # jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_p_rf.jpg",sep=""),width=1800, height=800)
  jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_p_rf_no_MAC1.jpg",sep=""),width=1800, height=800)
      par(oma=c(3,3,3,3),cex=1.4)
    barplot(barplot7,
      beside=T,
      main="DAF in all populations",
      names.arg=all_pop_DAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[c(1,2,3,5,6,7),1],ylim=c(0,60))
      mtext(1, text = "DAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Relative Frequency (N sites/Tot sites per category)(%)", line = 4,cex=1.4)
      legend("top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[c(1,2,3,5,6,7),1],legend=all_cols[c(1,2,3,5,6,7),2], ncol=6)
  dev.off()


########################################################################
#PLOT 1A: same as plot 1 but based on MAF: for the DAF categories, we plot the maf spectrum
#upload data for each population:
chr <- "22"

pops <- c("CEU","TSI","VBI","FVG","CARL")

#wrapped in this for-loop, I'll write a table for each chr, so we can easily import the thing to plot, after...
for (chr in 1:22) {
  base_folder <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/CHR",chr,sep="")
  print(chr)
  all_pop_MAF <- NULL

  for (pop in pops) {

    # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
    pop_table_file <- paste(base_folder,"/",pop,".chr",chr,".tab.gz",sep="")
    # for files with not mac=1
    # pop_table_file <- paste(pop,".chr",chr,".not_fixed.not_MAC1.tab.gz",sep="")

    pop_table <- read.table(pop_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
    colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
    pop_table$DAF <- as.numeric(as.character(pop_table$DAF))

    #remove nas
    #we want to keep the NA's because those are the novel variants!
    # pop_table <- pop_table[which(!is.na(pop_table$DAF)),]
    
    #remove fixed variants (the ones with DAF == 0 or == 1) ---> those are also the ones with MAF 0, basically!
    pop_table <- pop_table[which(pop_table$DAF != 0),]
    pop_table <- pop_table[which(pop_table$DAF != 1),]

    #remove MAC = 1
    # pop_table <- pop_table[which(pop_table$MAC > 1),]

    # current_hist <- hist(pop_table$MAF,plot=F,breaks=20)

    # #use relative frequency sites per bin count/total variants number
    # current_pop_MAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pop,length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS)),chr_length=(length(pop_table$POS))))
    
    # all_pop_MAF <- rbind(all_pop_MAF,current_pop_MAF)
    all_pop_MAF <- list(all_pop_MAF,pop_table$MAF)

  }

  #################################################################TEST TO PLOT##########################
   base_folder <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/CHR",chr,sep="")
   all_pop_MAF <- NULL

    # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
    VBI_table_file <- paste(base_folder,"/","VBI",".chr",chr,".tab.gz",sep="")
    FVG_table_file <- paste(base_folder,"/","FVG",".chr",chr,".tab.gz",sep="")
    CARL_table_file <- paste(base_folder,"/","CARL",".chr",chr,".tab.gz",sep="")
    CEU_table_file <- paste(base_folder,"/","CARL",".chr",chr,".tab.gz",sep="")
    TSI_table_file <- paste(base_folder,"/","CARL",".chr",chr,".tab.gz",sep="")
    # for files with not mac=1
    # pop_table_file <- paste(pop,".chr",chr,".not_fixed.not_MAC1.tab.gz",sep="")

    VBI_pop_table <- read.table(VBI_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
    FVG_pop_table <- read.table(FVG_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
    CARL_pop_table <- read.table(CARL_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
    colnames(VBI_pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
    colnames(FVG_pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
    colnames(CARL_pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
    VBI_pop_table$DAF <- as.numeric(as.character(VBI_pop_table$DAF))
    FVG_pop_table$DAF <- as.numeric(as.character(FVG_pop_table$DAF))
    CARL_pop_table$DAF <- as.numeric(as.character(CARL_pop_table$DAF))

    #remove nas
    #we want to keep the NA's because those are the novel variants!
    # pop_table <- pop_table[which(!is.na(pop_table$DAF)),]
    
    #remove fixed variants (the ones with DAF == 0 or == 1) ---> those are also the ones with MAF 0, basically!
    VBI_pop_table <- VBI_pop_table[which(VBI_pop_table$DAF != 0),]
    FVG_pop_table <- FVG_pop_table[which(FVG_pop_table$DAF != 0),]
    CARL_pop_table <- CARL_pop_table[which(CARL_pop_table$DAF != 0),]
    VBI_pop_table <- VBI_pop_table[which(VBI_pop_table$DAF != 1),]
    FVG_pop_table <- FVG_pop_table[which(FVG_pop_table$DAF != 1),]
    CARL_pop_table <- CARL_pop_table[which(CARL_pop_table$DAF != 1),]

    #remove MAC = 1
    all_pop_MAF <- list(VBI_pop_table$MAF,FVG_pop_table$MAF,CARL_pop_table$MAF)



  require(plotrix)
  base_folder <- getwd()
  jpeg(paste(base_folder,"/1_all_pop_MAF_plotrix.jpg",sep=""),width=1800, height=800)
  multhist(all_pop_MAF, col=c("purple", "green", "red" ), freq=FALSE, ylab="relative frequency",xlab="MAF", breaks=20, ylim=c(0, 40), main="Novel Sites ") 
  legend ("topright", c("CAR", "FVG", "VBI") , ncol=3 ,  col=c("purple", "green", "red" ), pch=16 , bty='n') 
  dev.off()
  
##################################################################################################################

  all_pop_MAF$breaks <- as.numeric(as.character(all_pop_MAF$breaks))
  all_pop_MAF$counts <- as.numeric(as.character(all_pop_MAF$counts))
  all_pop_MAF$rel_count <- as.numeric(as.character(all_pop_MAF$rel_count))
  all_pop_MAF$chr_length <- as.numeric(as.character(all_pop_MAF$chr_length))
  all_pop_MAF$pop <- as.character(all_pop_MAF$pop)

  all_pop_MAF_table <- all_pop_MAF[which(all_pop_MAF$pop == "CEU"),1]
  all_pop_MAF_table_names <- c("breaks")

  for(pop_i in pops){
    actual_pop <- all_pop_MAF[which(all_pop_MAF$pop == pop_i),]
    all_pop_MAF_table_names <- c(all_pop_MAF_table_names,paste(pop_i,"count",sep="_"),paste(pop_i,"rel_count",sep="_"),paste(pop_i,"chr_length",sep="_"))

    all_pop_MAF_table <- cbind(all_pop_MAF_table,cbind(actual_pop$counts,actual_pop$rel_count*100,actual_pop$chr_length))
  }
  all_pop_MAF_table <- as.data.frame(all_pop_MAF_table)
  colnames(all_pop_MAF_table) <- all_pop_MAF_table_names

  write.table(all_pop_MAF_table,file=paste("all_pop_MAF_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
  # write.table(all_pop_MAF_table,file=paste("all_pop_MAF_count_no_MAC1_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

}
#In order to plot the genomewide spectrum, we need to upload each file and sum things together
all_chr <- NULL
for (chrom in 1:22) {
  current_file <- paste("all_pop_MAF_count_chr",chrom,".txt",sep="")
  current_chr <- read.table(current_file,sep="\t",header=T)
  #remove the relative count
  current_chr <- current_chr[,c(grep("rel_count",colnames(current_chr),invert=T))]
  # colnames(current_chr)[which(colnames(current_chr) != "breaks")]
  all_chr <- as.data.frame(c(all_chr,current_chr))
}

#remove all unnecessary breaks 
all_chr <- all_chr[,c(grep("breaks.",colnames(all_chr),invert=T))]

#now sum all together the count and the total length for each population
for (pop in pops) {
all_chr$current_pop_tot_count <- apply(all_chr[,c(grep(paste(pop,"_count",sep=""),colnames(all_chr)))],1,sum)
all_chr$current_pop_tot_length <- apply(all_chr[,c(grep(paste(pop,"_chr_length",sep=""),colnames(all_chr)))],1,sum)
colnames(all_chr)[which(colnames(all_chr)=="current_pop_tot_count")] <- paste(pop,"_tot_count",sep="")
colnames(all_chr)[which(colnames(all_chr)=="current_pop_tot_length")] <- paste(pop,"_tot_length",sep="")
}

#now we want to keep only the tot count stuff and the breaks column
all_chr <- cbind(breaks=all_chr$breaks,all_chr[,c(grep("tot_count",colnames(all_chr)))],all_chr[,c(grep("tot_length",colnames(all_chr)))])

#now we need to create the table with the relative count
#we can do it strait away
#print the table, too!!
write.table(all_chr,file=paste("all_pop_MAF_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
# all_pop_MAF_table <- cbind(breaks=all_chr$breaks,all_chr[,c(grep("tot_count",colnames(all_chr)))],(all_chr[,c(grep("tot_count",colnames(all_chr)))]/all_chr[,c(grep("tot_length",colnames(all_chr)))])*100)

#absolute count(even columns)
# barplot3 <- as.matrix(t(all_pop_MAF_table[,c(2,4,6,8,10)]))
# barplot3 <- as.matrix(t(all_pop_MAF_table[,seq(2, ncol(all_pop_MAF_table), by = 2)]))
barplot3 <- as.matrix(t(all_chr[,c(grep("tot_count",colnames(all_chr)))]))

#relative site count (odd columns)
# barplot4 <- as.matrix(t(all_pop_MAF_table[,c(3,5,7,9,11)]))
# barplot4 <- as.matrix(t(all_pop_MAF_table[,seq(3, ncol(all_pop_MAF_table), by = 2)]))
all_chr_rel <- (all_chr[,c(grep("tot_count",colnames(all_chr)))]/all_chr[,c(grep("tot_length",colnames(all_chr)))])*100)
colnames(all_chr_rel) <- gsub("tot_count","rel_freq",colnames(all_chr_rel))

barplot4 <- as.matrix(t(all_chr_rel))

all_cols <- col_pop(pops)

base_folder <-getwd()

# jpeg(paste(base_folder,"/1_all_pop_MAF_",chr,".jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder,"/1_all_pop_MAF.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_no_MAC1.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot3,
    beside=T,
    main="MAF in all populations",
    # names.arg=all_pop_MAF_table$breaks,
    names.arg=all_chr$breaks,
    xlab="",
    ylab="",col=all_cols[,1])
    legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
    mtext(1, text = "MAF", line = 4,cex=1.4)
    mtext(2, text = "Site count", line = 4,cex=1.4)
dev.off()

# jpeg(paste(base_folder,"/1_all_pop_MAF_rf_",chr,".jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder,"/1_all_pop_MAF_rf.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_rf_no_MAC1.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot4,
    beside=T,
    main="MAF in all populations",
    # names.arg=all_pop_MAF_table$breaks,
    names.arg=all_chr$breaks,
    xlab="",
    ylab="",col=all_cols[,1])
    legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
    mtext(1, text = "MAF", line = 4,cex=1.4)
    mtext(2, text = "Relative Frequency (N sites/Tot sites in freq bin)(%)", line = 4,cex=1.4)
dev.off()
ACTUALLY
################################################################################################################
#WE NEED TO UPLOAD DATA ALSO FOR NOVEL SITES from here:
base_novel_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/NOVEL"
pops_ingi_novel <- c("VBI_n","FVG_n","CARL_n")
pops_ingi <- c("VBI","FVG","CARL")

#we work by chr and by population
#upload data for each population about novel sites
for (pop in pops_ingi) {
  current_pop_novel <- NULL
  for (chr in 1:22) {
    #read the current file for this population and this chromosome
    current_chr_novel <- read.table(paste(base_novel_folder,"/",pop,".",chr,".n.maf",sep=""), sep="\t",header=F)
    colnames(current_chr_novel) <- c("CHROM","POZ","POS","TYPE","CEU","TSI","VBI","FVG","CARL")
    current_pop_novel <-rbind(current_pop_novel,current_chr_novel)
  }

 current_hist <- hist(pop_table$MAF,plot=F,breaks=20)

 #use relative frequency sites per bin count/total variants number
 current_pop_MAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pop,length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS)),chr_length=(length(pop_table$POS))))
    
 all_pop_MAF <- rbind(all_pop_MAF,current_pop_MAF)

}
################################################################################################################
#now upload the private and shared plots and add them to the previous plot
# merged_daf <- read.table(paste0(base_folder,"INPUT_FILES/INGI_chr",chr,".merged_daf.tab.gz"),skip=1,header=F)
# colnames(merged_daf) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","FVG")

# merged_daf$CEU <- as.numeric(as.character(merged_daf$CEU))
# merged_daf$TSI <- as.numeric(as.character(merged_daf$TSI))
# merged_daf$VBI <- as.numeric(as.character(merged_daf$VBI))
# merged_daf$FVG <- as.numeric(as.character(merged_daf$FVG))

#we need to use the files with MAF info, but created after filtering by DAF: we need to put all chr together
pop_VBI_private <- NULL
pop_FVG_private <- NULL
pop_CARL_private <- NULL
pop_VBI_shared <- NULL
pop_FVG_shared <- NULL
pop_CARL_shared <- NULL

for (chr in 1:22) {
# for (chr in 18:22) {
print(chr)
base_folder_res <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/CHR",chr,sep="")
current_VBI_private <- read.table(paste(base_folder_res,"/VBI.chr",chr,".tab.gz.private.tab.gz",sep=""),header=F,skip=1)
current_FVG_private <- read.table(paste(base_folder_res,"/FVG.chr",chr,".tab.gz.private.tab.gz",sep=""),header=F,skip=1)
current_CARL_private <- read.table(paste(base_folder_res,"/CARL.chr",chr,".tab.gz.private.tab.gz",sep=""),header=F,skip=1)

colnames(current_VBI_private) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
colnames(current_FVG_private) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
colnames(current_CARL_private) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  
current_VBI_private$DAF <- as.numeric(as.character(current_VBI_private$DAF))
current_VBI_private$MAF <- as.numeric(as.character(current_VBI_private$MAF))
current_VBI_private$DAC <- as.numeric(as.character(current_VBI_private$DAC))
current_VBI_private$MAC <- as.numeric(as.character(current_VBI_private$MAC))
current_FVG_private$DAF <- as.numeric(as.character(current_FVG_private$DAF))
current_FVG_private$MAF <- as.numeric(as.character(current_FVG_private$MAF))
current_FVG_private$DAC <- as.numeric(as.character(current_FVG_private$DAC))
current_FVG_private$MAC <- as.numeric(as.character(current_FVG_private$MAC))
current_CARL_private$DAF <- as.numeric(as.character(current_CARL_private$DAF))
current_CARL_private$MAF <- as.numeric(as.character(current_CARL_private$MAF))
current_CARL_private$DAC <- as.numeric(as.character(current_CARL_private$DAC))
current_CARL_private$MAC <- as.numeric(as.character(current_CARL_private$MAC))


current_VBI_shared <- read.table(paste(base_folder_res,"/VBI.chr",chr,".tab.gz.shared.tab.gz",sep=""),header=F,skip=1)
current_FVG_shared <- read.table(paste(base_folder_res,"/FVG.chr",chr,".tab.gz.shared.tab.gz",sep=""),header=F,skip=1)
current_CARL_shared <- read.table(paste(base_folder_res,"/CARL.chr",chr,".tab.gz.shared.tab.gz",sep=""),header=F,skip=1)

colnames(current_VBI_shared) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
colnames(current_FVG_shared) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
colnames(current_CARL_shared) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")

current_VBI_shared$DAF <- as.numeric(as.character(current_VBI_shared$DAF))
current_VBI_shared$MAF <- as.numeric(as.character(current_VBI_shared$MAF))
current_VBI_shared$DAC <- as.numeric(as.character(current_VBI_shared$DAC))
current_VBI_shared$MAC <- as.numeric(as.character(current_VBI_shared$MAC))
current_FVG_shared$DAF <- as.numeric(as.character(current_FVG_shared$DAF))
current_FVG_shared$MAF <- as.numeric(as.character(current_FVG_shared$MAF))
current_FVG_shared$DAC <- as.numeric(as.character(current_FVG_shared$DAC))
current_FVG_shared$MAC <- as.numeric(as.character(current_FVG_shared$MAC))
current_CARL_shared$DAF <- as.numeric(as.character(current_CARL_shared$DAF))
current_CARL_shared$MAF <- as.numeric(as.character(current_CARL_shared$MAF))
current_CARL_shared$DAC <- as.numeric(as.character(current_CARL_shared$DAC))
current_CARL_shared$MAC <- as.numeric(as.character(current_CARL_shared$MAC))

pop_VBI_private <- as.data.frame(rbind(pop_VBI_private,current_VBI_private))
pop_FVG_private <- as.data.frame(rbind(pop_FVG_private,current_FVG_private))
pop_CARL_private <- as.data.frame(rbind(pop_CARL_private,current_CARL_private))
pop_VBI_shared <- as.data.frame(rbind(pop_VBI_shared,current_VBI_shared))
pop_FVG_shared <- as.data.frame(rbind(pop_FVG_shared,current_FVG_shared))
pop_CARL_shared <- as.data.frame(rbind(pop_CARL_shared,current_CARL_shared))

}

#now we should be able to proceed as before...

#remove MAC = 1
# pop_VBI_private$MAC <- pop_VBI_private[which(pop_VBI_private$MAC > 1),]
# pop_FVG_private$MAC <- pop_FVG_private[which(pop_FVG_private$MAC > 1),]

summary(pop_VBI_private)
summary(pop_FVG_private)
summary(pop_CARL_private)
# summary(pop_VBI_shared)
# summary(pop_FVG_shared)

dim(pop_VBI_private)
dim(pop_FVG_private)
dim(pop_CARL_private)
# dim(pop_VBI_shared)
# dim(pop_FVG_shared)

# after removal of MAC=1
# VBI_private : 21170
# FVG_private : 22631
# VBI_shared : 109341
# FVG_shared : 104707


#this is useless now...
# pop_VBI_private <- pop_VBI_private[which(pop_VBI_private$DAF != 1),]
# pop_FVG_private <- pop_FVG_private[which(pop_FVG_private$DAF != 1),]
# pop_VBI_shared <- pop_VBI_shared[which(pop_VBI_shared$DAF != 1),]
# pop_FVG_shared <- pop_FVG_shared[which(pop_FVG_shared$DAF != 1),]

pops_ingi <- c("VBI_p","FVG_p","CARL_p","VBI_s","FVG_s","CARL_s")
# pops_ingi <- c("VBI_p","FVG_p")
pops_ingi_files <- c("pop_VBI_private","pop_FVG_private","pop_CARL_private","pop_VBI_shared","pop_FVG_shared","pop_CARL_shared")
# pops_ingi_files <- c("pop_VBI_private","pop_FVG_private")

all_pop_MAF_sp <- NULL

for (k in 1:length(pops_ingi_files)){
  current_pop <- get(pops_ingi_files[k])

  current_hist <- hist(current_pop$MAF,plot=F, breaks=20)

  #added relative frequency sites per bin count/total variants number for the current category! shared or private!!
  current_pop_MAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pops_ingi[k],length(current_hist$counts)),rel_count=(current_hist$counts/length(current_pop$POS))))
  current_pop_MAF$breaks <- as.numeric(as.character(current_pop_MAF$breaks))
  current_pop_MAF$counts <- as.numeric(as.character(current_pop_MAF$counts))
  current_pop_MAF$rel_count <- as.numeric(as.character(current_pop_MAF$rel_count))
  current_pop_MAF$pop <- as.character(current_pop_MAF$pop)

  #add again to all_pop_DAF to have a complete plot
  # all_pop_MAF <- rbind(all_pop_MAF,current_pop_MAF)
  #but now we want to plot shared and private separately from all the others..
  all_pop_MAF_sp <- rbind(all_pop_MAF_sp,current_pop_MAF)
  
}

#overwrite all_pop_MAF with all_pop_MAF_sp in order to keep the code clean...
all_pop_MAF <- all_pop_MAF_sp
# all_pop_MAF_table_sp <- all_pop_MAF[which(all_pop_MAF$pop == "CEU"),1]
all_pop_MAF_table_sp <- all_pop_MAF[which(all_pop_MAF$pop == "VBI_p"),1]
all_pop_MAF_table_sp_names <- c("breaks")

#let's work only with ingi data at first
# pops2 <- sort(c(pops,pops_ingi))
pops2 <- sort(pops_ingi)

for(pop_i2 in pops2){
  actual_pop2 <- all_pop_MAF[which(all_pop_MAF$pop == pop_i2),]
  all_pop_MAF_table_sp_names <- c(all_pop_MAF_table_sp_names,paste(pop_i2,"count",sep="_"),paste(pop_i2,"rel_count",sep="_"))

  all_pop_MAF_table_sp <- cbind(all_pop_MAF_table_sp,cbind(actual_pop2$counts,actual_pop2$rel_count*100))
}

all_pop_MAF_table_sp <- as.data.frame(all_pop_MAF_table_sp)
colnames(all_pop_MAF_table_sp) <- all_pop_MAF_table_sp_names

# write.table(all_pop_MAF_table_sp,file=paste("All_pop_MAF_sp_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
# without MAC=1
# write.table(all_pop_MAF_table_sp,file=paste("All_pop_MAF_sp_count_no_MAC1_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

# now upload data from all pop maf
all_pop_MAF_old <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/all_pop_MAF_count.txt",sep="\t",header=T)
all_pop_MAF_old[,c(grep("tot_count",colnames(all_pop_MAF_old)))]

#absolute count
# barplot5 <- as.matrix(t(all_pop_MAF_table_sp[,c(2,4,6,8,10,12,14,16)]))
barplot5 <- as.matrix(t(cbind(all_pop_MAF_table_sp[,seq(2, ncol(all_pop_MAF_table_sp), by = 2)],all_pop_MAF_old[,c(grep("tot_count",colnames(all_pop_MAF_old)))])))

#relative site count
# barplot6 <- as.matrix(t(all_pop_MAF_table_sp[,c(3,5,7,9,11,13,15,17)]))
barplot6 <- as.matrix(t(cbind(all_pop_MAF_table_sp[,seq(3, ncol(all_pop_MAF_table_sp), by = 2)],(all_pop_MAF_old[,c(grep("tot_length",colnames(all_pop_MAF_old)))]/all_pop_MAF_old[,c(grep("tot_count",colnames(all_pop_MAF_old)))])*100)))
# barplot6 <- as.matrix(t(all_pop_MAF_table_sp[,seq(3, ncol(all_pop_MAF_table_sp), by = 2)]))

#relative count but only for private sites
# barplot7 <- as.matrix(t(all_pop_MAF_table_sp[,c(3,5,9,11,13,15)]))
barplot7 <- as.matrix(t(cbind(all_pop_MAF_table_sp[,c(3,7,11)],(all_pop_MAF_old[,c(grep("tot_length",colnames(all_pop_MAF_old)))]/all_pop_MAF_old[,c(grep("tot_count",colnames(all_pop_MAF_old)))]))))

all_cols <- NULL
pops2 <- sort(c(pops,pops_ingi))

for(i in 1:length(pops2)){
  if(pops2[i] == "CEU"){
    cur_col <- colors()[130]
  }
  if(pops2[i] == "FVG"){
    cur_col <- colors()[517]
  }
  if(pops2[i] == "FVG_p"){
    cur_col <- colors()[50]
  }
  if(pops2[i] == "FVG_s"){
    cur_col <- colors()[81]
  }
  if(pops2[i] == "TSI"){
    cur_col <- colors()[30]
  }
  if(pops2[i] == "VBI"){
    cur_col <- colors()[421]
  }
  if(pops2[i] == "VBI_p"){
    cur_col <- colors()[33]
  }
  if(pops2[i] == "VBI_s"){
    cur_col <- colors()[36]
  }
    if(pops2[i] == "CARL"){
    cur_col <- colors()[95]
  }
  if(pops2[i] == "CARL_p"){
    cur_col <- colors()[98]
  }
  if(pops2[i] == "CARL_s"){
    cur_col <- colors()[99]
  }
  pop_col <- cbind(cur_col,pops2[i])
  all_cols <- rbind(all_cols,pop_col)
}

# no MAC=1
base_folder_res <- getwd()
# jpeg(paste(base_folder_res,"/1_all_pop_MAF_sp_",chr,".jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder_res,"/1_all_pop_MAF_sp.jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder_res,"/1_all_5pop_MAF_sp.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder_res,"/1_all_pop_MAF_sp_no_MAC1.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_sp.jpg",sep=""),width=1800, height=800)
    par(oma=c(3,3,3,3),cex=1.4)
    barplot(barplot5,
      beside=T,
      main="MAF in all populations",
      names.arg=all_pop_MAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[,1]
      )
      mtext(1, text = "MAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Site count", line = 4,cex=1.4)
      legend(x="top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=8)
dev.off()

# jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_sp_rf.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder_res,"/1_all_pop_MAF_sp_rf_",chr,".jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder_res,"/1_all_pop_MAF_sp_rf.jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder_res,"/1_all_5pop_MAF_sp_rf.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder_res,"/1_all_pop_MAF_sp_rf_no_MAC1.jpg",sep=""),width=1800, height=800)
    par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot6,
    beside=T,
    main="MAF in all populations",
    names.arg=all_pop_MAF_table_sp$breaks*100,
    xlab="",
    ylab="",col=all_cols[,1],ylim=c(0,60))
    mtext(1, text = "MAF (%)", line = 4,cex=1.4)
    mtext(2, text = "Relative Frequency (N sites/Tot sites per category)(%)", line = 4,cex=1.4)
    legend("top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=8)
dev.off()

#plot only private, not shared
# jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_p_rf.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder_res,"/1_all_pop_MAF_p_rf.jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder_res,"/1_all_pop_3MAF_p_rf.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_p_rf_no_MAC1.jpg",sep=""),width=1800, height=800)
    par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot7,
    beside=T,
    main="MAF in all populations",
    names.arg=all_pop_MAF_table_sp$breaks*100,
    xlab="",
    # ylab="",col=all_cols[c(1,2,3,5,6,7),1],ylim=c(0,85))
    # ylab="",col=all_cols[c(1,3,5),1],ylim=c(0,85))
    ylab="",col=all_cols[c(2,6,10,6,8,9,5,1),1],ylim=c(0,85))
    mtext(1, text = "MAF (%)", line = 4,cex=1.4)
    mtext(2, text = "Relative Frequency (N sites/Tot sites per category)(%)", line = 4,cex=1.4)
    # legend("top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[c(1,2,3,5,6,7),1],legend=all_cols[c(1,2,3,5,6,7),2], ncol=6)
    legend("top",pch =c(rep(19,length(all_cols[c(2,6,10,6,8,9,5,1),1]))),col=all_cols[c(2,6,10,6,8,9,5,1),1],legend=all_cols[c(2,6,10,6,8,9,5,1),2], ncol=6)
    # legend("top",pch =c(rep(19,3)),col=all_cols[c(1,3,5),1],legend=all_cols[c(1,3,5),2], ncol=6)
dev.off()


###########################################################################################
#Plot 2: venn diagram with overlap between all populations and categories

require(gplots)
require(VennDiagram)

#We need to get everything that has DAF > 0 for each population, this means that in each row we want to keep stuff if there is at least a DAF > 0
#But since we're trying to show how ALL variants split in the 4 populations, we NEED TO DO NOTHING HERE!!!WE NEED ALL VARIANTS!!!

# #we are going to remove first the fixed sites from isolate populations
# merged_daf_diag_not_fixed_iso <- merged_daf[which(merged_daf$FVG != 1 & merged_daf$VBI != 1 ),]
# #now we'll remove from this dataset all the fixed sites with DAF == 0
# merged_daf_diag_not_fixed_iso <- merged_daf_diag_not_fixed_iso[which(merged_daf_diag_not_fixed_iso$FVG != 0 | merged_daf_diag_not_fixed_iso$VBI != 0 ),]

# merged_daf_diag_not_fixed <- merged_daf[which(merged_daf$CEU != 1 & merged_daf$FVG != 1 & merged_daf$TSI != 1 & merged_daf$VBI != 1 ),]
# merged_daf_diag_not_fixed <- merged_daf_diag_not_fixed[which(merged_daf_diag_not_fixed$CEU != 0 | merged_daf_diag_not_fixed$FVG != 0 | merged_daf_diag_not_fixed$TSI != 0 | merged_daf_diag_not_fixed$VBI != 0 ),]

# merged_daf_diag <- merged_daf_diag_not_fixed[,-c(1:4)]

# merged_daf_diag <- merged_daf_diag_not_fixed_iso[,-c(1:4)]
merged_daf_diag <- merged_daf[,-c(1:4)]

merged_daf_diag[merged_daf_diag>0]=1
aree=colSums(merged_daf_diag,na.rm=T)

# #easy venn diagramm EVER
# res.all <- list(which(merged_daf_diag[,1]>0),
#   which(merged_daf_diag[,2]>0),
#   which(merged_daf_diag[,3]>0),
#   which(merged_daf_diag[,4]>0))
# names(res.all) <- colnames(merged_daf_diag)

#we need to do all the cases for intersections
res.all2 <- NULL
for(i in c(0,1)){
  for(j in c(0,1)){
    for(k in c(0,1)){
      for(h in c(0,1)){
        if(any(c(i,j,k,h)==1)){
          res=length(which(merged_daf_diag$CEU==i & merged_daf_diag$TSI==j & merged_daf_diag$VBI==k & merged_daf_diag$FVG==h))
          res.all2=rbind(res.all2,c(res,paste(which(c(i,j,k,h)==1),collapse="")))
        }
      }
    }
  }
}

res.all3 <- as.data.frame(res.all2[order(res.all2[,2]),])
res.all3$V1 <- as.numeric(as.character(res.all3$V1))
res.all3$V2 <- as.character(res.all3$V2)
res.all3$V3 <- as.character(res.all3$V2)
colnames(res.all3) <- c("n_sites","cat_code","pop")

for(i in 1:length(res.all3$cat_code)){
  for(j in names(aree)){
    cat <- res.all3$pop[i]
    newcat <- sub(as.character(which(names(aree)==j)),paste(j,";",sep=""),cat)
    res.all3$pop[i] <- newcat
  }
}

all_cols <- NULL
cats <- colnames(merged_daf_diag)

for(i in 1:length(cats)){
  if(cats[i] == "CEU"){
    cur_col <- colors()[130]
  }
  if(cats[i] == "FVG"){
    cur_col <- colors()[517]
  }
  if(cats[i] == "FVG_p"){
    cur_col <- colors()[50]
  }
  if(cats[i] == "FVG_s"){
    cur_col <- colors()[81]
  }
  if(cats[i] == "TSI"){
    cur_col <- colors()[30]
  }
  if(cats[i] == "VBI"){
    cur_col <- colors()[421]
  }
  if(cats[i] == "VBI_p"){
    cur_col <- colors()[33]
  }
  if(cats[i] == "VBI_s"){
    cur_col <- colors()[36]
  }
  
  pop_col <- cbind(cur_col,cats[i])
  all_cols <- rbind(all_cols,pop_col)
}

#write the table:
write.table(res.all3,file=paste(base_folder,"RESULTS/VENN/All_pop_intersect_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#find all indexes
idx1=grep("1",res.all2[,2])
idx2=grep("2",res.all2[,2])
idx3=grep("3",res.all2[,2])
idx4=grep("4",res.all2[,2])

v <- draw.quad.venn(area1=as.numeric(sum(as.numeric(res.all2[idx1,1]))),
  area2=as.numeric(sum(as.numeric(res.all2[idx2,1]))),
  area3=as.numeric(sum(as.numeric(res.all2[idx3,1]))),
  area4=as.numeric(sum(as.numeric(res.all2[idx4,1]))),
  n12=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx2),1]))),
  n13=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx3),1]))),
  n14=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx4),1]))),
  n23=as.numeric(sum(as.numeric(res.all2[intersect(idx2,idx3),1]))),
  n24=as.numeric(sum(as.numeric(res.all2[intersect(idx2,idx4),1]))),
  n34=as.numeric(sum(as.numeric(res.all2[intersect(idx3,idx4),1]))),
  n123=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx2),idx3),1]))),
  n124=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx2),idx4),1]))),
  n134=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx3),idx4),1]))),
  n234=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx2,idx3),idx4),1]))),
  n1234=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(idx1,idx2),idx3),idx4),1]))),
  ind=FALSE, category=colnames(merged_daf_diag),fill=all_cols[,1])

jpeg(paste(base_folder,"PLOTS/2_all_pop_DAF_VENN.jpg",sep=""),width=800, height=800,pointsize = 20)
  grid.draw(v)
dev.off()


###########################################################################################
#Plot 3 barplot for privte and shared variants for each INGI population by category
#use also here the relative frequency

funk_cat <-read.csv(paste(base_folder,"INPUT_FILES/consequences.list",sep=""),header=F)
funk_cat$V1 <- as.character(funk_cat$V1)

chr <- "22"
pop_is <- c("VBI","FVG")

for (pop in pop_is) {

  conseq_count_tot <- NULL

  pop_table_file <- paste(base_folder,"INPUT_FILES/",pop,".chr",chr,".tab.gz",sep="")

  pop_table <- read.table(pop_table_file,header=F,skip=1)
  colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
  current_private <- get(paste("pop",pop,"private",sep="_"))
  current_shared <- get(paste("pop",pop,"shared",sep="_"))

  for(cat in funk_cat$V1){

    current_grep <- pop_table[(grep(cat,pop_table$INFO)),]

    #merge the grepped data with the private
    current_grep_private <- merge(current_private,current_grep,by.x="POS",by.y="POS")

    #merge the grepped data with the shared
    current_grep_shared <- merge(current_shared,current_grep,by.x="POS",by.y="POS")

    current_conseq <- cbind(category=as.character(cat),private=length(current_grep_private$POS),shared=length(current_grep_shared$POS))
    conseq_count_tot <- rbind(conseq_count_tot,current_conseq)

  }

  conseq_count_tot <- as.data.frame(conseq_count_tot)
  conseq_count_tot$category <- as.character(conseq_count_tot$category)
  conseq_count_tot$private <- as.numeric(as.character(conseq_count_tot$private))
  conseq_count_tot$shared <- as.numeric(as.character(conseq_count_tot$shared))
  conseq_count_tot$N_shared <- rep(length(current_shared$POS),length(conseq_count_tot$category))
  conseq_count_tot$N_private <- rep(length(current_private$POS),length(conseq_count_tot$category))
  conseq_count_tot$rel_private <- (conseq_count_tot$private/conseq_count_tot$N_private)*100
  conseq_count_tot$rel_shared <- (conseq_count_tot$shared/conseq_count_tot$N_shared)*100

  #lets check if the differences are significative with a chisquare test
  for(i in 1:length(conseq_count_tot$category)){
    print(conseq_count_tot$category[i])
    A=matrix(c(conseq_count_tot$shared[i], conseq_count_tot$N_shared[i], conseq_count_tot$private[i], conseq_count_tot$N_private[i]), nrow=2, ncol=2, byrow=T)
    pippo=chisq.test(A)
    conseq_count_tot$p[i] <- pippo$p.value
  }


  write.table(conseq_count_tot,file=paste(base_folder,"RESULTS/CONSEQUENCES/",pop,"_consequences_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
  
  #site count
  barplot1 <- as.matrix(t(conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),c(2:3)]))

  #relative site count
  barplot2 <- as.matrix(t(conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),c(6:7)]))

  jpeg(paste(base_folder,"PLOTS/3_",pop,"_all_conseq_DAF.jpg",sep=""),width=1800, height=800)
    par(oma=c(9,3,3,3),cex=1.4)
    barplot(barplot1,
      beside=T,
      col=c("red","blue"),
      legend=c(paste("Private (tot:",conseq_count_tot$N_private[1],")",sep=""),paste("Shared (tot:",conseq_count_tot$N_shared[1],")",sep="")),
      main=paste(pop," DAF in private/shared sites categories for different consequences annotation classes",sep=""),
      names.arg=conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),1],
      xlab="",las=2,
      ylab="")
      mtext(1, text = "Annotations", line = 11)
      mtext(2, text = "Frequency", line = 4)
  dev.off()

  jpeg(paste(base_folder,"PLOTS/3_",pop,"_all_conseq_DAF_rf.jpg",sep=""),width=1800, height=800)
    par(oma=c(9,3,3,3),cex=1.4)
    barplot(barplot2,
      beside=T,
      col=c("red","blue"),
      legend=c("Private","Shared"),
      main=paste(pop," DAF in private/shared sites categories for different consequences annotation classes",sep=""),
      names.arg=conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),1],
      xlab="",las=2,
      ylab="")
      mtext(1, text = "Annotations", line = 11)
      mtext(2, text = "Relative Frequency (N sites/Tot sites in category)(%)", line = 4)
  dev.off()
}

###########################################################################################
#Plot 4 For each population (VBI FVG)  density distribution + wilcoxon test  (shared e private)  GERP SCORE  [stratify for  genic/intergenic? functional categories?]
# We'll plot the gerp score distribution for shared and private
pop_is <- c("VBI","FVG")

for (pop in pop_is) {

  conseq_count_tot <- NULL

  pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/INPUT_FILES/",pop,".chr",chr,".tab",sep="")

  pop_table <- read.table(pop_table_file,header=T)
  # colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
  current_private <- get(paste("pop",pop,"private",sep="_"))
  current_shared <- get(paste("pop",pop,"shared",sep="_"))

  #we need to get the information of the gerp score
  current_grep_GERP <- pop_table[(grep("GERP",pop_table$INFO)),]
  
  #merge the grepped data with the private
  current_grep_GERP_private <- merge(current_private,current_grep_GERP,by.x="POS",by.y="POS")
  
  #merge the grepped data with the shared
  current_grep_GERP_shared <- merge(current_shared,current_grep_GERP,by.x="POS",by.y="POS")
  
  current_conseq <- cbind(category=as.character(cat),private=length(current_grep_private$POS),shared=length(current_grep_shared$POS))
  conseq_count_tot <- rbind(conseq_count_tot,current_conseq)

}

###########################################################################################
#Plot 5 SIFT Polyphen  https://faculty.washington.edu/wjs18/GS561/cSNPs_lab.html     [stratify for  genic/intergenic? functional categories?]

###########################################################################################


###########################################################################################
#Plot 6 ROH
chr <- "10"
pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")

# input_format <- "PLINK"
input_format <- "BEAGLE"
xmax <- NULL
# LOD <- 4
LOD <- 5
# base_folder <- getwd() #only if we're lazy
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    # pop_roh_file <- paste(pop,"roh.length.hbd",sep=".")
    # pop_roh_file <- paste(pop,"roh.length.4.hbd",sep=".")
    pop_roh_file <- paste("CHR",chr,"/",pop,".roh.length.",LOD,".hbd",sep="")
    pop_roh_table <- read.table(pop_roh_file,header=FALSE,sep=" ")
    colnames(pop_roh_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_roh_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_roh_table <- read.table(pop_roh_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_roh_table <- pop_roh_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_roh_table$IID2 <- as.character(pop_roh_table$IID2)
  pop_roh_table$IID1 <- as.character(pop_roh_table$IID1)
  current_pop_table_name <- paste(pop,"roh",sep="_")
  assign(current_pop_table_name,pop_roh_table)

  tot_roh <- tapply(pop_roh_table$LENGTH,pop_roh_table$IID1,sum)
  tot_roh <- as.data.frame(tot_roh, row.names=NULL)
  tot_roh$ID <- rownames(tot_roh)
  colnames(tot_roh) <- c("ROH_tot","ID")
  tot_roh$ID <- as.character(tot_roh$ID)
  tot_roh$ROH_tot <- as.numeric(as.character(tot_roh$ROH_tot))
  tot_roh$ROH_tot <- tot_roh$ROH_tot/1000000

  assign(paste(pop,"tot_roh",sep="_"),tot_roh)
  
  assign(paste("M",pop,sep="_"),ecdf(tot_roh$ROH_tot))
  xmax <- c(xmax,summary(ecdf(tot_roh$ROH_tot))[6])
  # jpeg(paste("6_roh_",pop".jpg",sep=""),width=800, height=800)
  # plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity")
  # lines(M_FVG,FVG_tot_roh$ROH_tot,col="red")
  # lines(M_TSI,TSI_tot_roh$ROH_tot,col="green")
  # lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue")
  # legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
  # dev.off()

}

q_CEU <- quantile(CEU_tot_roh$ROH_tot,c(.95))
q_TSI <- quantile(TSI_tot_roh$ROH_tot,c(.95))
q_FVG <- quantile(FVG_tot_roh$ROH_tot,c(.95))
q_VBI <- quantile(VBI_tot_roh$ROH_tot,c(.95))
q_CARL <- quantile(CARL_tot_roh$ROH_tot,c(.95))
q_Sauris <- quantile(Sauris_tot_roh$ROH_tot,c(.95))
q_Erto <- quantile(Erto_tot_roh$ROH_tot,c(.95))
q_Illegio <- quantile(Illegio_tot_roh$ROH_tot,c(.95))
q_Resia <- quantile(Resia_tot_roh$ROH_tot,c(.95))

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.335094 
# TSI: 95% -> 1.018848 
# FVG: 95% -> 4.769859 
# VBI: 95% -> 1.972573 

# W50 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU : 95% -> 1.253247
# TSI : 95% -> 1.147975
# FVG : 95% -> 4.582939
# VBI : 95% -> 2.059055

# W50K O7K (IBDseq, filtered at LOD = 5)
# q_CEU : 95% -> 0.4587313 
# q_TSI : 95% -> 0.9422786 
# q_FVG : 95% -> 7.134104 
# q_VBI : 95% -> 6.582102

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# q_CEU: 95% -> 1.269719
# q_TSI: 95% -> 1.138477
# q_FVG: 95% -> 4.863319
# q_VBI: 95% -> 2.069055
# q_CARL: 95% -> 2.658259

#W24 (Beagle RefinedIBD, filtered at LOD = 4)
# q_CEU: 95% -> 0.191057
# q_TSI: 95% -> 0.2552902
# q_FVG: 95% -> 1.918711
# q_VBI: 95% -> 0.7033285
# q_CARL: 95% -> 0.4669344

#W24 (Beagle RefinedIBD, filtered at LOD = 4) 01/08/2014
# q_CEU: 95% -> 4.118258 
# q_TSI: 95% -> 3.991706 
# q_FVG: 95% -> 15.47934 
# q_VBI: 95% -> 5.637363 
# q_CARL: 95% -> 7.155192 
# q_Sauris: 95% -> 13.92794 
# q_Erto: 95% -> 10.56745 
# q_Illegio: 95% -> 17.53563 
# q_Resia: 95% -> 23.23969

#W24 (Beagle RefinedIBD, filtered at LOD = 5) 01/08/2014
# q_CEU: 95% -> 3.407679 
# q_TSI: 95% -> 3.136137 
# q_FVG: 95% -> 15.98222 
# q_VBI: 95% -> 5.134627 
# q_CARL: 95% -> 6.092222 
# q_Sauris: 95% -> 13.14537 
# q_Erto: 95% -> 10.32173 
# q_Illegio: 95% -> 16.37155 
# q_Resia: 95% -> 22.91865 

# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all_IBDseq.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all.jpg",sep=""),width=800, height=800)
base_folder <- getwd()

jpeg(paste(base_folder,"/6_roh_all_5POP",chr,".jpg",sep=""),width=1000, height=1000)
jpeg(paste(base_folder,"/6_roh_all_5POP_lod5_",chr,".jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=1.5)
  plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  lines(M_TSI,TSI_tot_roh$ROH_tot,col="green", verticals=TRUE, pch=46)
  lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue", verticals=TRUE, pch=46)
  lines(M_FVG,FVG_tot_roh$ROH_tot,col="red", verticals=TRUE, pch=46)
  lines(M_CARL,CARL_tot_roh$ROH_tot,col="yellow", verticals=TRUE, pch=46)
  lines(M_Erto,Erto_tot_roh$ROH_tot,col="darkred", verticals=TRUE, pch=46)
  lines(M_Illegio,Illegio_tot_roh$ROH_tot,col="indianred", verticals=TRUE, pch=46)
  lines(M_Resia,Resia_tot_roh$ROH_tot,col="mediumvioletred", verticals=TRUE, pch=46)
  lines(M_Sauris,Sauris_tot_roh$ROH_tot,col="orangered", verticals=TRUE, pch=46)
  legend("bottomright",pch =c(rep(19,length(pops))),legend=pops,col=c("black","green","blue","red","yellow","darkred","indianred","mediumvioletred","orangered"),ncol=length(pops))
dev.off()

# jpeg(paste(base_folder,"PLOTS/6_roh_iso_IBDseq.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_iso.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"/6_roh_iso_3POP.jpg",sep=""),width=800, height=800)
  par(lwd=4,cex=2)
  plot(M_FVG,FVG_tot_roh$ROH_tot,verticals=TRUE, pch=46,main="",xlab="Total ROH homozigosity (Mb)",col="red")
  lines(M_VBI,VBI_tot_roh$ROH_tot,verticals=TRUE, pch=46,col="blue")
  lines(M_CARL,CARL_tot_roh$ROH_tot,verticals=TRUE, pch=46,col="yellow")
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("FVG","VBI","CARL"),col=c("red","blue","yellow"),ncol=3)
dev.off()

# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/7_roh.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"PLOTS/7_roh_density.jpg",sep=""),width=800, height=800)
  plot(density(VBI_tot_roh$ROH_tot),main="",xlab="Total ROH homozigosity", col="blue")
  lines(density(FVG_tot_roh$ROH_tot),col="red")
  lines(density(TSI_tot_roh$ROH_tot),col="green")
  lines(density(CEU_tot_roh$ROH_tot))
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()

###########################################################################################
#Plot 6_A_ROH GENOMEWIDE BY POPULATION
###########################################################################################

pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")

# input_format <- "PLINK"
input_format <- "BEAGLE"
xmax <- NULL
# LOD <- 4
LOD <- 5
# base_folder <- getwd() #only if we're lazy
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    # pop_roh_file <- paste(pop,"roh.length.hbd",sep=".")
    # pop_roh_file <- paste(pop,"roh.length.4.hbd",sep=".")
    pop_roh_file <- paste(pop,"/",pop,".WG.roh.length.",LOD,".hbd",sep="")
    pop_roh_table <- read.table(pop_roh_file,header=FALSE,sep=" ")
    colnames(pop_roh_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_roh_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_roh_table <- read.table(pop_roh_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_roh_table <- pop_roh_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_roh_table$IID2 <- as.character(pop_roh_table$IID2)
  pop_roh_table$IID1 <- as.character(pop_roh_table$IID1)
  current_pop_table_name <- paste(pop,"roh",sep="_")
  assign(current_pop_table_name,pop_roh_table)

  tot_roh <- tapply(pop_roh_table$LENGTH,pop_roh_table$IID1,sum)
  tot_roh <- as.data.frame(tot_roh, row.names=NULL)
  tot_roh$ID <- rownames(tot_roh)
  colnames(tot_roh) <- c("ROH_tot","ID")
  tot_roh$ID <- as.character(tot_roh$ID)
  tot_roh$ROH_tot <- as.numeric(as.character(tot_roh$ROH_tot))
  tot_roh$ROH_tot <- tot_roh$ROH_tot/1000000

  assign(paste(pop,"tot_roh",sep="_"),tot_roh)
  
  assign(paste("M",pop,sep="_"),ecdf(tot_roh$ROH_tot))
  xmax <- c(xmax,summary(ecdf(tot_roh$ROH_tot))[6])
  # jpeg(paste("6_roh_",pop".jpg",sep=""),width=800, height=800)
  # plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity")
  # lines(M_FVG,FVG_tot_roh$ROH_tot,col="red")
  # lines(M_TSI,TSI_tot_roh$ROH_tot,col="green")
  # lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue")
  # legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
  # dev.off()

}

q_CEU <- quantile(CEU_tot_roh$ROH_tot,c(.95))
q_TSI <- quantile(TSI_tot_roh$ROH_tot,c(.95))
q_FVG <- quantile(FVG_tot_roh$ROH_tot,c(.95))
q_VBI <- quantile(VBI_tot_roh$ROH_tot,c(.95))
q_CARL <- quantile(CARL_tot_roh$ROH_tot,c(.95))
q_Sauris <- quantile(Sauris_tot_roh$ROH_tot,c(.95))
q_Erto <- quantile(Erto_tot_roh$ROH_tot,c(.95))
q_Illegio <- quantile(Illegio_tot_roh$ROH_tot,c(.95))
q_Resia <- quantile(Resia_tot_roh$ROH_tot,c(.95))

summary_CEU <- summary(CEU_tot_roh$ROH_tot)
summary_TSI <- summary(TSI_tot_roh$ROH_tot)
summary_FVG <- summary(FVG_tot_roh$ROH_tot)
summary_VBI <- summary(VBI_tot_roh$ROH_tot)
summary_CARL <- summary(CARL_tot_roh$ROH_tot)
summary_Sauris <- summary(Sauris_tot_roh$ROH_tot)
summary_Erto <- summary(Erto_tot_roh$ROH_tot)
summary_Illegio <- summary(Illegio_tot_roh$ROH_tot)
summary_Resia <- summary(Resia_tot_roh$ROH_tot)

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.335094 
# TSI: 95% -> 1.018848 
# FVG: 95% -> 4.769859 
# VBI: 95% -> 1.972573 

# W50 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU : 95% -> 1.253247
# TSI : 95% -> 1.147975
# FVG : 95% -> 4.582939
# VBI : 95% -> 2.059055

# W50K O7K (IBDseq, filtered at LOD = 5)
# q_CEU : 95% -> 0.4587313 
# q_TSI : 95% -> 0.9422786 
# q_FVG : 95% -> 7.134104 
# q_VBI : 95% -> 6.582102

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# q_CEU: 95% -> 1.269719
# q_TSI: 95% -> 1.138477
# q_FVG: 95% -> 4.863319
# q_VBI: 95% -> 2.069055
# q_CARL: 95% -> 2.658259

#W24 (Beagle RefinedIBD, filtered at LOD = 4)
# q_CEU: 95% -> 0.191057
# q_TSI: 95% -> 0.2552902
# q_FVG: 95% -> 1.918711
# q_VBI: 95% -> 0.7033285
# q_CARL: 95% -> 0.4669344

#W24 (Beagle RefinedIBD, filtered at LOD = 4) 01/08/2014
# q_CEU: 95% -> 4.118258 
# q_TSI: 95% -> 3.991706 
# q_FVG: 95% -> 15.47934 
# q_VBI: 95% -> 5.637363 
# q_CARL: 95% -> 7.155192 
# q_Sauris: 95% -> 13.92794 
# q_Erto: 95% -> 10.56745 
# q_Illegio: 95% -> 17.53563 
# q_Resia: 95% -> 23.23969

#W24 (Beagle RefinedIBD, filtered at LOD = 5) 01/08/2014
# q_CEU: 95% -> 3.407679 
# q_TSI: 95% -> 3.136137 
# q_FVG: 95% -> 15.98222 
# q_VBI: 95% -> 5.134627 
# q_CARL: 95% -> 6.092222 
# q_Sauris: 95% -> 13.14537 
# q_Erto: 95% -> 10.32173 
# q_Illegio: 95% -> 16.37155 
# q_Resia: 95% -> 22.91865 

#W24 (Beagle RefinedIBD, filtered at LOD = 5) 01/08/2014 GENOMEWIDE per pop
# q_CEU: 95% -> 53.80378 
# q_TSI: 95% -> 46.95559 
# q_FVG: 95% -> 143.3815 
# q_VBI: 95% -> 74.93171 
# q_CARL: 95% -> 79.10135 
# q_Sauris: 95% -> 186.39 
# q_Erto: 95% -> 174.3446 
# q_Illegio: 95% -> 210.5669 
# q_Resia: 95% -> 188.7471 

#W24 (Beagle RefinedIBD, filtered at LOD = 5) ROH 05/08/2014 GENOMEWIDE per pop (subsampled at 46 individuals)
# q_CEU: 95% -> 99.67381 
# q_TSI: 95% -> 94.80848 
# q_FVG: 95% -> 235.589 
# q_VBI: 95% -> 189.19 
# q_CARL: 95% -> 118.7458 
# q_Sauris: 95% -> 211.1938 
# q_Erto: 95% -> 163.8415 
# q_Illegio: 95% -> 209.7505 
# q_Resia: 95% -> 206.2326 

# ROH           
# POP Min.  1st_Qu. Median  Mean  3rd_Qu. Max.
# summary_CEU 62.32 77.84 83.75 85.06 89.52 169.2
# summary_TSI 62.41 74.63 79.08 80.66 84.92 120.5
# summary_FVG 24.97 125.9 150.3 152.8 174.1 294.2
# summary_VBI 14.74 85.72 94.04 110.5 114.1 404.5
# summary_CARL  30.14 71.37 83.46 83.92 93.15 126
# summary_Sauris  77.59 101 139.9 136.3 159.7 231.8
# summary_Erto  66.56 101.1 118.4 121.4 141.9 212.2
# summary_Illegio 17.23 100.9 121.3 126.1 153.8 253
# summary_Resia 19.13 127.6 147.1 147.4 172.5 229.7



# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all_IBDseq.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"PLOTS/6_roh_all.jpg",sep=""),width=800, height=800)
pop_colors <- col_pop(pops)
base_folder <- getwd()

# jpeg(paste(base_folder,"/6_roh_all_5POP",chr,".jpg",sep=""),width=1000, height=1000)
jpeg(paste(base_folder,"/6_A_roh_all_5POP_lod5.WG.jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=1.5)
  # plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  plot(M_CEU,CEU_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "CEU"),1],main="",xlab="Total ROH homozigosity (Mb)", xlim=c(0,450), verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_TSI,TSI_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "TSI"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_VBI,VBI_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "VBI"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_FVG,FVG_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "FVG"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_CARL,CARL_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "CARL"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Erto,Erto_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Erto"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Illegio,Illegio_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Illegio"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Resia,Resia_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Resia"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  lines(M_Sauris,Sauris_tot_roh$ROH_tot,col=pop_colors[which(pop_colors$pop == "Sauris"),1], verticals=TRUE, pch=46,yaxs='i',col.01line='black')
  abline(h=0.95,col='grey',lty='dashed')
  #define parameters for legend
  
  leg_txt <- c(pop_colors[which(pop_colors$pop == "CEU"),2],pop_colors[which(pop_colors$pop == "TSI"),2],pop_colors[which(pop_colors$pop == "VBI"),2],pop_colors[which(pop_colors$pop == "FVG"),2],pop_colors[which(pop_colors$pop == "CARL"),2],pop_colors[which(pop_colors$pop == "Erto"),2],pop_colors[which(pop_colors$pop == "Illegio"),2],pop_colors[which(pop_colors$pop == "Resia"),2],pop_colors[which(pop_colors$pop == "Sauris"),2])
  bkg <- c(pop_colors[which(pop_colors$pop == "CEU"),1],pop_colors[which(pop_colors$pop == "TSI"),1],pop_colors[which(pop_colors$pop == "VBI"),1],pop_colors[which(pop_colors$pop == "FVG"),1],pop_colors[which(pop_colors$pop == "CARL"),1],pop_colors[which(pop_colors$pop == "Erto"),1],pop_colors[which(pop_colors$pop == "Illegio"),1],pop_colors[which(pop_colors$pop == "Resia"),1],pop_colors[which(pop_colors$pop == "Sauris"),1])

  legend("bottomright",pch =c(rep(22,length(pops))),legend=leg_txt, pt.lwd=2,pt.cex=2,pt.bg=bkg,col=c(rep('black',length(pops))),ncol=4,bty="n")
dev.off()

######################################################################################################
# Plot 7: IBD
###############################################################################
#Same as roh, but we consider samples's pairs:
# input_format <- "PLINK"
chr <- "10"
pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
LOD <- 5

# base_folder <- getwd() #only if we're lazy
input_format <- "BEAGLE"
xmax <- NULL
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    # pop_ibd_file <- paste(pop,"roh.length.ibd",sep=".")
    # pop_ibd_file <- paste(pop,"roh.length.4.ibd",sep=".")
    pop_ibd_file <- paste("CHR",chr,"/",pop,".roh.length.",LOD,".ibd",sep="")
    pop_ibd_table <- read.table(pop_ibd_file,header=FALSE,sep=" ")
    colnames(pop_ibd_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","PAIR","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_ibd_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_ibd_table <- read.table(pop_ibd_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_ibd_table <- pop_ibd_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_ibd_table$IID2 <- as.character(pop_ibd_table$IID2)
  pop_ibd_table$IID1 <- as.character(pop_ibd_table$IID1)
  pop_ibd_table$PAIR <- as.character(pop_ibd_table$PAIR)
  current_pop_table_name <- paste(pop,"ibd",sep="_")
  assign(current_pop_table_name,pop_ibd_table)

  tot_ibd <- tapply(pop_ibd_table$LENGTH,pop_ibd_table$PAIR,sum)
  tot_ibd <- as.data.frame(tot_ibd, row.names=NULL)
  tot_ibd$ID <- rownames(tot_ibd)
  colnames(tot_ibd) <- c("IBD_tot","ID")
  tot_ibd$ID <- as.character(tot_ibd$ID)
  tot_ibd$IBD_tot <- as.numeric(as.character(tot_ibd$IBD_tot))
  tot_ibd$IBD_tot <- tot_ibd$IBD_tot/1000000

  assign(paste(pop,"tot_ibd",sep="_"),tot_ibd)
  
  assign(paste("M_ibd",pop,sep="_"),ecdf(tot_ibd$IBD_tot))
  xmax <- c(xmax,summary(ecdf(tot_ibd$IBD_tot))[6])

}

q_CEU_ibd <- quantile(CEU_tot_ibd$IBD_tot,c(.95))
q_TSI_ibd <- quantile(TSI_tot_ibd$IBD_tot,c(.95))
q_FVG_ibd <- quantile(FVG_tot_ibd$IBD_tot,c(.95))
q_VBI_ibd <- quantile(VBI_tot_ibd$IBD_tot,c(.95))
q_CARL_ibd <- quantile(CARL_tot_ibd$IBD_tot,c(.95))
q_Sauris_ibd <- quantile(Sauris_tot_ibd$IBD_tot,c(.95))
q_Erto_ibd <- quantile(Erto_tot_ibd$IBD_tot,c(.95))
q_Illegio_ibd <- quantile(Illegio_tot_ibd$IBD_tot,c(.95))
q_Resia_ibd <- quantile(Resia_tot_ibd$IBD_tot,c(.95))



#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.49898 
# TSI: 95% -> 1.242603 
# FVG: 95% -> 5.354979 
# VBI: 95% -> 2.017785 

#W50 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.473412
# TSI: 95% -> 1.32598
# FVG: 95% -> 5.255203
# VBI: 95% -> 1.994098

#W24 (Beagle RefinedIBD, filtered at LOD = 4)
# q_CEU_ibd : 95% -> 0.2779214 
# q_TSI_ibd : 95% -> 0.3357718 
# q_FVG_ibd : 95% -> 1.440876 
# q_VBI_ibd : 95% -> 0.6909859 
# q_CARL_ibd : 95% -> 0.5372137 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# q_CEU_ibd : 95% ->5.505728 
# q_TSI_ibd : 95% ->4.427433 
# q_FVG_ibd : 95% ->16.18415 
# q_VBI_ibd : 95% ->5.645495 
# q_CARL_ibd : 95% ->7.41046 
# q_Sauris_ibd : 95% ->35.33875 
# q_Erto_ibd : 95% ->25.70715 
# q_Illegio_ibd : 95% ->32.76607 
# q_Resia_ibd : 95% ->35.3291 


# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
# jpeg(paste(base_folder,"/8_ibd.jpg",sep=""),width=800, height=800)
#   par(lwd=3,cex=1.8)
#   plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
#   lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col="red", verticals=TRUE, pch=46)
#   lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col="green", verticals=TRUE, pch=46)
#   lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col="blue", verticals=TRUE, pch=46)
#   lines(M_ibd_CARL,CARL_tot_ibd$IBD_tot,col="yellow", verticals=TRUE, pch=46)
#   legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI","CARL"),col=c("black","red","green","blue","yellow"),ncol=5)
# dev.off()


jpeg(paste(base_folder,"/8_ibd_all_5POP_lod5_",chr,".jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=1.5)
  plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col="green", verticals=TRUE, pch=46)
  lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col="blue", verticals=TRUE, pch=46)
  lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col="red", verticals=TRUE, pch=46)
  lines(M_ibd_CARL,CARL_tot_ibd$IBD_tot,col="yellow", verticals=TRUE, pch=46)
  lines(M_ibd_Erto,Erto_tot_ibd$IBD_tot,col="darkred", verticals=TRUE, pch=46)
  lines(M_ibd_Illegio,Illegio_tot_ibd$IBD_tot,col="indianred", verticals=TRUE, pch=46)
  lines(M_ibd_Resia,Resia_tot_ibd$IBD_tot,col="mediumvioletred", verticals=TRUE, pch=46)
  lines(M_ibd_Sauris,Sauris_tot_ibd$IBD_tot,col="orangered", verticals=TRUE, pch=46)
  legend("bottomright",pch =c(rep(19,length(pops))),legend=pops,col=c("black","green","blue","red","yellow","darkred","indianred","mediumvioletred","orangered"),ncol=length(pops))
dev.off()

#plot the ROH cumulative dstribution (for length ?)

######################################################################################################
# Plot 7_A_IBD GENOMEWIDE for each POPULATION
###############################################################################
#Same as roh, but we consider samples's pairs:
# input_format <- "PLINK"
pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
LOD <- 5

# base_folder <- getwd() #only if we're lazy
input_format <- "BEAGLE"
xmax <- NULL
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    # pop_ibd_file <- paste(pop,"roh.length.ibd",sep=".")
    # pop_ibd_file <- paste(pop,"roh.length.4.ibd",sep=".")
    pop_ibd_file <- paste(pop,"/",pop,".WG.roh.length.",LOD,".ibd",sep="")
    pop_ibd_table <- read.table(pop_ibd_file,header=FALSE,sep=" ")
    colnames(pop_ibd_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","PAIR","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_ibd_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_ibd_table <- read.table(pop_ibd_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_ibd_table <- pop_ibd_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_ibd_table$IID2 <- as.character(pop_ibd_table$IID2)
  pop_ibd_table$IID1 <- as.character(pop_ibd_table$IID1)
  pop_ibd_table$PAIR <- as.character(pop_ibd_table$PAIR)
  current_pop_table_name <- paste(pop,"ibd",sep="_")
  assign(current_pop_table_name,pop_ibd_table)

  tot_ibd <- tapply(as.numeric(pop_ibd_table$LENGTH),pop_ibd_table$PAIR,sum)
  tot_ibd <- as.data.frame(tot_ibd, row.names=NULL)
  tot_ibd$ID <- rownames(tot_ibd)
  colnames(tot_ibd) <- c("IBD_tot","ID")
  tot_ibd$ID <- as.character(tot_ibd$ID)
  tot_ibd$IBD_tot <- as.numeric(as.character(tot_ibd$IBD_tot))
  tot_ibd$IBD_tot <- tot_ibd$IBD_tot/1000000

  assign(paste(pop,"tot_ibd",sep="_"),tot_ibd)
  
  assign(paste("M_ibd",pop,sep="_"),ecdf(tot_ibd$IBD_tot))
  xmax <- c(xmax,summary(ecdf(tot_ibd$IBD_tot))[6])

}

q_CEU_ibd <- quantile(CEU_tot_ibd$IBD_tot,c(.95))
q_TSI_ibd <- quantile(TSI_tot_ibd$IBD_tot,c(.95))
q_FVG_ibd <- quantile(FVG_tot_ibd$IBD_tot,c(.95))
q_VBI_ibd <- quantile(VBI_tot_ibd$IBD_tot,c(.95))
q_CARL_ibd <- quantile(CARL_tot_ibd$IBD_tot,c(.95))
q_Sauris_ibd <- quantile(Sauris_tot_ibd$IBD_tot,c(.95))
q_Erto_ibd <- quantile(Erto_tot_ibd$IBD_tot,c(.95))
q_Illegio_ibd <- quantile(Illegio_tot_ibd$IBD_tot,c(.95))
q_Resia_ibd <- quantile(Resia_tot_ibd$IBD_tot,c(.95))


summary_CEU_ibd <- summary(CEU_tot_ibd$IBD_tot)
summary_TSI_ibd <- summary(TSI_tot_ibd$IBD_tot)
summary_FVG_ibd <- summary(FVG_tot_ibd$IBD_tot)
summary_VBI_ibd <- summary(VBI_tot_ibd$IBD_tot)
summary_CARL_ibd <- summary(CARL_tot_ibd$IBD_tot)
summary_Sauris_ibd <- summary(Sauris_tot_ibd$IBD_tot)
summary_Erto_ibd <- summary(Erto_tot_ibd$IBD_tot)
summary_Illegio_ibd <- summary(Illegio_tot_ibd$IBD_tot)
summary_Resia_ibd <- summary(Resia_tot_ibd$IBD_tot)


#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.49898 
# TSI: 95% -> 1.242603 
# FVG: 95% -> 5.354979 
# VBI: 95% -> 2.017785 

#W50 (Beagle RefinedIBD, filtered at LOD = 5)
# CEU: 95% -> 1.473412
# TSI: 95% -> 1.32598
# FVG: 95% -> 5.255203
# VBI: 95% -> 1.994098

#W24 (Beagle RefinedIBD, filtered at LOD = 4)
# q_CEU_ibd : 95% -> 0.2779214 
# q_TSI_ibd : 95% -> 0.3357718 
# q_FVG_ibd : 95% -> 1.440876 
# q_VBI_ibd : 95% -> 0.6909859 
# q_CARL_ibd : 95% -> 0.5372137 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)
# q_CEU_ibd : 95% ->5.505728 
# q_TSI_ibd : 95% ->4.427433 
# q_FVG_ibd : 95% ->16.18415 
# q_VBI_ibd : 95% ->5.645495 
# q_CARL_ibd : 95% ->7.41046 
# q_Sauris_ibd : 95% ->35.33875 
# q_Erto_ibd : 95% ->25.70715 
# q_Illegio_ibd : 95% ->32.76607 
# q_Resia_ibd : 95% ->35.3291 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)GENOME WIDE 1/8/2014
# q_CEU_ibd: 95% -> 80.43721 
# q_TSI_ibd: 95% -> 61.59844 
# q_FVG_ibd: 95% -> 246.3632 
# q_VBI_ibd: 95% -> 61.37638 
# q_CARL_ibd: 95% -> 91.72759 
# q_Sauris_ibd: 95% -> 443.422 
# q_Erto_ibd: 95% -> 363.6143 
# q_Illegio_ibd: 95% -> 449.9667 
# q_Resia_ibd: 95% -> 412.6811 

#W24 (Beagle RefinedIBD, filtered at LOD = 5)GENOME WIDE 5/8/2014 ONLY 46 SAMPLES
# q_CEU_ibd: 95% -> 136.0912 
# q_TSI_ibd: 95% -> 122.4737 
# q_FVG_ibd: 95% -> 328.0547 
# q_VBI_ibd: 95% -> 162.723 
# q_CARL_ibd: 95% -> 157.1352 
# q_Sauris_ibd: 95% -> 461.4193 
# q_Erto_ibd: 95% -> 337.5515 
# q_Illegio_ibd: 95% -> 444.8172 
# q_Resia_ibd: 95% -> 429.6575 

# IBD           
# POP Min.  1st_Qu. Median  Mean  3rd_Qu. Max.
# summary_CEU_ibd 87.09 111.7 118.8 119.2 125.4 312.7
# summary_TSI_ibd 85.98 102.1 107.4 108.2 113.3 144.4
# summary_FVG_ibd 36.17 117.7 127.6 158.5 147.2 933.8
# summary_VBI_ibd 61.24 125.2 133.8 134.8 143.9 519.1
# summary_CARL_ibd  50.79 83.86 94.51 104 109.5 612.5
# summary_Sauris_ibd  84.33 212.1 258 275.3 315 1242
# summary_Erto_ibd  95.73 171.1 226.7 229.3 269 1151
# summary_Illegio_ibd 49.36 176.4 256.8 261 320 1101
# summary_Resia_ibd 89.86 247.4 292.5 301.8 336.6 2521


pop_colors <- col_pop(pops)
base_folder <- getwd()
jpeg(paste(base_folder,"/7_A_ibd_all_5POP_lod5_WG_1500.jpg",sep=""),width=1000, height=1000)
  par(lwd=4,cex=1.5)
  # plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple (Mb)", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "CEU"),1],main="",xlab="IBD segments length per couple (Mb)", xlim=c(0,1000), verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "TSI"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "VBI"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "FVG"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_CARL,CARL_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "CARL"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Erto,Erto_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Erto"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Illegio,Illegio_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Illegio"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Resia,Resia_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Resia"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  lines(M_ibd_Sauris,Sauris_tot_ibd$IBD_tot,col=pop_colors[which(pop_colors$pop == "Sauris"),1], verticals=TRUE, pch=46,col.01line='black',yaxs='i')
  abline(h=0.95,col='grey',lty='dashed')
  leg_txt <- c(pop_colors[which(pop_colors$pop == "CEU"),2],pop_colors[which(pop_colors$pop == "TSI"),2],pop_colors[which(pop_colors$pop == "VBI"),2],pop_colors[which(pop_colors$pop == "FVG"),2],pop_colors[which(pop_colors$pop == "CARL"),2],pop_colors[which(pop_colors$pop == "Erto"),2],pop_colors[which(pop_colors$pop == "Illegio"),2],pop_colors[which(pop_colors$pop == "Resia"),2],pop_colors[which(pop_colors$pop == "Sauris"),2])
  bkg <- c(pop_colors[which(pop_colors$pop == "CEU"),1],pop_colors[which(pop_colors$pop == "TSI"),1],pop_colors[which(pop_colors$pop == "VBI"),1],pop_colors[which(pop_colors$pop == "FVG"),1],pop_colors[which(pop_colors$pop == "CARL"),1],pop_colors[which(pop_colors$pop == "Erto"),1],pop_colors[which(pop_colors$pop == "Illegio"),1],pop_colors[which(pop_colors$pop == "Resia"),1],pop_colors[which(pop_colors$pop == "Sauris"),1])
  legend("bottomright",pch =c(rep(22,length(pops))),legend=leg_txt, pt.lwd=2,pt.cex=2,pt.bg=bkg,col=c(rep('black',length(pops))),ncol=4,bty="n")
dev.off()


# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/7_roh.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"PLOTS/7_roh_density.jpg",sep=""),width=800, height=800)
  plot(density(VBI_tot_roh$ROH_tot),main="",xlab="Total ROH homozigosity", col="blue")
  lines(density(FVG_tot_roh$ROH_tot),col="red")
  lines(density(TSI_tot_roh$ROH_tot),col="green")
  lines(density(CEU_tot_roh$ROH_tot))
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()


##########################################################################################################################
# PLOT 8: DAF spectrum for different functional annotations

rm(list=ls())
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
#outbred
o_pops <- c("TSI","CEU")
#isolates
i_pops <- c("FVG","VBI")
all_pops <- c(i_pops,o_pops)

chr <- "22"
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
# maf_classes <- c(0,0.01,0.02,0.05,0.10,0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1)
maf_classes <- c(0,0.05,0.10,0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1)
conseq_classes <- paste(base_folder,"INPUT_FILES/consequences.list",sep="")
conseq <- read.table(conseq_classes,header=F)
conseq$V1 <- as.character(conseq$V1)
# conseq <- c("missense_variant","synonymous_variant")
# cons <- conseq$V1[1]
# pop <- "FVG"

# for (cons in conseq) {
all_pop_maf_classes <- NULL
all_pop_maf_classes_fract <- NULL

for (pop in all_pops){
  pop_all_maf_classes <- NULL
  pop_all_maf_classes_fract <- NULL
  pop_maf_name <- paste(pop,"_all_cons_resume",sep="")
  pop_maf_name_fract <- paste(pop,"_all_cons_resume_fract",sep="")

  for (cons in conseq$V1) {
    # file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_no_fixed_MAC2/",pop,"/",sep="")
    file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_private_no_fixed/",pop,"/",sep="")
    # file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_no_fixed/",pop,"/",sep="")
    # file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_fixed/",pop,"/",sep="")
    # file_name <- paste(pop,".",cons,".",chr,".tab.gz",sep="")
    file_name <- paste(pop,".private.",cons,".",chr,".tab.gz",sep="")
    pop_table <- read.table(paste(file_path,file_name,sep=""),header=T,stringsAsFactors=F, comment.char="")
    #replace header because we know already how it is and we need to have a column called "CHROM", it's MANDATORY!
    colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")

    if(length(pop_table$DAF) > 0) {

      pop_table$DAC <- as.numeric(as.character(pop_table$DAC))
      pop_table$MAC <- as.numeric(as.character(pop_table$MAC))
      pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
      pop_table$MAF <- as.numeric(as.character(pop_table$MAF))
      
      #remove all MAC = 1
      # pop_table <- pop_table[which(pop_table$MAC > 1),]
      
      #write a cute output
      outdir <- paste(pop,"_",chr,sep="")
      system(paste("mkdir -p ",pop,"_",chr,"/summaries",sep=""))

      #remove sites without the DAF info
      dim(pop_table)
      pop_table_no_na <- pop_table[which(!is.na(pop_table$DAF)),]
      dim(pop_table_no_na)
      summary(pop_table_no_na)
      pop_table_no_mono <- pop_table_no_na[which(pop_table_no_na$DAF !=1 | pop_table_no_na$DAF != 0),]
      all_maf_classes <- split_bins(maf_classes,pop_table_no_mono,file_name,"DAF",outdir)
      gc()
      
      sink(paste(outdir,"/summaries/",file_name,'_maf_bin_resume.txt',sep=""))
      print(all_maf_classes)
      sink()
      tot <- sum(all_maf_classes)
      current_pop_conseq <- cbind(all_maf_classes,tot,cons,pop)
      current_pop_conseq_fract <- cbind(all_maf_classes/tot,tot,cons,pop)
      #concat all resumes for plotting
      pop_all_maf_classes <- rbind(pop_all_maf_classes,current_pop_conseq)
      pop_all_maf_classes_fract <- rbind(pop_all_maf_classes_fract,current_pop_conseq_fract)
    }
  assign(pop_maf_name,pop_all_maf_classes)
  assign(pop_maf_name_fract,pop_all_maf_classes_fract)
  }
  all_pop_maf_classes <- rbind(all_pop_maf_classes,pop_all_maf_classes)
  all_pop_maf_classes_fract <- rbind(all_pop_maf_classes_fract,pop_all_maf_classes_fract)
}
#write the table: all data
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#private data
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#Data filtered by MAC > 1 
#write the table: all data
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_MAC2_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_MAC2_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#private data
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)


#now chose a couple of consequences and test plot
cons_to_plot <- c("missense_variant","synonymous_variant")

maf_to_plot <- c(0.05,0.10,0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1)
#read data for each population
#plot
all_cols <- NULL
# cats <- colnames(merged_daf_diag)

for(i in 1:length(all_pops)){
  if(all_pops[i] == "CEU"){
    cur_col <- colors()[130]
  }
  if(all_pops[i] == "FVG"){
    cur_col <- colors()[517]
  }
  if(all_pops[i] == "FVG_p"){
    cur_col <- colors()[50]
  }
  if(all_pops[i] == "FVG_s"){
    cur_col <- colors()[81]
  }
  if(all_pops[i] == "TSI"){
    cur_col <- colors()[30]
  }
  if(all_pops[i] == "VBI"){
    cur_col <- colors()[421]
  }
  if(all_pops[i] == "VBI_p"){
    cur_col <- colors()[33]
  }
  if(all_pops[i] == "VBI_s"){
    cur_col <- colors()[36]
  }
  
  pop_col <- cbind(cur_col,all_pops[i])
  all_cols <- rbind(all_cols,pop_col)
}
all_cols <- as.data.frame(all_cols)

jpeg(paste(base_folder,"PLOTS/10_conseq_daf.jpg",sep=""),width=800, height=800)
  par(lwd=2)
  for (cons in cons_to_plot) {
    # current_cat <- all_pop_maf_classes[which(all_pop_maf_classes$cons == "missense_variant"),]
    current_cat <- all_pop_maf_classes[which(all_pop_maf_classes$cons == cons),]

    plot(t(current_cat[which(current_cat$pop == "CEU"),1:11]),type="d",col=all_cols[which(all_cols$V2 == "CEU"),1])
    lines(t(current_cat[which(current_cat$pop == "FVG"),1:11]),type="l",col=all_cols[which(all_cols$V2 == "FVG"),1])
    lines(t(current_cat[which(current_cat$pop == "TSI"),1:11]),type="d",col=all_cols[which(all_cols$V2 == "TSI"),1])
    lines(t(current_cat[which(current_cat$pop == "VBI"),1:11]),type="l",col=all_cols[which(all_cols$V2 == "VBI"),1])
  }
  legend("topright",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=4)
dev.off()


barplot(barplot7,
      beside=T,
      main="DAF in all populations",
      names.arg=all_pop_DAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[c(1,2,3,5,6,7),1],ylim=c(0,60))
      mtext(1, text = "DAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Relative Frequency (N sites/Tot sites per category)(%)", line = 4,cex=1.4)
      legend("top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[c(1,2,3,5,6,7),1],legend=all_cols[c(1,2,3,5,6,7),2], ncol=6)


################################################################################################
# PLOT 8a: MAF spectrum for different functional annotations

rm(list=ls())
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
#outbred
o_pops <- c("TSI","CEU")
#isolates
i_pops <- c("FVG","VBI")
all_pops <- c(i_pops,o_pops)

chr <- "22"
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
# maf_classes <- c(0,0.01,0.02,0.05,0.10,0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1)
maf_classes <- c(0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.4,0.45,0.5)
conseq_classes <- paste(base_folder,"INPUT_FILES/consequences.list",sep="")
conseq <- read.table(conseq_classes,header=F)
conseq$V1 <- as.character(conseq$V1)
# conseq <- c("missense_variant","synonymous_variant")
# cons <- conseq$V1[1]
# pop <- "FVG"

# for (cons in conseq) {
all_pop_maf_classes <- NULL
all_pop_maf_classes_fract <- NULL

# for (pop in i_pops){
for (pop in all_pops){
  pop_all_maf_classes <- NULL
  pop_all_maf_classes_fract <- NULL
  pop_maf_name <- paste(pop,"_all_cons_resume",sep="")
  pop_maf_name_fract <- paste(pop,"_all_cons_resume_fract",sep="")

  for (cons in conseq$V1) {
    # file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_no_fixed_MAC2/",pop,"/",sep="")
    # file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_private_no_fixed/",pop,"/",sep="")
    file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_no_fixed/",pop,"/",sep="")
    # file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_fixed/",pop,"/",sep="")
    file_name <- paste(pop,".",cons,".",chr,".tab.gz",sep="")
    # file_name <- paste(pop,".private.",cons,".",chr,".tab.gz",sep="")
    pop_table <- read.table(paste(file_path,file_name,sep=""),header=T,stringsAsFactors=F, comment.char="")
    #replace header because we know already how it is and we need to have a column called "CHROM", it's MANDATORY!
    colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")

    if(length(pop_table$DAF) > 0) {

      pop_table$DAC <- as.numeric(as.character(pop_table$DAC))
      pop_table$MAC <- as.numeric(as.character(pop_table$MAC))
      pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
      pop_table$MAF <- as.numeric(as.character(pop_table$MAF))
      
      #remove all MAC = 1
      pop_table <- pop_table[which(pop_table$MAC > 1),]
      
      #write a cute output
      # outdir <- paste(pop,"_",chr,sep="")
      # outdir <- paste(pop,"_",chr,"_MAF",sep="")
      outdir <- paste(pop,"_",chr,"_MAF_MAC2",sep="")
      system(paste("mkdir -p ",outdir,"/summaries",sep=""))

      #remove sites without the DAF info
      dim(pop_table)
      pop_table_no_na <- pop_table[which(!is.na(pop_table$DAF)),]
      dim(pop_table_no_na)
      summary(pop_table_no_na)
      pop_table_no_mono <- pop_table_no_na[which(pop_table_no_na$DAF !=1 | pop_table_no_na$DAF != 0),]
      # for DAF
      # all_maf_classes <- split_bins(maf_classes,pop_table_no_mono,file_name,"DAF",outdir)
      # for MAF
      all_maf_classes <- split_bins(maf_classes,pop_table_no_mono,file_name,"MAF",outdir)
      gc()
      
      # for DAF
      # sink(paste(outdir,"/summaries/",file_name,'_daf_bin_resume.txt',sep=""))
      # for MAF
      sink(paste(outdir,"/summaries/",file_name,'_maf_bin_resume.txt',sep=""))
      print(all_maf_classes)
      sink()
      tot <- sum(all_maf_classes)
      current_pop_conseq <- cbind(all_maf_classes,tot,cons,pop)
      current_pop_conseq_fract <- cbind(all_maf_classes/tot,tot,cons,pop)
      #concat all resumes for plotting
      pop_all_maf_classes <- rbind(pop_all_maf_classes,current_pop_conseq)
      pop_all_maf_classes_fract <- rbind(pop_all_maf_classes_fract,current_pop_conseq_fract)
    }
  assign(pop_maf_name,pop_all_maf_classes)
  assign(pop_maf_name_fract,pop_all_maf_classes_fract)
  }
  all_pop_maf_classes <- rbind(all_pop_maf_classes,pop_all_maf_classes)
  all_pop_maf_classes_fract <- rbind(all_pop_maf_classes_fract,pop_all_maf_classes_fract)
}
###################################################
# DAF
#write the table: all data
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#private data
# DAF
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#Data filtered by MAC > 1 
#write the table: all data
# DAF
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_MAC2_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_MAC2_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#private data
# DAF
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_MAC2_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_daf_classes_private_MAC2_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)


###########################################################
# MAF
#write the table: all data
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#private data
# MAF
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_private_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_private_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#Data filtered by MAC > 1 
#write the table: all data
# MAF
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_MAC2_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_MAC2_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
#private data
# MAF
write.table(all_pop_maf_classes,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_private_MAC2_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(all_pop_maf_classes_fract,file=paste(base_folder,"RESULTS/CONSEQUENCES/all_pop_conseq_maf_classes_private_MAC2_fract_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

##########################################
#now chose a couple of consequences and test plot
cons_to_plot <- c("missense_variant","synonymous_variant")

maf_to_plot <- c(0.05,0.10,0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1)
#read data for each population
#plot
all_cols <- NULL
# cats <- colnames(merged_daf_diag)

for(i in 1:length(all_pops)){
  if(all_pops[i] == "CEU"){
    cur_col <- colors()[130]
  }
  if(all_pops[i] == "FVG"){
    cur_col <- colors()[517]
  }
  if(all_pops[i] == "FVG_p"){
    cur_col <- colors()[50]
  }
  if(all_pops[i] == "FVG_s"){
    cur_col <- colors()[81]
  }
  if(all_pops[i] == "TSI"){
    cur_col <- colors()[30]
  }
  if(all_pops[i] == "VBI"){
    cur_col <- colors()[421]
  }
  if(all_pops[i] == "VBI_p"){
    cur_col <- colors()[33]
  }
  if(all_pops[i] == "VBI_s"){
    cur_col <- colors()[36]
  }
  
  pop_col <- cbind(cur_col,all_pops[i])
  all_cols <- rbind(all_cols,pop_col)
}
all_cols <- as.data.frame(all_cols)

jpeg(paste(base_folder,"PLOTS/10_conseq_daf.jpg",sep=""),width=800, height=800)
  par(lwd=2)
  for (cons in cons_to_plot) {
    # current_cat <- all_pop_maf_classes[which(all_pop_maf_classes$cons == "missense_variant"),]
    current_cat <- all_pop_maf_classes[which(all_pop_maf_classes$cons == cons),]

    plot(t(current_cat[which(current_cat$pop == "CEU"),1:11]),type="d",col=all_cols[which(all_cols$V2 == "CEU"),1])
    lines(t(current_cat[which(current_cat$pop == "FVG"),1:11]),type="l",col=all_cols[which(all_cols$V2 == "FVG"),1])
    lines(t(current_cat[which(current_cat$pop == "TSI"),1:11]),type="d",col=all_cols[which(all_cols$V2 == "TSI"),1])
    lines(t(current_cat[which(current_cat$pop == "VBI"),1:11]),type="l",col=all_cols[which(all_cols$V2 == "VBI"),1])
  }
  legend("topright",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=4)
dev.off()



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

#we need to extract all cases of sharing in order to do the venn diagram
#CEU
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'
# #TSI
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'
# #VBI
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'
# #FVG
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'