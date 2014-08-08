########################################################################
#plot 1 DAF spectrun for each population
#upload data for each population:
#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"


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

