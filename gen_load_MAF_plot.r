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

  # all_pop_MAF$breaks <- as.numeric(as.character(all_pop_MAF$breaks))
  # all_pop_MAF$counts <- as.numeric(as.character(all_pop_MAF$counts))
  # all_pop_MAF$rel_count <- as.numeric(as.character(all_pop_MAF$rel_count))
  # all_pop_MAF$chr_length <- as.numeric(as.character(all_pop_MAF$chr_length))
  # all_pop_MAF$pop <- as.character(all_pop_MAF$pop)

  # all_pop_MAF_table <- all_pop_MAF[which(all_pop_MAF$pop == "CEU"),1]
  # all_pop_MAF_table_names <- c("breaks")

  # for(pop_i in pops){
  #   actual_pop <- all_pop_MAF[which(all_pop_MAF$pop == pop_i),]
  #   all_pop_MAF_table_names <- c(all_pop_MAF_table_names,paste(pop_i,"count",sep="_"),paste(pop_i,"rel_count",sep="_"),paste(pop_i,"chr_length",sep="_"))

  #   all_pop_MAF_table <- cbind(all_pop_MAF_table,cbind(actual_pop$counts,actual_pop$rel_count*100,actual_pop$chr_length))
  # }
  # all_pop_MAF_table <- as.data.frame(all_pop_MAF_table)
  # colnames(all_pop_MAF_table) <- all_pop_MAF_table_names

# write.table(all_pop_MAF_table,file=paste("all_pop_MAF_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
  # write.table(all_pop_MAF_table,file=paste("all_pop_MAF_count_no_MAC1_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

# }
#In order to plot the genomewide spectrum, we need to upload each file and sum things together
# all_chr <- NULL
# for (chrom in 1:22) {
#   current_file <- paste("all_pop_MAF_count_chr",chrom,".txt",sep="")
#   current_chr <- read.table(current_file,sep="\t",header=T)
#   #remove the relative count
#   current_chr <- current_chr[,c(grep("rel_count",colnames(current_chr),invert=T))]
#   # colnames(current_chr)[which(colnames(current_chr) != "breaks")]
#   all_chr <- as.data.frame(c(all_chr,current_chr))
# }

#remove all unnecessary breaks 
# all_chr <- all_chr[,c(grep("breaks.",colnames(all_chr),invert=T))]

#now sum all together the count and the total length for each population
# for (pop in pops) {
#   all_chr$current_pop_tot_count <- apply(all_chr[,c(grep(paste(pop,"_count",sep=""),colnames(all_chr)))],1,sum)
#   all_chr$current_pop_tot_length <- apply(all_chr[,c(grep(paste(pop,"_chr_length",sep=""),colnames(all_chr)))],1,sum)
#   colnames(all_chr)[which(colnames(all_chr)=="current_pop_tot_count")] <- paste(pop,"_tot_count",sep="")
#   colnames(all_chr)[which(colnames(all_chr)=="current_pop_tot_length")] <- paste(pop,"_tot_length",sep="")
# }

#now we want to keep only the tot count stuff and the breaks column
# all_chr <- cbind(breaks=all_chr$breaks,all_chr[,c(grep("tot_count",colnames(all_chr)))],all_chr[,c(grep("tot_length",colnames(all_chr)))])

#now we need to create the table with the relative count
#we can do it strait away
#print the table, too!!
# write.table(all_chr,file=paste("all_pop_MAF_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
# all_pop_MAF_table <- cbind(breaks=all_chr$breaks,all_chr[,c(grep("tot_count",colnames(all_chr)))],(all_chr[,c(grep("tot_count",colnames(all_chr)))]/all_chr[,c(grep("tot_length",colnames(all_chr)))])*100)

#absolute count(even columns)
# barplot3 <- as.matrix(t(all_pop_MAF_table[,c(2,4,6,8,10)]))
# barplot3 <- as.matrix(t(all_pop_MAF_table[,seq(2, ncol(all_pop_MAF_table), by = 2)]))
# barplot3 <- as.matrix(t(all_chr[,c(grep("tot_count",colnames(all_chr)))]))

#relative site count (odd columns)
# barplot4 <- as.matrix(t(all_pop_MAF_table[,c(3,5,7,9,11)]))
# barplot4 <- as.matrix(t(all_pop_MAF_table[,seq(3, ncol(all_pop_MAF_table), by = 2)]))
# all_chr_rel <- (all_chr[,c(grep("tot_count",colnames(all_chr)))]/all_chr[,c(grep("tot_length",colnames(all_chr)))])*100)
# colnames(all_chr_rel) <- gsub("tot_count","rel_freq",colnames(all_chr_rel))

# barplot4 <- as.matrix(t(all_chr_rel))

# all_cols <- col_pop(pops)

# base_folder <-getwd()

# # jpeg(paste(base_folder,"/1_all_pop_MAF_",chr,".jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder,"/1_all_pop_MAF.jpg",sep=""),width=1800, height=800)
# # jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_no_MAC1.jpg",sep=""),width=1800, height=800)
#   par(oma=c(3,3,3,3),cex=1.4)
#   barplot(barplot3,
#     beside=T,
#     main="MAF in all populations",
#     # names.arg=all_pop_MAF_table$breaks,
#     names.arg=all_chr$breaks,
#     xlab="",
#     ylab="",col=all_cols[,1])
#     legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
#     mtext(1, text = "MAF", line = 4,cex=1.4)
#     mtext(2, text = "Site count", line = 4,cex=1.4)
# dev.off()

# # jpeg(paste(base_folder,"/1_all_pop_MAF_rf_",chr,".jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder,"/1_all_pop_MAF_rf.jpg",sep=""),width=1800, height=800)
# # jpeg(paste(base_folder,"PLOTS/1_all_pop_MAF_rf_no_MAC1.jpg",sep=""),width=1800, height=800)
#   par(oma=c(3,3,3,3),cex=1.4)
#   barplot(barplot4,
#     beside=T,
#     main="MAF in all populations",
#     # names.arg=all_pop_MAF_table$breaks,
#     names.arg=all_chr$breaks,
#     xlab="",
#     ylab="",col=all_cols[,1])
#     legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(pops))),legend=all_cols[,2], ncol=2,bty="n")
#     mtext(1, text = "MAF", line = 4,cex=1.4)
#     mtext(2, text = "Relative Frequency (N sites/Tot sites in freq bin)(%)", line = 4,cex=1.4)
# dev.off()
ACTUALLY
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
# colnames(merged_daf) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","FVG")

# merged_daf$CEU <- as.numeric(as.character(merged_daf$CEU))
# merged_daf$TSI <- as.numeric(as.character(merged_daf$TSI))
# merged_daf$VBI <- as.numeric(as.character(merged_daf$VBI))
# merged_daf$FVG <- as.numeric(as.character(merged_daf$FVG))
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
#we need to use the files with MAF info, but created after filtering by DAF: we need to put all chr together
pops_ingi_class <- c("VBI_p","FVG_p","CARL_p","VBI_s","FVG_s","CARL_s")
pops_ingi <- c("VBI","FVG","CAR")
pops_ingi_files <- c("pop_VBI_private","pop_FVG_private","pop_CARL_private","pop_VBI_shared","pop_FVG_shared","pop_CARL_shared")
var_class <- c("private","shared")
all_pop_private_MAF <- NULL
all_pop_shared_MAF <- NULL

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
      colnames(current_chr_current_type) <- c("CHROM","POZ","POS","TYPE","CEU","TSI","VBI","FVG","CAR")
      current_pop_col <- grep(pop,colnames(current_chr_current_type))
      current_chr_current_type <- current_chr_current_type[which(current_chr_current_type[,current_pop_col] != 0),]
      
      current_pop_current_type <-rbind(current_pop_current_type,current_chr_current_type[,c(1,2,3,4,current_pop_col)])
    }
    
    current_pop_col <- grep(pop,colnames(current_pop_current_type))
    all_pop_current_type <- append(all_pop_current_type,list(current_pop_current_type[,current_pop_col]))

  }
  all_pop_current_type_table_name <- paste("all_pop_MAF_",type,sep=0)
  assign(all_pop_current_type_table_name,all_pop_current_type)
  
}


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

#################################################################TEST TO PLOT##########################
   # base_folder <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/CHR",chr,sep="")
   # all_pop_MAF <- NULL

   #  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
   #  VBI_table_file <- paste(base_folder,"/","VBI",".chr",chr,".tab.gz",sep="")
   #  FVG_table_file <- paste(base_folder,"/","FVG",".chr",chr,".tab.gz",sep="")
   #  CARL_table_file <- paste(base_folder,"/","CARL",".chr",chr,".tab.gz",sep="")
   #  CEU_table_file <- paste(base_folder,"/","CARL",".chr",chr,".tab.gz",sep="")
   #  TSI_table_file <- paste(base_folder,"/","CARL",".chr",chr,".tab.gz",sep="")
   #  # for files with not mac=1
   #  # pop_table_file <- paste(pop,".chr",chr,".not_fixed.not_MAC1.tab.gz",sep="")

   #  VBI_pop_table <- read.table(VBI_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
   #  FVG_pop_table <- read.table(FVG_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
   #  CARL_pop_table <- read.table(CARL_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
   #  colnames(VBI_pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
   #  colnames(FVG_pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
   #  colnames(CARL_pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")
   #  VBI_pop_table$DAF <- as.numeric(as.character(VBI_pop_table$DAF))
   #  FVG_pop_table$DAF <- as.numeric(as.character(FVG_pop_table$DAF))
   #  CARL_pop_table$DAF <- as.numeric(as.character(CARL_pop_table$DAF))

   #  #remove nas
   #  #we want to keep the NA's because those are the novel variants!
   #  # pop_table <- pop_table[which(!is.na(pop_table$DAF)),]
    
   #  #remove fixed variants (the ones with DAF == 0 or == 1) ---> those are also the ones with MAF 0, basically!
   #  VBI_pop_table <- VBI_pop_table[which(VBI_pop_table$DAF != 0),]
   #  FVG_pop_table <- FVG_pop_table[which(FVG_pop_table$DAF != 0),]
   #  CARL_pop_table <- CARL_pop_table[which(CARL_pop_table$DAF != 0),]
   #  VBI_pop_table <- VBI_pop_table[which(VBI_pop_table$DAF != 1),]
   #  FVG_pop_table <- FVG_pop_table[which(FVG_pop_table$DAF != 1),]
   #  CARL_pop_table <- CARL_pop_table[which(CARL_pop_table$DAF != 1),]

   #  #remove MAC = 1
   #  all_pop_MAF <- list(VBI_pop_table$MAF,FVG_pop_table$MAF,CARL_pop_table$MAF)

##################################################################################################################
