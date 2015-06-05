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



########################################################################
#PLOT 1A: same as plot 1 but based on DAF: for the DAF categories, we plot the maf spectrum
#upload data for each population:
# chr <- "22"
# chr <- "12"
# pop <- "Resia"
# conseq <- conseqs[1]
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
conseqs <- c("condel.deleterious","condel.neutral","miss","neut","polyphen.benign","polyphen.possibly.damaging","polyphen.probably.damaging","sift.deleterious","sift.tolerated","syn")
# in_folder <- paste(getwd(),conseq,sep="/")
in_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/INPUT_FILES/FIVE_POPS/WG"

#the order of all pops is due to the merge order
# pops <- c("Erto","Illegio","Resia","Sauris","CEU","TSI","VBI","FVG","CARL")
pops <- c("Erto","Illegio","Resia","Sauris","CEU","TSI","VBI","CARL")
#wrapped in this for-loop, I'll write a table for each chr, so we can easily import the thing to plot, after...
all_pop_MAF <- NULL
for (pop in pops) {
  print(pop)
  all_chr_MAF <- NULL
  for (chr in 1:22) {
  # for (chr in 20:22) {
  print(chr)
  base_folder <- paste(in_folder,"/CHR",chr,sep="")

    # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
    pop_table_file <- paste(base_folder,"/",pop,".chr",chr,".tab.gz",sep="")
    # for files with not mac=1
    # pop_table_file <- paste(pop,".chr",chr,".not_fixed.not_MAC1.tab.gz",sep="")

    pop_table <- read.table(pop_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="",colClasses=c("integer",rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
    # pop_table <- read.table(pop_table_file,header=F,skip=1,nrows=100000,stringsAsFactors=F, comment.char="",colClasses=c(rep("integer",3),rep("character",3),"NULL",rep("integer",4),rep("numeric",2)))
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
    
    all_chr_MAF <- rbind(all_chr_MAF,pop_table)

  }
  all_pop_MAF <- append(all_pop_MAF,list(all_chr_MAF$DAF))
}

#save the R object so we just need to reload this, eventually
save(all_pop_MAF,file="all_pop_DAF.RData")


pop_col <- col_pop(pops)
require(plotrix)
all_pop_hist <- multhist(all_pop_MAF,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_MAF_table <- as.data.frame(cbind((all_pop_hist[[1]]$mids),t(all_pop_hist[[2]])))
colnames(all_pop_MAF_table) <- c("breaks",pops)

write.table(all_pop_MAF_table,file=paste("all_pop_DAF_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)


########################################################################################################
#PLOT 1B: same as plot 1 but based on DAF: for the consequences categories, we plot the daf spectrum
#upload data for each population:
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")

# pops <- c("CEU","TSI","VBI","FVG","CARL")
# base_conseq_maf_folder <- '/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/MAF'
base_conseq_maf_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/46_SAMPLES/RESULTS/DAF/20150519"

# conseq <- system('for i in `ls -d /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/listsites/*/`;do echo \"${i%*/}\"| awk \'BEGIN{FS=\"/\"};{print $(NF)}\';done',intern=T)
# conseq <- c("condel","sift","polyphen")
conseq <- c("condel.deleterious","condel.neutral","polyphen.benign","polyphen.possibly.damaging","polyphen.probably.damaging","sift.deleterious","sift.tolerated","syn","neut","miss")
# conseq <- c("syn","miss","condel","sift","polyphen","neut")
pops <- c("VBI","CARL","Erto","Illegio","Resia","Sauris")

#we need to select different categories
# categories <- c("shared","private","novel")
categories <- c("shared","private")
for (cat in categories){
  if (cat=="shared"){
    pops_ingi_class <- c("VBI_s","CARL_s","Erto_s","Illegio_s","Resia_s","Sauris_s")
  } else if (cat=="private"){
    pops_ingi_class <- c("VBI_p","CARL_p","Erto_p","Illegio_p","Resia_p","Sauris_p")
  } else if (cat=="novel"){
    pops_ingi_class <- c("VBI_n","CARL_n","Erto_n","Illegio_n","Resia_n","Sauris_n")
  }

  for (con in conseq){
    current_con_path <- paste(base_conseq_maf_folder,con,sep="/")
  #read the current file for this population and this chromosome
  # current_chr_novel <- read.table(paste(base_novel_folder,"/",pop,".",chr,".n.maf",sep=""), sep="\t",header=F)
    #we assume all files are splitted by chr
    if (con == "condel"){
      subcats <- c("deleterious","neutral")
    } else if (con == "polyphen"){
      subcats <- c("benign","possibly_damaging","probably_damaging")
    } else if (con == "sift"){
      subcats <- c("deleterious","tolerated")
    } else if (con == "gerp"){
      # those don't have subcategories
      subcats <- NULL
    } else if (con == "neut"){
     # those don't have subcategories
      subcats <- NULL
    } else if (con == "miss"){
      # those don't have subcategories
      subcats <- NULL
    } else if (con == "syn"){
      # those don't have subcategories
      subcats <- NULL
    } else {subcats <- NULL}

    if (!is.null(subcats)){

      for(subcat in subcats){
        all_pop_current_consec_current_subcat <- NULL
        for (pop in pops) {
          current_pop_current_consec_current_subcat <- NULL
          for (chr in 1:22){
            current_chr <- try(read.table(paste(current_con_path,"/CHR",chr,"/",pop,"_",cat,"_chr",chr,".merged_frq.tab.gz",sep=""), sep="\t",header=F))
            if (!inherits(current_chr, 'try-error')){ 
              print(pop)
              print(chr)
              # colnames(current_chr_cat) <- c("CHROM","POZ","POS","TYPE","CEU","TSI","VBI","FVG","CAR")
              colnames(current_chr) <- c("CHROM","POS","VT","CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
              current_pop_col <- grep(pop,colnames(current_chr))
              current_chr <- current_chr[which(current_chr[,current_pop_col] != 0),]
              current_pop_current_consec_current_subcat <-rbind(current_pop_current_consec_current_subcat,current_chr[,c(1,2,3,current_pop_col)])
            }
          }
          
          all_pop_current_consec_current_subcat <- append(all_pop_current_consec_current_subcat,list(current_pop_current_consec_current_subcat[,4]))
        }
        save(all_pop_current_consec_current_subcat,file=paste("all_pop_",con,"_",subcat,"_",cat,".RData"))

        pop_col <- col_pop(pops_ingi_class)
        require(plotrix)

        # now plot everything together
        jpeg(paste(base_conseq_maf_folder,"/1B_all_pop_all_MAF_",subcat,"_",con,"_",cat,"_plotrix_texture.jpg",sep=""),width=1800, height=800)
          par(oma=c(3,3,3,3),cex=1.4)
          multhist(all_pop_current_consec_current_subcat,
           col=pop_col[pops_ingi_class],
           # density=all_cols$density*20,
           freq=FALSE,
           ylab="Proportion of sites",
           xlab="DAF",
           breaks=20,
           ylim=c(0, 40),
           main=paste("DAF in all populations for ",con,subcat,"in",cat,"sites"))
          legend("topright",pch =c(rep(22,length(pop_col))),pt.lwd=2,pt.cex=2,pt.bg=pop_col[pops_ingi_class],col=c(rep('black',length(pops_ingi_class))),legend=names(pop_col), ncol=2,bty="n")
          # legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(all_cols$color))),legend=all_cols[,2], ncol=2,bty="n",density=all_cols$density*20)
        dev.off()

      }
    }else{
      all_pop_current_consec_current_subcat <- NULL
      for (pop in pops) {
        print(con)
        print(pop)

        current_pop_current_consec_current_subcat <- NULL
        for (chr in 1:22){
          print(chr)
          # current_chr <- read.table(paste(current_con_path,"/",con,".",chr,".list.maf_file",sep=""),header=F)
          current_chr <- try(read.table(paste(current_con_path,"/CHR",chr,"/",pop,"_",cat,"_chr",chr,".merged_frq.tab.gz",sep=""), sep="\t",header=F))
          if (!inherits(current_chr, 'try-error')){ 
            colnames(current_chr) <- c("CHROM","POS","VT","CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
            current_pop_col <- grep(pop,colnames(current_chr))
            current_chr <- current_chr[which(current_chr[,current_pop_col] != 0),]
            current_pop_current_consec_current_subcat <-rbind(current_pop_current_consec_current_subcat,current_chr[,c(1,2,3,current_pop_col)])
          }
        }
        # current_pop_col <- grep(pop,colnames(current_pop_current_consec_current_subcat)) #this is always the fourth column
        all_pop_current_consec_current_subcat <- append(all_pop_current_consec_current_subcat,list(current_pop_current_consec_current_subcat[,4]))
      }
      save(all_pop_current_consec_current_subcat,file=paste0("all_pop_",con,"_",cat,".RData"))
      pop_col <- col_pop(pops_ingi_class)
      require(plotrix)

      # now plot everything together
      jpeg(paste(base_conseq_maf_folder,"/1B_all_pop_all_MAF_",con,"_",cat,"_plotrix_texture.jpg",sep=""),width=1800, height=800)
        par(oma=c(3,3,3,3),cex=1.4)
        multhist(all_pop_current_consec_current_subcat,
         col=pop_col[pops_ingi_class],
         # density=all_cols$density*20,
         freq=FALSE,
         ylab="Proportion of sites",
         xlab="DAF",
         breaks=20,
         ylim=c(0, 40),
         main=paste("DAF in all populations for",con,"in", cat,"sites",sep=" "))
        legend("topright",pch =c(rep(22,length(pop_col))),pt.lwd=2,pt.cex=2,pt.bg=pop_col[pops_ingi_class],col=c(rep('black',length(pops_ingi_class))),legend=names(pop_col), ncol=2,bty="n")
        # legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(all_cols$color))),legend=all_cols[,2], ncol=2,bty="n",density=all_cols$density*20)
      dev.off()
    }
  }
}
