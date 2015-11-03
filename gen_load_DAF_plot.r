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
  colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","DAF")
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
all_pop_DAF <- NULL
for (pop in pops) {
  print(pop)
  all_chr_DAF <- NULL
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
    colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","DAF")
    pop_table$DAF <- as.numeric(as.character(pop_table$DAF))

    #remove nas
    #we want to keep the NA's because those are the novel variants!
    # pop_table <- pop_table[which(!is.na(pop_table$DAF)),]
    
    #remove fixed variants (the ones with DAF == 0 or == 1) ---> those are also the ones with DAF 0, basically!
    pop_table <- pop_table[which(pop_table$DAF != 0),]
    pop_table <- pop_table[which(pop_table$DAF != 1),]

    #remove MAC = 1
    # pop_table <- pop_table[which(pop_table$MAC > 1),]

    # current_hist <- hist(pop_table$DAF,plot=F,breaks=20)

    # #use relative frequency sites per bin count/total variants number
    # current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pop,length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS)),chr_length=(length(pop_table$POS))))
    
    all_chr_DAF <- rbind(all_chr_DAF,pop_table)

  }
  all_pop_DAF <- append(all_pop_DAF,list(all_chr_DAF$DAF))
}

#save the R object so we just need to reload this, eventually
save(all_pop_DAF,file="all_pop_DAF.RData")


pop_col <- col_pop(pops)
require(plotrix)
all_pop_hist <- multhist(all_pop_DAF,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_DAF_table <- as.data.frame(cbind((all_pop_hist[[1]]$mids),t(all_pop_hist[[2]])))
colnames(all_pop_DAF_table) <- c("breaks",pops)

write.table(all_pop_DAF_table,file=paste("all_pop_DAF_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)


########################################################################################################
#PLOT 1B: same as plot 1 but based on DAF: for the consequences categories, we plot the daf spectrum
#upload data for each population:
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")

# pops <- c("CEU","TSI","VBI","FVG","CARL")
# base_conseq_maf_folder <- '/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/DAF'
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
        jpeg(paste(base_conseq_maf_folder,"/1B_all_pop_all_DAF_",subcat,"_",con,"_",cat,"_plotrix_texture.jpg",sep=""),width=1800, height=800)
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
      jpeg(paste(base_conseq_maf_folder,"/1B_all_pop_all_DAF_",con,"_",cat,"_plotrix_texture.jpg",sep=""),width=1800, height=800)
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


################################################################################################################
##### 11/09/2015 
##### Replot again with data from the joint call set
###### REPLOT with ggplot
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
require(plotrix)
require(ggplot2)
require(reshape2)
base_folder <- getwd()
#WE DONT't have CARL pop anymore
pops <- c("CEU","TSI","VBI","FVG-E","FVG-I","FVG-R","FVG-S")
pops_for_tables <- c("CEU","TSI","VBI","Erto","Illegio","Resia","Sauris")
pops_p <- c("CEU_p","TSI_p","VBI_p","FVG-E_p","FVG-I_p","FVG-R_p","FVG-S_p")
pops_s <- c("CEU_s","TSI_s","VBI_s","FVG-E_s","FVG-I_s","FVG-R_s","FVG-S_s")
pops_n <- c("CEU","TSI","VBI","FVG_E","FVG_I","FVG_R","FVG_S")
pops_n_ingi <- c("VBI","FVG_E","FVG_I","FVG_R","FVG_S")
pops_ingi_novel <- c("VBI_n","FVG_n")
pops_ingi_class <- c("VBI_p","FVG_p","VBI_s","FVG_s")


all_pops <- c(pops_ingi_class_n,pops_ingi_class_s)
# all_pops <- c(pops,pops_ingi_novel,pops_ingi_class)

all_pop_DAF <- NULL
for(chr in 1:22){
  print(chr)
  # current_chr_daf <- read.table(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/REVISION_201508/FIGURE1/df_all/",chr,".joint.daf.sp",sep=""),header=T)
  current_chr_daf <- read.table(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/max/20150809_REVISION/PLOTS/df_all/",chr,".joint.daf.sp",sep=""),header=T)
  # current_chr_daf$ps_CEU  <- NULL
  # current_chr_daf$ps_Erto <- NULL
  # current_chr_daf$ps_Illegio <- NULL
  # current_chr_daf$ps_Resia <- NULL
  # current_chr_daf$ps_Sauris <- NULL
  # current_chr_daf$ps_TSI <- NULL
  # current_chr_daf$ps_VBI <- NULL

  all_pop_DAF <- rbind(all_pop_DAF,current_chr_daf)
}
 
# save(all_pop_DAF,file="all_pop_DAF.RData")
load('all_pop_DAF.RData')
all_pop_DAF <-all_pop_DAF[-which(all_pop_DAF$CEU==1 & all_pop_DAF$TSI==1 & all_pop_DAF$VBI==1 & all_pop_DAF$Erto==1 & all_pop_DAF$Illegio==1 & all_pop_DAF$Resia==1 & all_pop_DAF$Sauris==1), ]

all_pop_DAF_private <- all_pop_DAF[which(all_pop_DAF$state == "private"),]
all_pop_DAF_shared <- all_pop_DAF[which(all_pop_DAF$state == "shared"),]

all_pop_DAF_2 <- all_pop_DAF
all_pop_DAF$CHR <- NULL
all_pop_DAF$POS <- NULL
all_pop_DAF$state <- NULL
all_pop_DAF$sum <- NULL

all_pop_DAF_private$CHR <- NULL
all_pop_DAF_private$POS <- NULL
all_pop_DAF_private$state <- NULL

all_pop_DAF_shared$CHR <- NULL
all_pop_DAF_shared$POS <- NULL
all_pop_DAF_shared$state <- NULL


save(all_pop_DAF_private,file="all_pop_DAF_private.RData")
save(all_pop_DAF_shared,file="all_pop_DAF_shared.RData")

load('all_pop_novel_DAF.RData')
load('all_pop_DAF_private_shared.RData')
#To do after uploading the RData
load("all_pop_DAF_novel.RData")
load("all_pop_DAF_shared.RData")
load("all_pop_DAF_count_rel.RData")

all_pop_DAF_novel
all_pop_DAF_shared
all_pop_DAF_count_rel
all_pop_DAF_count_rel_rf

all_pop_DAF_table_reshaped <- NULL
all_pop_DAF_table_reshaped_shared <- NULL
all_pop_DAF_table_reshaped_private <- NULL
for(pop in pops_for_tables){
sel_col <- c("CHR","POS",pop,"state")
pop_current <- (all_pop_DAF[,sel_col])
# pop_current_nomono <- (pop_current[which(pop_current[,pop] !=0 & pop_current[,pop] != 1 ),])
pop_current_nomono <- (pop_current[which(pop_current[,pop] !=0 ),])
pop_current_nomono$CHR <- NULL
pop_current_nomono$POS <- NULL

pop_current_nomono_shared <- pop_current_nomono[which(pop_current_nomono$state == "shared"),]
pop_current_nomono_private <- pop_current_nomono[which(pop_current_nomono$state == "private"),]

pop_current_nomono$state <- NULL
pop_current_nomono_shared$state <- NULL
pop_current_nomono_private$state <- NULL

# current_pop_hist <- multhist(pop_current_nomono,freq=TRUE,breaks=seq(0.0,0.5,by=0.02),plot=F)
current_pop_hist <- multhist(pop_current_nomono,freq=FALSE,breaks=seq(0.0,1,by=0.02),plot=F)
current_pop_hist_shared <- multhist(pop_current_nomono_shared,freq=FALSE,breaks=seq(0.0,1,by=0.02),plot=F)
current_pop_hist_private <- multhist(pop_current_nomono_private,freq=FALSE,breaks=seq(0.0,1,by=0.02),plot=F)

current_pop_DAF_table <- as.data.frame(cbind((current_pop_hist[[1]]$mids),t(current_pop_hist[[2]])*2))
current_pop_DAF_table_shared <- as.data.frame(cbind((current_pop_hist_shared[[1]]$mids),t(current_pop_hist_shared[[2]])*2))
current_pop_DAF_table_private <- as.data.frame(cbind((current_pop_hist_private[[1]]$mids),t(current_pop_hist_private[[2]])*2))

colnames(current_pop_DAF_table) <- c("breaks",pop)
colnames(current_pop_DAF_table_shared) <- c("breaks",pop)
colnames(current_pop_DAF_table_private) <- c("breaks",pop)
current_pop_DAF_table_reshaped <- melt(current_pop_DAF_table, id='breaks')
current_pop_DAF_table_reshaped_shared <- melt(current_pop_DAF_table_shared, id='breaks')
current_pop_DAF_table_reshaped_private <- melt(current_pop_DAF_table_private, id='breaks')
all_pop_DAF_table_reshaped <- rbind(all_pop_DAF_table_reshaped,current_pop_DAF_table_reshaped)
all_pop_DAF_table_reshaped_shared <- rbind(all_pop_DAF_table_reshaped_shared,current_pop_DAF_table_reshaped_shared)
all_pop_DAF_table_reshaped_private <- rbind(all_pop_DAF_table_reshaped_private,current_pop_DAF_table_reshaped_private)

}

all_pop_DAF_table_reshaped_shared$cat <- "Shared"
all_pop_DAF_table_reshaped_private$cat <- "Private"


# ks.test(CEU_novel_maf[which(CEU_novel_maf <= 0.03)],VBI_novel_maf[which(VBI_novel_maf <= 0.03)])
# ks.test(CEU_novel_maf[which(CEU_novel_maf <= 0.03)],CARL_novel_maf[which(CARL_novel_maf <= 0.03)])
# ks.test(CEU_novel_maf[which(CEU_novel_maf <= 0.03)],FVG_I_novel_maf[which(FVG_I_novel_maf <= 0.03)])
# ks.test(CEU_novel_maf[which(CEU_novel_maf <= 0.03)],FVG_E_novel_maf[which(FVG_E_novel_maf <= 0.03)])
# ks.test(CEU_novel_maf[which(CEU_novel_maf <= 0.03)],FVG_R_novel_maf[which(FVG_R_novel_maf <= 0.03)])
# ks.test(CEU_novel_maf[which(CEU_novel_maf <= 0.03)],FVG_S_novel_maf[which(FVG_S_novel_maf <= 0.03)])

# ks.test(TSI_novel_maf[which(TSI_novel_maf <= 0.03)],VBI_novel_maf[which(VBI_novel_maf <= 0.03)]) 
# ks.test(TSI_novel_maf[which(TSI_novel_maf <= 0.03)],CARL_novel_maf[which(CARL_novel_maf <= 0.03)]) 
# ks.test(TSI_novel_maf[which(TSI_novel_maf <= 0.03)],FVG_I_novel_maf[which(FVG_I_novel_maf <= 0.03)]) 
# ks.test(TSI_novel_maf[which(TSI_novel_maf <= 0.03)],FVG_E_novel_maf[which(FVG_E_novel_maf <= 0.03)]) 
# ks.test(TSI_novel_maf[which(TSI_novel_maf <= 0.03)],FVG_R_novel_maf[which(FVG_R_novel_maf <= 0.03)]) 
# ks.test(TSI_novel_maf[which(TSI_novel_maf <= 0.03)],FVG_S_novel_maf[which(FVG_S_novel_maf <= 0.03)]) 

#all
# all_pop_hist <- multhist(all_pop_DAF,freq=FALSE,breaks=seq(0.00,0.5,by=0.02),plot=F)
# all_pop_hist <- multhist(all_pop_DAF,freq=TRUE,breaks=seq(0.00,0.5,by=0.02),plot=F)
# all_pop_DAF_table <- as.data.frame(cbind((all_pop_hist[[1]]$mids),t(all_pop_hist[[2]])))
# colnames(all_pop_DAF_table) <- c("breaks",pops)

# all_pop_DAF_novel_table
# all_pop_DAF_shared_table
# #private
# # private_hist <- multhist(all_pop_DAF_private,freq=FALSE,breaks=seq(0.00,0.5,by=0.02),plot=F)
# private_hist <- multhist(all_pop_DAF_private,freq=TRUE,breaks=seq(0.00,0.5,by=0.02),plot=F)
# all_pop_private_DAF_table <- as.data.frame(cbind((private_hist[[1]]$mids),t(private_hist[[2]])))
# colnames(all_pop_private_DAF_table) <- c("breaks",pops_p)

# #shared
# # all_pop_DAF_shared_hist <- multhist(all_pop_DAF_shared,freq=FALSE,breaks=seq(0.00,0.5,by=0.02),plot=F)
# all_pop_DAF_shared_hist <- multhist(all_pop_DAF_shared,freq=TRUE,breaks=seq(0.00,0.5,by=0.02),plot=F)
# all_pop_DAF_shared_table <- as.data.frame(cbind((all_pop_DAF_shared_hist[[1]]$mids),t(all_pop_DAF_shared_hist[[2]])))
# colnames(all_pop_DAF_shared_table) <- c("breaks",pops_s)


# all_cols <-col_pop(all_pops)
pops_ingi <- c("CEU","TSI","VBI","FVG-E","FVG-I","FVG-R","FVG-S")
all_cols <-col_pop(pops_ingi)

# all_pop_DAF_table_reshaped <- melt(all_pop_DAF_table, id='breaks')

# all_pop_DAF_count_rel_reshaped <- melt(all_pop_DAF_count_rel, id='breaks')
# all_pop_DAF_count_rel_rf_reshaped <- melt(all_pop_DAF_count_rel_rf, id='breaks')

# all_pop_private_DAF_table_reshaped <-melt(all_pop_private_DAF_table,id='breaks')
# all_pop_DAF_shared_table_reshaped <-melt(all_pop_DAF_shared_table,id='breaks')
# all_pop_DAF_private_shared_table_reshaped <- rbind(all_pop_private_DAF_table_reshaped,all_pop_DAF_shared_table_reshaped)

# all_pop_DAF_private_shared_table_reshaped <- melt(all_pop_DAF_private_shared_table, id='breaks')

# all_pop_novel_DAF_table_reshaped <- melt(all_pop_novel_DAF_table, id='breaks')

#merge data together
# all_pop_all_DAF_table_reshaped <- rbind(all_pop_DAF_table_reshaped,all_pop_DAF_private_shared_table_reshaped,all_pop_novel_DAF_table_reshaped)
# all_pop_all_DAF_table_reshaped <- rbind(all_pop_DAF_private_shared_table_reshaped)
all_pop_all_DAF_table_reshaped <- rbind(all_pop_DAF_table_reshaped_shared,all_pop_DAF_table_reshaped_private)
all_pop_all_DAF_table_reshaped <- all_pop_DAF_table_reshaped
# all_pop_all_DAF_table_reshaped$variable <- factor(all_pop_all_DAF_table_reshaped$variable, levels = pops_ingi)
# all_pop_all_DAF_table_reshaped$variable <- as.character(all_pop_all_DAF_table_reshaped$variable)
# all_pop_all_DAF_table_reshaped$variable <- factor(all_pop_all_DAF_table_reshaped$variable, levels = all_pops)

# all_pop_DAF_count_rel_reshaped$variable <- factor(all_pop_DAF_count_rel_reshaped$variable, levels = pops_ingi)
# all_pop_DAF_count_rel_rf_reshaped$variable <- factor(all_pop_DAF_count_rel_rf_reshaped$variable, levels = pops_ingi)
# all_pop_DAF_novel_shared_table_reshaped$variable <- factor(all_pop_DAF_novel_shared_table_reshaped$variable, levels = pops_ingi)

# all_pop_all_DAF_table_reshaped$cat <- "total"
# # all_pop_all_DAF_table_reshaped[grep("_n",all_pop_all_DAF_table_reshaped$variable),]$cat <- "novel"
# all_pop_all_DAF_table_reshaped[grep("_s",all_pop_all_DAF_table_reshaped$variable),]$cat <- "Shared"
# all_pop_all_DAF_table_reshaped[grep("_p",all_pop_all_DAF_table_reshaped$variable),]$cat <- "Private"

# all_pop_all_DAF_table_reshaped$pop <- all_pop_all_DAF_table_reshaped$variable
# all_pop_all_DAF_table_reshaped$pop <- gsub("_s","",all_pop_all_DAF_table_reshaped$pop )
# all_pop_all_DAF_table_reshaped$pop <- gsub("_p","",all_pop_all_DAF_table_reshaped$pop )
# all_pop_all_DAF_table_reshaped$pop <- factor(all_pop_all_DAF_table_reshaped$pop, levels = pops_ingi)


# all_pop_DAF_count_rel_reshaped$cat <- "total"
# all_pop_DAF_count_rel_reshaped[grep("_n",all_pop_DAF_count_rel_reshaped$variable),]$cat <- "Novel"
# all_pop_DAF_count_rel_reshaped[grep("_s",all_pop_DAF_count_rel_reshaped$variable),]$cat <- "Shared"
# all_pop_DAF_count_rel_reshaped[grep("_n",all_pop_DAF_count_rel_reshaped$variable),]$cat <- "Private"


all_pop_all_DAF_table_reshaped$pop <- all_pop_all_DAF_table_reshaped$variable
all_pop_all_DAF_table_reshaped$pop <- factor(all_pop_all_DAF_table_reshaped$pop, levels = pops_ingi)
all_pop_all_DAF_table_reshaped[grep("Erto",all_pop_all_DAF_table_reshaped$variable),]$pop <- "FVG-E"
all_pop_all_DAF_table_reshaped[grep("Illegio",all_pop_all_DAF_table_reshaped$variable),]$pop <- "FVG-I"
all_pop_all_DAF_table_reshaped[grep("Resia",all_pop_all_DAF_table_reshaped$variable),]$pop <- "FVG-R"
all_pop_all_DAF_table_reshaped[grep("Sauris",all_pop_all_DAF_table_reshaped$variable),]$pop <- "FVG-S"

# all_pop_DAF_count_rel_rf_reshaped$cat <- "total"
# all_pop_DAF_count_rel_rf_reshaped[grep("_n",all_pop_DAF_count_rel_rf_reshaped$variable),]$cat <- "Novel"
# all_pop_DAF_count_rel_rf_reshaped[grep("_s",all_pop_DAF_count_rel_rf_reshaped$variable),]$cat <- "Shared"
# all_pop_DAF_count_rel_rf_reshaped[grep("_n",all_pop_DAF_count_rel_rf_reshaped$variable),]$cat <- "Private"

# all_pop_DAF_count_rel_rf_reshaped$pop <- "all"
# all_pop_DAF_count_rel_rf_reshaped[grep("CEU",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "CEU"
# all_pop_DAF_count_rel_rf_reshaped[grep("TSI",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "TSI"
# all_pop_DAF_count_rel_rf_reshaped[grep("VBI",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "VBI"
# all_pop_DAF_count_rel_rf_reshaped[grep("CARL",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "CARL"
# all_pop_DAF_count_rel_rf_reshaped[grep("Erto",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "FVG-E"
# all_pop_DAF_count_rel_rf_reshaped[grep("Illegio",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "FVG-I"
# all_pop_DAF_count_rel_rf_reshaped[grep("Resia",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "FVG-R"
# all_pop_DAF_count_rel_rf_reshaped[grep("Sauris",all_pop_DAF_count_rel_rf_reshaped$variable),]$pop <- "FVG-S"
# all_pop_DAF_count_rel_rf_reshaped$pop <- factor(all_pop_DAF_count_rel_rf_reshaped$pop, levels = pops)

# all_pop_DAF_novel_shared_table_reshaped <- cbind(all_pop_DAF_count_rel_rf,CEU=(all_pop_DAF_count_rel$CEU/sum(all_pop_DAF_count_rel$CEU))*100,TSI=(all_pop_DAF_count_rel$TSI/sum(all_pop_DAF_count_rel$TSI))*100)



# all_pop_DAF_novel_shared_table_reshaped$cat <- "total"
# all_pop_DAF_novel_shared_table_reshaped[grep("_n",all_pop_DAF_novel_shared_table_reshaped$variable),]$cat <- "novel"
# all_pop_DAF_novel_shared_table_reshaped[grep("_s",all_pop_DAF_novel_shared_table_reshaped$variable),]$cat <- "Shared"
# all_pop_DAF_novel_shared_table_reshaped[grep("_n",all_pop_DAF_novel_shared_table_reshaped$variable),]$cat <- "Private"

# all_pop_DAF_novel_shared_table_reshaped$pop <- "all"
# all_pop_DAF_novel_shared_table_reshaped[grep("CEU",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "CEU"
# all_pop_DAF_novel_shared_table_reshaped[grep("TSI",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "TSI"
# all_pop_DAF_novel_shared_table_reshaped[grep("VBI",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "VBI"
# all_pop_DAF_novel_shared_table_reshaped[grep("CARL",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "CARL"
# all_pop_DAF_novel_shared_table_reshaped[grep("Erto",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "FVG-E"
# all_pop_DAF_novel_shared_table_reshaped[grep("Illegio",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "FVG-I"
# all_pop_DAF_novel_shared_table_reshaped[grep("Resia",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "FVG-R"
# all_pop_DAF_novel_shared_table_reshaped[grep("Sauris",all_pop_DAF_novel_shared_table_reshaped$variable),]$pop <- "FVG-S"
# all_pop_DAF_novel_shared_table_reshaped$pop <- factor(all_pop_DAF_novel_shared_table_reshaped$pop, levels = pops)

all_pop_all_DAF_table_reshaped_1 <- all_pop_all_DAF_table_reshaped[which(all_pop_all_DAF_table_reshaped$breaks <= 0.1) ,]
all_pop_all_DAF_table_reshaped_4 <- all_pop_all_DAF_table_reshaped[which(all_pop_all_DAF_table_reshaped$breaks >= 0.45) ,]
all_pop_all_DAF_table_reshaped_to_4 <- all_pop_all_DAF_table_reshaped[which(all_pop_all_DAF_table_reshaped$breaks <= 0.5) ,]
all_pop_all_DAF_table_reshaped_14 <- rbind(all_pop_all_DAF_table_reshaped_1,all_pop_all_DAF_table_reshaped_4) 
all_pop_all_DAF_table_reshaped_2 <- all_pop_all_DAF_table_reshaped[which(all_pop_all_DAF_table_reshaped$breaks <= 0.23) ,]

# all_pop_DAF_count_rel_reshaped_1 <- all_pop_DAF_count_rel_reshaped[which(all_pop_DAF_count_rel_reshaped$breaks <= 0.24) ,]
# all_pop_DAF_count_rel_reshaped_4 <- all_pop_DAF_count_rel_reshaped[which(all_pop_DAF_count_rel_reshaped$breaks >= 0.45) ,]
# all_pop_DAF_count_rel_reshaped_14 <- rbind(all_pop_DAF_count_rel_reshaped_1,all_pop_DAF_count_rel_reshaped_4) 
# all_pop_DAF_count_rel_reshaped_2 <- all_pop_DAF_count_rel_reshaped[which(all_pop_DAF_count_rel_reshaped$breaks <= 0.23) ,]


# all_pop_DAF_count_rel_rf_reshaped_1 <- all_pop_DAF_count_rel_rf_reshaped[which(all_pop_DAF_count_rel_rf_reshaped$breaks <= 0.24) ,]
# all_pop_DAF_count_rel_rf_reshaped_4 <- all_pop_DAF_count_rel_rf_reshaped[which(all_pop_DAF_count_rel_rf_reshaped$breaks >= 0.45) ,]
# all_pop_DAF_count_rel_rf_reshaped_14 <- rbind(all_pop_DAF_count_rel_rf_reshaped_1,all_pop_DAF_count_rel_rf_reshaped_4) 
# all_pop_DAF_count_rel_rf_reshaped_2 <- all_pop_DAF_count_rel_rf_reshaped[which(all_pop_DAF_count_rel_rf_reshaped$breaks <= 0.23) ,]

# all_pop_DAF_novel_shared_table_reshaped_1 <- all_pop_DAF_novel_shared_table_reshaped[which(all_pop_DAF_novel_shared_table_reshaped$breaks <= 0.24) ,]
# all_pop_DAF_novel_shared_table_reshaped_4 <- all_pop_DAF_novel_shared_table_reshaped[which(all_pop_DAF_novel_shared_table_reshaped$breaks >= 0.45) ,]
# all_pop_DAF_novel_shared_table_reshaped_14 <- rbind(all_pop_DAF_novel_shared_table_reshaped_1,all_pop_DAF_novel_shared_table_reshaped_4) 
# all_pop_DAF_novel_shared_table_reshaped_2 <- all_pop_DAF_novel_shared_table_reshaped[which(all_pop_DAF_novel_shared_table_reshaped$breaks <= 0.23) ,]


#plot 
# pl <- ggplot(all_pop_DAF_count_rel_reshaped_14)
# pl <- ggplot(all_pop_DAF_count_rel_reshaped_2)
# pl <- ggplot(all_pop_DAF_novel_shared_table_reshaped)
# pl <- ggplot(all_pop_DAF_novel_shared_table_reshaped_14)
# pl <- ggplot(all_pop_DAF_novel_shared_table_reshaped_2)

# pl <- pl + geom_bar(stat="identity",width=0.5,colour="black")
# pl <- pl + geom_bar(stat="bin",width=0.5, position = position_dodge(width=0.8),colour="black")
# pl <- pl + aes(x = factor(breaks), y = value, fill=variable)
# pl <- ggplot(all_pop_all_DAF_table_reshaped)
# pl <- ggplot(all_pop_all_DAF_table_reshaped_2)
# pl <- pl + geom_bar(stat="identity",width=0.5, position = position_dodge(width=0.8),colour="black")
# pl <- pl + aes(x = factor(breaks), y = value, fill=pop)
# pl <- pl + xlab("DAF")
# pl <- pl + ylab("Proportion of sites")
# # pl <- pl + ylab("Sites count")
# # pl <- pl + guides(shape = guide_legend(override.aes = list(colour = "pink")))
# # pl <- pl + scale_fill_manual("Cohorts", values=all_cols)
# # pl <- pl + scale_shape_manual(values=c(11,11))
# # pl <- pl + ggtitle(main)
# pl <- pl + scale_fill_manual("", values=all_cols)
# pl <- pl + facet_wrap( ~ cat, ncol=1)
# pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
# pl <- pl + theme_bw()

# pl <- pl + theme(strip.text.x = element_text(size = 20))
# pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
# pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
# pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
# pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
   
# # jpeg(paste(base_folder,"/",chr,"_point_dens.jpg",sep=""),width=1800, height=800)
# # ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# # ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_1.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# # ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_20150525.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# # ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_14_20150525.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# # ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_2_20150525.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# # ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_20150525_freq.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# # ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_14_20150525_freq.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# # png(paste(base_folder,"/1_all_pop_DAF_ggplot_2_20150525_freq.png",sep=""),width=1400, height=500,res=300)
# # png(paste(base_folder,"/1_all_pop_DAF_ggplot_14_20150525_freq.png",sep=""),width=1400, height=500,res=300)
# # png(paste(base_folder,"/1_all_pop_DAF_ggplot_20150525_freq.png",sep=""),width=1400, height=500,res=300)
# # png(paste(base_folder,"/1_all_pop_DAF_ggplot_20150525.png",sep=""),width=1400, height=500,res=300)
# # png(paste(base_folder,"/1_all_pop_DAF_ggplot_14_20150525.png",sep=""),width=1400, height=500,res=300)
# # png(paste(base_folder,"/1_all_pop_DAF_ggplot_2_20150525.png",sep=""),width=1400, height=500,res=300)
# png(paste(base_folder,"/1_all_pop_DAF_ggplot_2_20150911_freq.png",sep=""),width=1400, height=500,res=300)
# print(pl)
# dev.off()

# ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_2_20150525_freq.jpeg",sep=""),width=12, height=7,units="px",dpi=300,plot=pl)
# ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_2_20150911_freq.jpeg",sep=""),width=12, height=7,units="in",dpi=300,plot=pl)
# ggsave(filename=paste(base_folder,"/1_all_pop_DAF_ggplot_20150911_freq.jpeg",sep=""),width=12, height=7,units="in",dpi=300,plot=pl)


maf_sets_1 <- c("all_pop_DAF_count_rel_reshaped","all_pop_DAF_count_rel_reshaped_14","all_pop_DAF_count_rel_reshaped_2")
maf_sets_2 <- c("all_pop_DAF_count_rel_rf_reshaped","all_pop_DAF_count_rel_rf_reshaped_14","all_pop_DAF_count_rel_rf_reshaped_2")
maf_sets_3 <- c("all_pop_all_DAF_table_reshaped","all_pop_all_DAF_table_reshaped_2","all_pop_all_DAF_table_reshaped_4","all_pop_all_DAF_table_reshaped_to_4")
# for(set in maf_sets_2){
for(set in maf_sets_3){
  current_set <- get(set)
  pl <- ggplot(current_set)

  pl <- pl + geom_bar(stat="identity",width=0.5, position = position_dodge(width=0.8),colour="black")
  pl <- pl + aes(x = factor(breaks), y = value, fill=pop)
  pl <- pl + xlab("DAF")
  pl <- pl + ylab("Proportion of sites (%)")
  # pl <- pl + ylab("Site count")
  pl <- pl + scale_fill_manual("", values=all_cols)
  pl <- pl + facet_wrap( ~ cat, ncol=1)
  pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
  pl <- pl + theme_bw()

  pl <- pl + theme(strip.text.x = element_text(size = 20))
  pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
  pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
  pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
  pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
  
  # png(paste(base_folder,"/1_",set,"_all_fixed1_20150915.png",sep=""),width=1400, height=500)
  png(paste(base_folder,"/1_",set,"_split_fixed1_20150915.png",sep=""),width=1400, height=500)
  # png(paste(base_folder,"/1_",set,"_all_20150915.png",sep=""),width=1400, height=500)
  # png(paste(base_folder,"/1_",set,"_20150915.png",sep=""),width=1400, height=500)
  # png(paste(base_folder,"/1_",set,"_20150912_count.png",sep=""),width=1400, height=500)
  print(pl)
  dev.off()
  # ggsave(filename=paste(base_folder,"/1_",set,"_20150525.jpeg",sep=""),width=12, height=7,units="in",dpi=300,plot=pl)
}





