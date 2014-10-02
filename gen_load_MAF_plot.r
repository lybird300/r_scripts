#r script to plot for enza
rm(list=ls())
# source("/home/max/Work/script/r_scripts/col_pop.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
in_folder <- getwd()
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
  base_folder <- paste(in_folder,"/CHR",chr,sep="")

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

#save the R object so we just need to reload this, eventually
save(all_pop_MAF,file="all_pop_MAF.RData")


pop_col <- col_pop(pops)

all_pop_hist <- multhist(all_pop_MAF,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_MAF_table <- as.data.frame(cbind((all_pop_hist[[1]]$mids),t(all_pop_hist[[2]])))
colnames(all_pop_MAF_table) <- c("breaks",pops)

write.table(all_pop_MAF_table,file=paste("all_pop_MAF_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
  
####################ONLY FOR PLOTTING POURPOSE!!!!#################
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
base_novel_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/NOVEL/novelmaf"
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

#save the R object so we just need to reload this, eventually
save(all_pop_novel_MAF,file="all_pop_novel_MAF.RData")

require(plotrix)
pop_col <- col_pop(pops_ingi_novel)
base_folder <- getwd()

novel_hist <- multhist(all_pop_novel_MAF,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_novel_MAF_table <- as.data.frame(cbind((novel_hist[[1]]$mids),t(novel_hist[[2]])))
colnames(all_pop_novel_MAF_table) <- c("breaks",pops_ingi_novel)

write.table(all_pop_novel_MAF_table,file=paste("all_ingi_pop_MAF_novel_count.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

####################ONLY FOR PLOTTING POURPOSE!!!!#################
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
      base_type_folder <- paste(in_folder,"/CHR",chr,sep="")
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

#save the R object so we just need to reload this, eventually
save(all_pop_MAF_private_shared,file="all_pop_MAF_private_shared.RData")

#now we should be able to proceed as before...and plot
pop_col <- col_pop(pops_ingi_class)
base_folder <- getwd()
require(plotrix)

all_pop_MAF_private_shared_hist <- multhist(all_pop_MAF_private_shared,
   freq=FALSE,
   breaks=20,
   plot=F)

all_pop_MAF_private_shared_table <- as.data.frame(cbind((all_pop_MAF_private_shared_hist[[1]]$mids),t(all_pop_MAF_private_shared_hist[[2]])))
colnames(all_pop_MAF_private_shared_table) <- c("breaks",pops_ingi_class)

write.table(all_pop_MAF_private_shared_table,file=paste("all_pop_MAF_private_shared_table.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)


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
  legend("topright",pch =c(rep(22,length(pop_col[,1]))),pt.lwd=2,pt.cex=2,pt.bg=pop_col[,1],col=c(rep('black',length(pops_ingi_class))),legend=pop_col[,2], ncol=2,bty="n")
dev.off()

#now we need to put all together to plot the complete spectum for all classes!
all_pop_all_MAF <- append(append(all_pop_MAF,all_pop_novel_MAF),all_pop_MAF_private_shared)

#we need pops in the same order as the list
all_pops <- c(pops,pops_ingi_novel,pops_ingi_class)

all_cols <- col_pop(all_pops)
require(plotrix)

all_pop_MAF_count_rel <- NULL
for(i in 1:length(all_pops)){
  current_pop <- all_pops[i]
  current_pop_MAF_hist <- hist(all_pop_all_MAF[[i]],
    breaks=20,
    plot=F)
  assign(paste(current_pop,"_maf_hist",sep=""),current_pop_MAF_hist)
  all_pop_MAF_count_rel <- cbind(all_pop_MAF_count_rel,current_pop_MAF_hist$counts) 
}

all_pop_MAF_count_rel <- as.data.frame(all_pop_MAF_count_rel)
#add colnames (we can use the same order of the population data)
colnames(all_pop_MAF_count_rel) <- all_pops
#add a column for breaks
all_pop_MAF_count_rel$breaks <- CEU_maf_hist$breaks[2:length(CEU_maf_hist$breaks)]

all_ingi_MAF_not_all <- NULL
#we need to add also the column with all - the sum of categories, for the ingi pops
for(inpop in pops_ingi){
  current_pop_col <- grep(inpop,colnames(all_pop_MAF_count_rel))
  current_pop_not_all <- all_pop_MAF_count_rel[current_pop_col[1]] - apply(all_pop_MAF_count_rel[current_pop_col[2:length(current_pop_col)]],1,sum)
  all_ingi_MAF_not_all <- c(all_ingi_MAF_not_all,current_pop_not_all)
}
all_ingi_MAF_not_all <- as.data.frame(all_ingi_MAF_not_all)
colnames(all_ingi_MAF_not_all) <- paste(colnames(all_ingi_MAF_not_all),"no_all",sep="_")

#add this set to the complete set
all_pop_MAF_count_rel <- cbind(all_pop_MAF_count_rel,all_ingi_MAF_not_all)

#now create a relatove freq set
all_pop_MAF_count_rel_rf <- NULL
#we need to add also the column with all - the sum of categories, for the ingi pops
for(inpop in pops_ingi){
  current_pop_col_rf <- grep(inpop,colnames(all_pop_MAF_count_rel))
  current_pop_not_all_rf <- ((all_pop_MAF_count_rel[current_pop_col_rf[2:length(current_pop_col_rf)]])/sum(all_pop_MAF_count_rel[current_pop_col_rf[1]]))*100
  all_pop_MAF_count_rel_rf <- c(all_pop_MAF_count_rel_rf,current_pop_not_all_rf)
}
all_pop_MAF_count_rel_rf <- as.data.frame(all_pop_MAF_count_rel_rf)

#add CEU and TSI to this
all_pop_MAF_count_rel_rf <- cbind(all_pop_MAF_count_rel_rf,CEU=(all_pop_MAF_count_rel$CEU/sum(all_pop_MAF_count_rel$CEU))*100,TSI=(all_pop_MAF_count_rel$TSI/sum(all_pop_MAF_count_rel$TSI))*100)

all_pop_MAF_count_rel_rf <- as.data.frame(all_pop_MAF_count_rel_rf)

barplot1 <- as.matrix(t(all_pop_MAF_count_rel_rf[,c(13,14,3,7,11,2,6,10,1,5,9)]))

all_cols <- col_pop(rownames(barplot1))
all_pop_no_INGI_MAF <- list(all_pop_all_MAF[[1]],all_pop_all_MAF[[2]],all_pop_all_MAF[[6]],all_pop_all_MAF[[7]],all_pop_all_MAF[[8]],all_pop_all_MAF[[9]],all_pop_all_MAF[[10]],all_pop_all_MAF[[11]],all_pop_all_MAF[[12]],all_pop_all_MAF[[13]],all_pop_all_MAF[[14]])

# all_cols <- col_pop(all_pops)


# now plot everything together
# jpeg(paste(base_folder,"/1_all_pop_all_MAF_plotrix.jpg",sep=""),width=1800, height=800)
# jpeg(paste(base_folder,"/1_all_pop_all_MAF_plotrix_texture.jpg",sep=""),width=1800, height=800)
jpeg(paste(base_folder,"/1_all_pop_all_MAF_plotrix_texture_no_INGI.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  multhist(all_pop_no_INGI_MAF,
  # multhist(all_pop_all_MAF,
  # multhist(all_pop_all_MAF_cat,
  # multhist(barplot1,
  # barplot(barplot1,
   col=all_cols$color,
   density=all_cols$density*20,
   freq=FALSE,
   names.arg=all_pop_MAF_count_rel$breaks,
   ylab="Relative Frequency (N sites/Tot sites in freq bin)(%)",
   xlab="MAF",
   breaks=20,
   ylim=c(0, 50),
   main="MAF in all populations",beside=T)
  legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(all_pops))),legend=all_cols[,2], ncol=2,bty="n",density=all_cols$density*20)
dev.off()

#################################################################################
#PLOT 1B: same as plot 1 but based on MAF: for the consequences categories, we plot the maf spectrum
#upload data for each population:
source("/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r")
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")

pops <- c("CEU","TSI","VBI","FVG","CARL")

base_conseq_maf_folder <- '/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/MAF'

# conseq <- system('for i in `ls -d /lustre/scratch113/projects/esgi-vbseq/20140430_purging/enza/listsites/*/`;do echo \"${i%*/}\"| awk \'BEGIN{FS=\"/\"};{print $(NF)}\';done',intern=T)
# conseq <- c("condel","sift","polyphen")
conseq <- c("syn","miss","gerp")

for (con in conseq){
  current_con_path <- paste(base_conseq_maf_folder,con,sep="/")
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
  } else if (con == "miss"){
    # those don't have subcategories
    subcats <- NULL
  } else if (con == "syn"){
    # those don't have subcategories
    subcats <- NULL
  }

  if (!is.null(subcats)){

    for(subcat in subcats){
      all_pop_current_consec_current_subcat <- NULL
      for (pop in pops) {
        current_pop_current_consec_current_subcat <- NULL
        for (chr in 1:22){
          current_chr <- read.table(paste(current_con_path,"/",con,".",subcat,".",chr,".list.maf_file",sep=""),header=F)
          colnames(current_chr) <- c("CHROM","POS","TYPE","CEU","TSI","VBI","FVG","CARL")
          current_pop_col <- grep(pop,colnames(current_chr))
          current_chr <- current_chr[which(current_chr[,current_pop_col] != 0),]
          current_pop_current_consec_current_subcat <-rbind(current_pop_current_consec_current_subcat,current_chr[,c(1,2,3,current_pop_col)])
        }
        
        current_pop_col <- grep(pop,colnames(current_pop_current_consec_current_subcat))
        all_pop_current_consec_current_subcat <- append(all_pop_current_consec_current_subcat,list(current_pop_current_consec_current_subcat[,current_pop_col]))
      }
      all_cols <- col_pop(pops)
      require(plotrix)

      # now plot everything together
      jpeg(paste(base_conseq_maf_folder,"/1B_all_pop_all_MAF_",subcat,"_",con,"_plotrix_texture.jpg",sep=""),width=1800, height=800)
        par(oma=c(3,3,3,3),cex=1.4)
        multhist(all_pop_current_consec_current_subcat,
         col=all_cols$color,
         density=all_cols$density*20,
         freq=FALSE,
         ylab="Relative Frequency (N sites/Tot sites in freq bin)(%)",
         xlab="MAF",
         breaks=20,
         ylim=c(0, 40),
         main=paste("MAF in all populations for ",con," ",subcat))
        legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(all_cols$color))),legend=all_cols[,2], ncol=2,bty="n",density=all_cols$density*20)
      dev.off()

    }
  }else{
    all_pop_current_consec_current_subcat <- NULL
    for (pop in pops) {
      current_pop_current_consec_current_subcat <- NULL
      for (chr in 1:22){
        current_chr <- read.table(paste(current_con_path,"/",con,".",chr,".list.maf_file",sep=""),header=F)
        colnames(current_chr) <- c("CHROM","POS","TYPE","CEU","TSI","VBI","FVG","CARL")
        current_pop_col <- grep(pop,colnames(current_chr))
        current_chr <- current_chr[which(current_chr[,current_pop_col] != 0),]
        current_pop_current_consec_current_subcat <-rbind(current_pop_current_consec_current_subcat,current_chr[,c(1,2,3,current_pop_col)])
      }
      
      current_pop_col <- grep(pop,colnames(current_pop_current_consec_current_subcat))
      all_pop_current_consec_current_subcat <- append(all_pop_current_consec_current_subcat,list(current_pop_current_consec_current_subcat[,current_pop_col]))
    }
    all_cols <- col_pop(pops)
    require(plotrix)

    # now plot everything together
    jpeg(paste(base_conseq_maf_folder,"/1B_all_pop_all_MAF_",con,"_plotrix_texture.jpg",sep=""),width=1800, height=800)
      par(oma=c(3,3,3,3),cex=1.4)
      multhist(all_pop_current_consec_current_subcat,
       col=all_cols$color,
       density=all_cols$density*20,
       freq=FALSE,
       ylab="Relative Frequency (N sites/Tot sites in freq bin)(%)",
       xlab="MAF",
       breaks=20,
       ylim=c(0, 40),
       main=paste("MAF in all populations for",con,sep=" "))
      legend("topright",pch =c(rep(22,length(all_cols[,1]))),pt.lwd=2,pt.cex=2,pt.bg=all_cols[,1],col=c(rep('black',length(all_cols$color))),legend=all_cols[,2], ncol=2,bty="n",density=all_cols$density*20)
    dev.off()
  }
}

###################################################################################
###### REPLOT with ggplot
require(ggplot2)
require(reshape2)
all_cols <-col_pop(all_pops)


ylab <- "Relative Frequency (N sites/Tot sites in freq bin)(%)"
xlab <- "MAF"
main <- "MAF in all populations"

all_pop_MAF_table_reshaped <- melt(all_pop_MAF_table, id='breaks')
all_pop_MAF_private_shared_table_reshaped <- melt(all_pop_MAF_private_shared_table, id='breaks')
all_pop_novel_MAF_table_reshaped <- melt(all_pop_novel_MAF_table, id='breaks')

#merge data together
all_pop_all_MAF_table_reshaped <- rbind(all_pop_MAF_table_reshaped,all_pop_MAF_private_shared_table_reshaped,all_pop_novel_MAF_table_reshaped)
all_pop_all_MAF_table_reshaped$variable <- as.character(all_pop_all_MAF_table_reshaped$variable)
all_pop_all_MAF_table_reshaped$variable <- factor(all_pop_all_MAF_table_reshaped$variable, levels = all_pops)

all_pop_all_MAF_table_reshaped$cat <- "total"
all_pop_all_MAF_table_reshaped[grep("_n",all_pop_all_MAF_table_reshaped$variable),]$cat <- "novel"
all_pop_all_MAF_table_reshaped[grep("_s",all_pop_all_MAF_table_reshaped$variable),]$cat <- "shared"
all_pop_all_MAF_table_reshaped[grep("_p",all_pop_all_MAF_table_reshaped$variable),]$cat <- "private"

all_pop_all_MAF_table_reshaped_1 <- all_pop_all_MAF_table_reshaped[which(all_pop_all_MAF_table_reshaped$breaks <= 0.1) ,]
all_pop_all_MAF_table_reshaped_4 <- all_pop_all_MAF_table_reshaped[which(all_pop_all_MAF_table_reshaped$breaks >= 0.45) ,]
all_pop_all_MAF_table_reshaped_14 <- rbind(all_pop_all_MAF_table_reshaped_1,all_pop_all_MAF_table_reshaped_4) 
all_pop_all_MAF_table_reshaped_2 <- all_pop_all_MAF_table_reshaped[which(all_pop_all_MAF_table_reshaped$breaks <= 0.23) ,]

# pl <- ggplot(all_pop_all_MAF_table_reshaped)
# pl <- ggplot(all_pop_all_MAF_table_reshaped_1)
pl <- ggplot(all_pop_all_MAF_table_reshaped_2)
# pl <- ggplot(all_pop_all_MAF_table_reshaped_14)
pl <- pl + geom_bar(stat="identity",position="dodge",colour="black")
# pl <- pl + geom_bar(stat="identity",position="dodge",colour="white")
pl <- pl + aes(x = factor(breaks), y = value, fill=variable)
pl <- pl + xlab(xlab)
pl <- pl + ylab(ylab)
pl <- pl + guides(fill=guide_legend(title="Cohorts"))
# pl <- pl + scale_x_continuous(breaks=c(all_pop_all_MAF_table_reshaped_1$breaks,all_pop_all_MAF_table_reshaped_4$breaks),limits=c(0,0.1,0.4,0.5))
# pl <- pl + scale_fill_manual("Cohorts", values=all_cols$color)
pl <- pl + ggtitle(main)
pl <- pl + theme_bw()
# pl <- pl + facet_wrap(~cats)
  
# jpeg(paste(base_folder,"/",chr,"_point_dens.jpg",sep=""),width=1800, height=800)
# ggsave(filename=paste(base_folder,"/1_all_pop_MAF_ggplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(base_folder,"/1_all_pop_MAF_ggplot_1.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
ggsave(filename=paste(base_folder,"/1_all_pop_MAF_ggplot_2.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(base_folder,"/1_all_pop_MAF_ggplot_14.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(base_folder,"/1_all_pop_MAF_ggplot_panels.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(base_folder,"/1_all_pop_MAF_ggplot_panels_1.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(base_folder,"/1_all_pop_MAF_ggplot_panels_2.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)



