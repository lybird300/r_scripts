#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"

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


###############################################################
# PLOT 8b: count load of mutation for each sample by category

rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")

#INGI
# pops <- sort(c("VBI","FVG","CARL"))
categs <- c("private","shared","novel")
#OUTBRED
pops <- c("CEU","TSI","CAR","VBI","FVG")
# categs <- c("all")
# csqs <- c("CHR22","miss","syn")

csqs <- c("CHR22","miss","syn","polyphen.benign","polyphen.possibly.damaging","polyphen.probably.damaging","sift.deleterious","sift.tolerated")


for(csq in csqs){
  for(cat in categs){
    print(cat)
    print(csq)
    base_folder <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/CONSEQUENCES/",csq,"/",cat,"/VCF",sep="")
    print(base_folder)
    setwd(base_folder)
    #upload files with format: CHR POS REF ALT AC AN [SAMPLES/GT]
    # we need to convert all genotypes to 1 0 2, we can create a function to apply to our data for genotype conversion...
    for (pop in pops) {
      current_region_current_current_pop <- NULL
      for(chr in 1:22) {
        if (cat == "all"){
          current_region_file <- paste(base_folder,"/",chr,".",pop,".chr",chr,".tab.gz.",csq,".regions.vcf.gz.tab",sep="")
        } else {
          # 1.INGI_chr1.merged_maf.tab.gz.VBI.shared.tab.gz.miss.regions.vcf.gz.tab
          if (pop == "CAR" && cat != "novel"){
            pop <- "CARL"
          }
          print(pop)
          current_region_file <- paste(base_folder,"/",chr,".INGI_chr",chr,".merged_maf.tab.gz.",pop,".",cat,".tab.gz.",csq,".regions.vcf.gz.tab",sep="")
        }
        print(current_region_file)
        if (file.exists(current_region_file)){
          current_region_current_chr_current_pop <- read.table(current_region_file,header=T)
          if (dim(current_region_current_chr_current_pop)[1] > 0 ){
            #we need to convert our genotypes, now!
            for(j in 7:length(colnames(current_region_current_chr_current_pop))){
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "0|0"] <- "0" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "0|1"] <- "1" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "1|0"] <- "1" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "1|1"] <- "2" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "2|0"] <- "1" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "0|2"] <- "1" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "2|1"] <- "1" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "1|2"] <- "1" 
                levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "2|2"] <- "2" 
                current_region_current_chr_current_pop[,j] <- as.numeric(as.character(current_region_current_chr_current_pop[,j]))
                
            }
            # EGAN00001172162
            current_region_current_chr_current_pop_samples <- current_region_current_chr_current_pop[,7:length(colnames(current_region_current_chr_current_pop))]
            current_region_current_current_pop <- rbind(current_region_current_current_pop,current_region_current_chr_current_pop_samples)
          }
        }
      }
      #current_region_current_current_pop contains ALL SNP variants for that category for each sample across all chromosome
      #so if I want to caculate the frequency instead of a simple count, I need to divide by this number of variants

      if (!is.null(current_region_current_current_pop)){
        #sum across al variants on all chromosomes
        current_pop_current_cat_current_type <- apply(current_region_current_current_pop,2,sum)
        #   FORMAT IN THIS WAY: sample num pop
        current_pop_current_cat_current_type_df <- as.data.frame(current_pop_current_cat_current_type)
        current_pop_current_cat_current_type_df$samples <- row.names(current_pop_current_cat_current_type_df)
        if (pop == "CARL"){
          current_pop_current_cat_current_type_df$pop <- "CAR"
        }else{
          current_pop_current_cat_current_type_df$pop <- pop
        }
        if(csq == "miss") {
          current_pop_current_cat_current_type_df$cons <- "Missense"
        }else if(csq == "syn"){
          current_pop_current_cat_current_type_df$cons <- "Synonymous"
            
        }else {
          current_pop_current_cat_current_type_df$cons <- csq
          
        }
        current_pop_current_cat_current_type_df$cat <- cat
        colnames(current_pop_current_cat_current_type_df)[1] <- "value"
        row.names(current_pop_current_cat_current_type_df) <- NULL

        #since I want a frequecy, now I've to divide all samples count by the number of variants for that category
        tot_current_cat_var_num <- dim(current_region_current_current_pop)[1]
        current_pop_current_cat_current_type_df$freq <- current_pop_current_cat_current_type_df$value/tot_current_cat_var_num
        
        if (pop == "CARL"){
            pop <- "CAR"
        }
        assign(paste(pop,"_",cat,"_",csq,"_df",sep=""),current_pop_current_cat_current_type_df)
        assign(paste(pop,"_",cat,"_",csq,sep=""),current_pop_current_cat_current_type)
        write.table(current_pop_current_cat_current_type_df,file=paste(base_folder,pop,"_",cat,"_",csq,sep=""),sep="\t",col.names=T,quote=F,row.names=F)
      }
    }
  }
}

#box plot
###### REPLOT with ggplot
require(ggplot2)
require(reshape2)
all_cols <-col_pop(pops)

ylab <- "Number of mutations per individual"

#for ALL the variant classes together
# all_pop_merged <- NULL
# for(pop in pops){
#   current_syn <- get(paste(pop,"_all_syn_df",sep=''))
#   current_miss <- get(paste(pop,"_all_miss_df",sep=''))
#   current_merged <- rbind(current_miss,current_syn)
#   assign(paste(pop,"_all_merged_df",sep=""),current_merged)
#   all_pop_merged <- rbind(all_pop_merged,current_merged)
# }

#for splitted variant classes
shared_private_all_pop_merged <- NULL

for(pop in pops){
  all_csq_all_cat_current_pop <- NULL
  for (cat in categs){
    all_csq_current_cat_current_pop <- NULL
    for(csq in csqs){
      # if (pop == "CARL" && cat == "novel"){
      #   pop <- "CAR"
      # }
      all_obj_name <- paste(pop,"_",cat,"_",csq,"_df",sep="")
      # alt_obj_name <- paste(pop,"_alt_",cat,"_",csq,"_df",sep="")
      if(exists(all_obj_name)){
        current_csq_current_cat_current_pop <- get(all_obj_name)
        all_csq_current_cat_current_pop <- rbind(all_csq_current_cat_current_pop,current_csq_current_cat_current_pop)
      }
    }
    all_csq_all_cat_current_pop <- rbind(all_csq_all_cat_current_pop,all_csq_current_cat_current_pop)
  }
  assign(paste(pop,"_all_csq_all_cat_merged_df",sep=""),all_csq_all_cat_current_pop)
  shared_private_all_pop_merged <- rbind(shared_private_all_pop_merged,all_csq_all_cat_current_pop)
}

data_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/CONSEQUENCES/"
# write.table(shared_private_all_pop_merged,file=paste(data_folder,"shared_private_all_pop_merged.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(shared_private_all_pop_merged,file=paste(data_folder,"shared_private_all_pop_merged_others.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

shared_private_all_pop_merged$pop2 <- factor(shared_private_all_pop_merged$pop,pops)

pl <- ggplot(shared_private_all_pop_merged)
pl <- pl + geom_boxplot()
# pl <- pl + aes(x = factor(pop2), y = value, fill=pop2)
pl <- pl + aes(x = factor(pop2), y = freq, fill=pop2)
pl <- pl + ylab(ylab)
pl <- pl + xlab("")
pl <- pl + guides(fill=guide_legend(title="Cohorts"))
pl <- pl + scale_fill_manual("Cohorts", values=all_cols)
pl <- pl + theme_bw(18)
pl <- pl + facet_grid(cat~cons, scales="free")
ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_ggplot_freq.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_ggplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)


########################################################
#alternative plot with sum private and novel
pops_ingi <- sort(c("VBI","FVG","CARL"))
shared_private_plus_novel_all_pop_merged <- NULL

for(pop in pops){
  all_csq_all_cat_current_pop <- NULL
  # for (cat in categs){
    # all_csq_current_cat_current_pop <- NULL
    for(csq in csqs){
      shared_obj_name <- paste(pop,"_shared_",csq,"_df",sep="")
      private_obj_name <- paste(pop,"_private_",csq,"_df",sep="")
      novel_obj_name <- paste(pop,"_novel_",csq,"_df",sep="")
      # alt_obj_name <- paste(pop,"_alt_",cat,"_",csq,"_df",sep="")
      if(exists(shared_obj_name)){
        current_csq_current_shared_current_pop <- get(shared_obj_name)
      }
      if(exists(private_obj_name)){
        current_csq_current_private_current_pop <- get(private_obj_name)
      }
      if(exists(novel_obj_name)){
        current_csq_current_novel_current_pop <- get(novel_obj_name)
      }
      if (exists(private_obj_name) && exists(novel_obj_name)){
        #merge private and novel and sum value columns
        current_csq_current_private_novel_current_pop <- merge(current_csq_current_private_current_pop,current_csq_current_novel_current_pop[,c(1,2)],by="samples",sort=F)
        #now sum values and remove useless columns
        current_csq_current_private_novel_current_pop$value <- current_csq_current_private_novel_current_pop$value.x + current_csq_current_private_novel_current_pop$value.y
        current_csq_current_private_novel_current_pop$value.x <- NULL
        current_csq_current_private_novel_current_pop$value.y <- NULL
        current_csq_current_private_novel_current_pop$cat <- "private+novel"
      }else{
        current_csq_current_private_novel_current_pop <- NULL
      }
      #now rbind shared and private+novel
      current_csq_current_cat_current_pop <- rbind(current_csq_current_shared_current_pop,current_csq_current_private_novel_current_pop)
      all_csq_all_cat_current_pop <- rbind(all_csq_all_cat_current_pop,current_csq_current_cat_current_pop)
      # }
    }
    # all_csq_all_cat_current_pop <- rbind(all_csq_all_cat_current_pop,all_csq_current_cat_current_pop)
  # }
  assign(paste(pop,"_all_csq_all_cat_merged_df",sep=""),all_csq_all_cat_current_pop)
  shared_private_plus_novel_all_pop_merged <- rbind(shared_private_plus_novel_all_pop_merged,all_csq_all_cat_current_pop)
}


data_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/CONSEQUENCES/"
# write.table(shared_private_plus_novel_all_pop_merged,file=paste(data_folder,"shared_private_plus_novel_all_pop_merged.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
write.table(shared_private_plus_novel_all_pop_merged,file=paste(data_folder,"shared_private_plus_novel_all_pop_merged_others.txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#sort factors for plot in right order
shared_private_plus_novel_all_pop_merged$pop2 <- factor(shared_private_plus_novel_all_pop_merged$pop,pops)

pl <- ggplot(shared_private_plus_novel_all_pop_merged)
pl <- pl + geom_boxplot()
pl <- pl + aes(x = factor(pop2), y = value, fill=pop2)
pl <- pl + ylab(ylab)
pl <- pl + xlab("")
pl <- pl + guides(fill=guide_legend(title="Cohorts"))
pl <- pl + scale_fill_manual("Cohorts", values=all_cols)
pl <- pl + theme_bw(18)
pl <- pl + facet_grid(cat~cons, scales="free")
ggsave(filename=paste(data_folder,"/8b_shared_private_plus_novel_pop_conseq_carriers_ggplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

###############################################
#we should split all villages in fvg
#read keeplist, than change population name and replot
fvg_pops <- c("Erto","Illegio","Resia","Sauris")
pop_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop"
shared_private_plus_novel_all_pop_merged_fvg_split <- shared_private_plus_novel_all_pop_merged
shared_private_all_pop_merged_fvg_split <- shared_private_all_pop_merged

for(f_pop in fvg_pops){
  current_pop_list <- read.table(paste(pop_folder,"/",f_pop,"_unrelated.list",sep=""))
  colnames(current_pop_list) <- c("samples")
  if (f_pop =="Resia"){
    current_pop_list$samples <- gsub("(^\\d+)","X\\1",current_pop_list$samples)
  }
  if(f_pop=="Illegio"){ c_pop <- "FVI"}
  if(f_pop=="Resia"){ c_pop <- "FVR"}
  if(f_pop=="Erto"){ c_pop <- "FVE"}
  if(f_pop=="Sauris"){ c_pop <- "FVS"}
  shared_private_all_pop_merged_fvg_split[which(shared_private_all_pop_merged_fvg_split$samples %in% current_pop_list$samples),]$pop <- c_pop
  # shared_private_plus_novel_all_pop_merged_fvg_split[which(shared_private_plus_novel_all_pop_merged_fvg_split$samples %in% current_pop_list$samples),]$pop <- c_pop
}

all_pops <- c("CEU","TSI","CAR","VBI","FVE","FVI","FVR","FVS")
all_cols <-col_pop(all_pops)

csqs1 <- c("CHR22","Missense","Synonymous")
csqs2 <- c("CHR22","Missense","polyphen.possibly.damaging","polyphen.probably.damaging","sift.deleterious","Synonymous","polyphen.benign","sift.tolerated")
csqs3 <- c("polyphen.benign","polyphen.possibly.damaging","polyphen.probably.damaging","sift.deleterious","sift.tolerated")

shared_private_all_pop_merged_fvg_split$pop2 <- factor(shared_private_all_pop_merged_fvg_split$pop,all_pops)

#extract data for different categories
shared_private_all_pop_merged_fvg_split_csqs1 <- shared_private_all_pop_merged_fvg_split[shared_private_all_pop_merged_fvg_split$cons %in% csqs1,]
shared_private_all_pop_merged_fvg_split_csqs2 <- shared_private_all_pop_merged_fvg_split[shared_private_all_pop_merged_fvg_split$cons %in% csqs2,]
shared_private_all_pop_merged_fvg_split_csqs3 <- shared_private_all_pop_merged_fvg_split[shared_private_all_pop_merged_fvg_split$cons %in% csqs3,]

#fix the factor order
shared_private_all_pop_merged_fvg_split_csqs1$pop2 <- factor(shared_private_all_pop_merged_fvg_split_csqs1$pop,all_pops)
shared_private_all_pop_merged_fvg_split_csqs2$pop2 <- factor(shared_private_all_pop_merged_fvg_split_csqs2$pop,all_pops)
shared_private_all_pop_merged_fvg_split_csqs3$pop2 <- factor(shared_private_all_pop_merged_fvg_split_csqs3$pop,all_pops)

#fix the freq by dividing by 2
shared_private_all_pop_merged_fvg_split_csqs1$freq <- (shared_private_all_pop_merged_fvg_split_csqs1$freq)/2
shared_private_all_pop_merged_fvg_split_csqs2$freq <- (shared_private_all_pop_merged_fvg_split_csqs2$freq)/2
shared_private_all_pop_merged_fvg_split_csqs3$freq <- (shared_private_all_pop_merged_fvg_split_csqs3$freq)/2
# shared_private_plus_novel_all_pop_merged_fvg_split$pop2 <- factor(shared_private_plus_novel_all_pop_merged_fvg_split$pop,all_pops)

#first plot
pl <- ggplot(shared_private_all_pop_merged_fvg_split)
pl <- pl + geom_boxplot()
# pl <- pl + aes(x = factor(pop2), y = value, fill=pop2)
pl <- pl + aes(x = factor(pop2), y = freq, fill=pop2)
pl <- pl + ylab(ylab)
pl <- pl + xlab("")
pl <- pl + guides(fill=guide_legend(title="Cohorts"))
pl <- pl + scale_fill_manual("Cohorts", values=all_cols)
pl <- pl + theme_bw(18)
pl <- pl + facet_grid(cat~cons, scales="free")
pl <- pl + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_fvg_split_ggplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_fvg_split_ggplot_chr22.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_fvg_split_ggplot_others.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_fvg_split_ggplot_chr22_freq.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_fvg_split_ggplot_others_freq.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)

#plot for the different set of categories
ylab <- "Frequency of mutations per individual"
for(i in 1:3){
  current_to_plot <- get(paste("shared_private_all_pop_merged_fvg_split_csqs",i,sep=""))
  pl <- ggplot(current_to_plot)
  pl <- pl + geom_boxplot()
  # pl <- pl + aes(x = factor(pop2), y = value, fill=pop2)
  pl <- pl + aes(x = factor(pop2), y = freq, fill=pop2)
  pl <- pl + ylab(ylab)
  pl <- pl + xlab("")
  pl <- pl + guides(fill=guide_legend(title="Cohorts"))
  pl <- pl + scale_fill_manual("Cohorts", values=all_cols)
  pl <- pl + theme_bw(14)
  pl <- pl + facet_grid(cat~cons, scales="free")
  pl <- pl + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename=paste(data_folder,"/8b_shared_private_novel_pop_conseq_carriers_fvg_split_ggplot_freq_csqs",i,".jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
  
}

# VBI
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),])

# syn
vbi_novels_s <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),]
vbi_private_s <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),]

merged_vbi_np_s <- merge(vbi_novels_s,vbi_private_s,by="samples")
merged_vbi_np_s$np_c <- merged_vbi_np_s$value.x + merged_vbi_np_s$value.y
summary(merged_vbi_np_s)

# miss
vbi_novels_m <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),]
vbi_private_m <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),]

merged_vbi_np_m <- merge(vbi_novels_m,vbi_private_m,by="samples")
merged_vbi_np_m$np_c <- merged_vbi_np_m$value.x + merged_vbi_np_m$value.y
summary(merged_vbi_np_m)

# chr22
vbi_novels_22 <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),]
vbi_private_22 <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "VBI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),]

merged_vbi_np_22 <- merge(vbi_novels_22,vbi_private_22,by="samples")
merged_vbi_np_22$np_c <- merged_vbi_np_22$value.x + merged_vbi_np_22$value.y
summary(merged_vbi_np_22)

# CAR
car_novels <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),]
car_private <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),]

merged_car_np <- merge(car_novels,car_private,by="samples")
merged_car_np$np_c <- merged_car_np$value.x + merged_car_np$value.y
summary(merged_car_np)

car_novels_m <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),]
car_private_m <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),]

merged_car_np_m <- merge(car_novels_m,car_private_m,by="samples")
merged_car_np_m$np_c <- merged_car_np_m$value.x + merged_car_np_m$value.y
summary(merged_car_np_m)

summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),])

# chr22
car_novels_22 <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),]
car_private_22 <- shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CAR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),]

merged_car_np_22 <- merge(car_novels_22,car_private_22,by="samples")
merged_car_np_22$np_c <- merged_car_np_22$value.x + merged_car_np_22$value.y
summary(merged_car_np_22)

# FVG
summary(rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),]))

summary(rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),]))

summary(rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),]))

# syn
fvg_novels_s <- rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),])

fvg_private_s <- rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),])


merged_fvg_np_s <- merge(fvg_novels_s,fvg_private_s,by="samples")
merged_fvg_np_s$np_c <- merged_fvg_np_s$value.x + merged_fvg_np_s$value.y
summary(merged_fvg_np_s)

# miss
fvg_novels_m <- rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),])

fvg_private_m <- rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),])

merged_fvg_np_m <- merge(fvg_novels_m,fvg_private_m,by="samples")
merged_fvg_np_m$np_c <- merged_fvg_np_m$value.x + merged_fvg_np_m$value.y
summary(merged_fvg_np_m)

# chr22
fvg_novels_22 <- rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="novel" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),])

fvg_private_22 <- rbind(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVE" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVR" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),],
shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "FVS" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="private" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),])


merged_fvg_np_22 <- merge(fvg_novels_22,fvg_private_22,by="samples")
merged_fvg_np_22$np_c <- merged_fvg_np_22$value.x + merged_fvg_np_22$value.y
summary(merged_fvg_np_22)

# CEU
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CEU" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CEU" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "CEU" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),])

# TSI
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "TSI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Synonymous"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "TSI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="Missense"),])
summary(shared_private_all_pop_merged_fvg_split_csqs1[which(shared_private_all_pop_merged_fvg_split_csqs1$pop == "TSI" & shared_private_all_pop_merged_fvg_split_csqs1$cat=="shared" & shared_private_all_pop_merged_fvg_split_csqs1$cons=="CHR22"),])


#second plot
pl <- ggplot(shared_private_plus_novel_all_pop_merged_fvg_split)
pl <- pl + geom_boxplot()
pl <- pl + aes(x = factor(pop2), y = value, fill=pop2)
pl <- pl + ylab(ylab)
pl <- pl + xlab("")
pl <- pl + guides(fill=guide_legend(title="Cohorts"))
pl <- pl + scale_fill_manual("Cohorts", values=all_cols)
pl <- pl + theme_bw(18)
pl <- pl + facet_grid(cat~cons, scales="free")
pl <- pl + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename=paste(data_folder,"/8b_shared_private_plus_novel_pop_conseq_carriers_fvg_split_ggplot.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
# ggsave(filename=paste(data_folder,"/8b_shared_private_plus_novel_pop_conseq_carriers_fvg_split_ggplot_chr22.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
ggsave(filename=paste(data_folder,"/8b_shared_private_plus_novel_pop_conseq_carriers_fvg_split_ggplot_others.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)


