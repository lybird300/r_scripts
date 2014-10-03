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

# pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
# chr <- "10"
# pop <- "CARL"
#INGI
# pops <- c("VBI","FVG","CARL")
# categs <- c("private","shared","all")
#OUTBRED
pops <- c("CEU","TSI","VBI","FVG","CARL")
categs <- c("all")
csqs <- c("miss","syn")

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
        }else{
          # 1.INGI_chr1.merged_maf.tab.gz.VBI.shared.tab.gz.miss.regions.vcf.gz.tab
          current_region_file <- paste(base_folder,"/",chr,".INGI_chr",chr,".merged_maf.tab.gz.",pop,".",cat,".tab.gz.",csq,".regions.vcf.gz.tab",sep="")
        }
        print(current_region_file)
        current_region_current_chr_current_pop <- read.table(current_region_file,header=T)
        #we need to convert our genotypes, now!
        for(j in 7:length(colnames(current_region_current_chr_current_pop))){
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "0|0"] <- "0" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "0|1"] <- "1" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "1|0"] <- "1" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "1|1"] <- "1" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "2|0"] <- "1" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "0|2"] <- "1" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "2|1"] <- "1" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "1|2"] <- "1" 
          levels(current_region_current_chr_current_pop[,j])[levels(current_region_current_chr_current_pop[,j]) == "2|2"] <- "1" 
          current_region_current_chr_current_pop[,j] <- as.numeric(as.character(current_region_current_chr_current_pop[,j]))
        }
        # EGAN00001172162
        current_region_current_chr_current_pop_samples <- current_region_current_chr_current_pop[,7:length(colnames(current_region_current_chr_current_pop))]
        current_region_current_current_pop <- rbind(current_region_current_current_pop,current_region_current_chr_current_pop_samples)
      }
      current_pop_current_cat_current_type <- apply(current_region_current_current_pop,2,sum)
      assign(paste(pop,"_",cat,"_",csq,sep=""),current_pop_current_cat_current_type)
      write.table(current_pop_current_cat_current_type,file=paste(base_folder,pop,"_",cat,"_",csq,sep=""),sep="\t",col.names=T,quote=F,row.names=F)
    }
  }
}
