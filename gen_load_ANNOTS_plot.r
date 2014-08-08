#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"




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
