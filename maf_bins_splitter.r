#function to count how many sites we have for each maf class
#Each bin includes the upper bound, so if we have:
#maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
#we will have numbers for:
#N_SNPS for freq in (0,0.01]
#N_SNPS for freq in (0.01,0.02]
#N_SNPS for freq in (0.02,0.05]
#N_SNPS for freq in (0.05,0.1]
#N_SNPS for freq in (0.1,0.2]
#N_SNPS for freq in (0.2,0.3]
#N_SNPS for freq in (0.3,0.4]
#N_SNPS for freq in (0.4,0.5]

split_bins <- function(maf_classes,snp_set,population,col_name){

  #snp_set <- dataframe of snps. We assume there is a column named MAF
  #
  column <- which(colnames(snp_set) %in% col_name)
  
  total_maf_count <- NULL
  # snp_set$MAF <- as.numeric(as.character(snp_set$MAF))
  snp_set[,column] <- as.numeric(as.character(snp_set[,column]))
  #se ser as a character a column CHROM if it exists
  if ("CHROM" %in% colnames(snp_set)){
  	snp_set$CHROM <- as.character(snp_set$CHROM)
  }
  #we want to set D_MAF column to nomeric, if this column exists...
  if ("D_MAF" %in% colnames(snp_set)){
	snp_set$D_MAF <- as.numeric(as.character(snp_set$D_MAF))
  }
  #count how many snps we have in each maf range
  for(i in 2:(length(maf_classes) + 1)){
    class_maf_name <- paste(population,'_maf_lte_',maf_classes[i],sep='')
    
    if (maf_classes[i - 1] == maf_classes[length(maf_classes)]){
      if (maf_classes[length(maf_classes)] != 0.5) {
        class_maf_count <- snp_set[snp_set[,column] > maf_classes[i-1],]
        class_maf_name <- paste(population,'_maf_gte_',maf_classes[i-1],sep='')
      }else{ break}
    }else{
      if(maf_classes[i] == maf_classes[2]){
        class_maf_count <- snp_set[snp_set[,column] <= maf_classes[i],]
      }else {
        class_maf_count <- snp_set[snp_set[,column] > maf_classes[i-1],]
        class_maf_count <- class_maf_count[class_maf_count[,column] <= maf_classes[i],]
      }
    }

    assign(paste(class_maf_name,'summary',sep='_'), summary(class_maf_count))
    
    assign(class_maf_name, length(class_maf_count[,column]))
    total_maf_count <- cbind(total_maf_count, get(class_maf_name))
    #if(chr != NULL){
    #  write.table(summary(class_maf_count),file=paste(class_maf_name,'chr',chr,'summary.txt',sep='_'), sep="\t", row.names=FALSE, col.names=TRUE, quote=F)
    #}else{
      write.table(summary(class_maf_count),file=paste(class_maf_name,'summary.txt',sep='_'), sep="\t", row.names=FALSE, col.names=TRUE, quote=F)
      write.table(class_maf_count,file=paste(class_maf_name,'table.txt',sep='_'), sep="\t", row.names=FALSE, col.names=TRUE, quote=F)

    #}

    rm(class_maf_count)
    gc()
    gc()
    gc()
    gc()
  }
gc()

total_maf_count <- as.data.frame(total_maf_count)
maf_classes <- gsub("0.","lt",as.character(maf_classes),fixed=T)
colnames(total_maf_count) <- as.character(maf_classes[2:length(maf_classes)])

return (total_maf_count)
}

