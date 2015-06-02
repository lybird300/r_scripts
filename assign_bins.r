#function to use with apply working on rows,to add a class based on numerical comparison

assign_bins <- function(data,col_name,bin_classes){

  #snp_set : dataframe of snps. We assume there is a column named MAF
  #
  column <- which(colnames(data) %in% col_name)

  data[,column] <- as.numeric(as.character(data[,column]))
  
  #return the cerresponding bin for each bin class
  for (i in 1:length(bin_classes)){
    if ( i == 1 ){
      data[data[,column] <= bin_classes[i],]
      bin <- bin_classes[i]
    }else{
      data[data[,column] > bin_classes[i] & data[,column] <= bin_classes[i-1],]
      bin <- bin_classes[i]
    }
  } 
return (bin)
}

