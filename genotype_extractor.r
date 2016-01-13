#script to extract genotypes for samples by snp

library(buRlo)
file_path <- getwd()
file_path <- setwd("/home/cocca/breast_cancer/PAPER_RESULTS/700K/")
villages <- c("clauzetto","erto","illegio","resia","sauris","smc")

for (village in villages) {
current_ped <- paste(file_path,village,".ped")
current_map <- paste(file_path,village,".map")
current_out <- paste(file_path,village,".ped")

convert.snp.ped(pedfile, mapfile, outfile, format = "premakeped", traits = 1, strand = "u", bcast = 10000000, wslash=FALSE, mapHasHeaderLine=FALSE)
}

