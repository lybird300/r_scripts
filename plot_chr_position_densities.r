#09/19/2012

#R script for plotting not overlapping snps...
#
rm(list=ls())

#read commandline args
args <- commandArgs(trailing=TRUE)
 
#Args list:
#args[[1]]: complete chr map file path
#args[[2]]: partial chr map file path
#args[[3]]: chr number

complete_chr_filepath <- args[[1]]
partial_chr_filepath <- args[[2]]
chr <- args[[3]]

#Upload chr data:
complete_chr <- read.table(complete_chr_filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="")
partial_chr <- read.table(partial_chr_filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="")

#set colnames
colnames(complete_chr) <- c('CHROM','POS','ID','AN','AC','DP')
#the partial subset only has HWE column
colnames(partial_chr) <- c('CHROM','POS','ID','AN','AC','DP','HWE')

#plot first all chr pos vs DP, then add partial chr pos vs DP

jpeg(paste("chr",chr,"not_overlapping_dp.jpg",sep="_"),width=1000, height=1000)
par(lab=c(10,10,12),pch=20,cex=2)
plot(complete_chr$POS,complete_chr$DP, main=paste("Chr",chr,sep=" "), xlab="Position",ylab="Read depth",type='p',col=c('blue'))
lines(partial_chr$POS,partial_chr$DP, type='p', col=c('red'))
dev.off()

#now calculate the AF based on AC/AN
complete_chr <- cbind(complete_chr,(complete_chr$AC/complete_chr$AN))
colnames(complete_chr) <- c('CHROM','POS','ID','AN','AC','DP','AF')

#now for partial chr
partial_chr <- cbind(partial_chr,(partial_chr$AC/partial_chr$AN))
colnames(partial_chr) <- c('CHROM','POS','ID','AN','AC','DP','HWE','AF')

#NOW the maf
#complete
complete_chr_maf <- complete_chr[which(complete_chr$AF <= 0.5),]
complete_chr_af <- complete_chr[which(complete_chr$AF > 0.5),]
complete_chr_af$AF <- (1-complete_chr_af$AF)
complete_chr_maf <- rbind(complete_chr_maf,complete_chr_af)
complete_chr_maf <- complete_chr_maf[order(complete_chr_maf$POS),]
colnames(complete_chr_maf) <- c('CHROM','POS','ID','AN','AC','DP','MAF')

rm(complete_chr_af)
gc()


#partial
partial_chr_maf <- partial_chr[which(partial_chr$AF <= 0.5),]
partial_chr_af <- partial_chr[which(partial_chr$AF > 0.5),]
partial_chr_af$AF <- (1-partial_chr_af$AF)
partial_chr_maf <- rbind(partial_chr_maf,partial_chr_af)
partial_chr_maf <- partial_chr_maf[order(partial_chr_maf$POS),]
colnames(partial_chr_maf) <- c('CHROM','POS','ID','AN','AC','DP','HWE','MAF')

rm(partial_chr_af)
gc()

#write a summary for each chr
write.table(summary(partial_chr_maf),file=paste("summary_chr",chr,"not_overlap.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)

#plot pos vs maf
jpeg(paste("chr",chr,"not_overlapping_maf.jpg",sep="_"),width=1000, height=1000)
par(lab=c(10,10,12),pch=20,cex=2)
plot(partial_chr_maf$POS,partial_chr_maf$MAF, main=paste("Chr",chr,sep=" "), xlab="Position",ylab="MAF",type='p',col=c('blue'))
dev.off()

#plot pos vs hwe
jpeg(paste("chr",chr,"not_overlapping_hwe.jpg",sep="_"),width=1000, height=1000)
par(lab=c(10,10,12),pch=20,cex=2)
plot(partial_chr_maf$POS,-log10(partial_chr_maf$HWE), main=paste("Chr",chr,sep=" "), xlab="Position",ylab="-Log(HWE pval)",type='p',col=c('green'))
dev.off()
q(save="no")
