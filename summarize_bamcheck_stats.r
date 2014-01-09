#Summarize stats from bamcheck for single bam file

#read commandline args
args <- commandArgs(trailing=TRUE)


#Args:
#args[[1]]= input file
#Every file name is intended to be used with ABSOLUTE PATH

input_file <- args[[1]]

#start with indels plot. Plot indels per cycle count for fwd and rev
#example commmands
basepath <- "/nfs/users/nfs_m/mc14/lustre109_home/INGI_FVG/BAM_QC"
bgi_path <- "BGI/STATS/BAMCHECK_STATS/TEST"

#read_indel <- read.table(paste(basepath,bgi_path,"591041.dedup.realn.recal.bam.bamchek.stats.indel_per_cycle",sep="/"),sep="\t",header=TRUE)
read_indel <- read.table(paste(input_file,"indel_per_cycle",sep="."),sep="\t",header=TRUE)

min_values_1 <- c(min(read_indel$insertion_fwd_n),min(read_indel$insertion_rev_n),min(read_indel$deletion_fwd_n),min(read_indel$deletion_rev_n))
max_values_1 <- c(max(read_indel$insertion_fwd_n),max(read_indel$insertion_rev_n),max(read_indel$deletion_fwd_n),max(read_indel$deletion_rev_n))

jpeg(paste("Indel_count_",basename(input_file),".jpg",sep=""),height=1000,width=1000)
plot(read_indel$cycle,read_indel$insertion_fwd_n,type="l",ylim=c(min(min_values_1),max(max_values_1)), main=paste("Indel count",basename(input_file),sep=" "),xlab="Read cycle",ylab="Indels count",col="red")
lines(read_indel$cycle,read_indel$insertion_rev_n,col="black" )
lines(read_indel$cycle,read_indel$deletion_fwd_n,col="green" )
lines(read_indel$cycle,read_indel$deletion_rev_n,col="blue" )
legend(20,max(max_values_1),c("Insertion fwd","Insertion rev","Deletion fwd","Deletion rev"),lty=c(1,1),col=c("red","black","green","blue"))
grid(5,5)
dev.off()

#now read indel count vs length
#read_indel_count <- read.table(paste(basepath,bgi_path,"591041.dedup.realn.recal.bam.bamchek.stats.indel_dist",sep="/"),sep="\t",header=TRUE)
read_indel_count <- read.table(paste(input_file,"indel_dist",sep="."),sep="\t",header=TRUE)
read_indel_count <- as.data.frame(cbind(read_indel_count,indel_ratio=read_indel_count$insertion_n/read_indel_count$deletion_n))

min_values_2 <- c(min(read_indel_count$insertion_n),min(read_indel_count$deletion_n))
max_values_2 <- c(max(log(read_indel_count$insertion_n)),max(log(read_indel_count$deletion_n)))

#and plot
jpeg(paste("Indel_length_",basename(input_file),".jpg",sep=""),height=1000,width=1000)
plot(read_indel_count$length,log(read_indel_count$insertion_n),type="l",ylim=c(min(min_values_2),max(max_values_2)), main=paste("Indel count",basename(input_file),sep=" "),xlab="Length",ylab="Indels count (log)",col="red")
lines(read_indel_count$length,log(read_indel_count$deletion_n),col="black" )
#lines(read_indel_count$length,read_indel_count$indel_ratio,col="green" )
#legend(80,max(max_values_2),c("Insertion","Deletion","Insertion/Deletion ratio"),lty=c(1,1),col=c("red","black","green"))
legend(80,max(max_values_2),c("Insertion","Deletion"),lty=c(1,1),col=c("red","black"))
grid(5,5)
dev.off()

#read ATGC counts
read_ATGC <- read.table(paste(input_file,"atgc",sep="."),sep="\t",header=TRUE)

#and plot
jpeg(paste("ATGC_",basename(input_file),".jpg",sep=""),height=1000,width=1000)
plot(read_ATGC$cycle,read_ATGC$A_count,type="l", ylim=c(0,100),main=paste("ATGC content",basename(input_file),sep=" "),xlab="Read cycle",ylab="Base Content [%]",col="green")
lines(read_ATGC$cycle,read_ATGC$C_count,col="red" )
lines(read_ATGC$cycle,read_ATGC$G_count,col="black" )
lines(read_ATGC$cycle,read_ATGC$T_count,col="blue" )
legend(80,80,c("A","C","G","T"),lty=c(1,1),col=c("green","red","black","blue"))
grid(5,5)
dev.off

#read the coverage distribution plot
read_coverage <- read.table(paste(input_file,"cov_dist",sep="."),sep="\t",header=TRUE)
colnames(read_coverage) <- c("cov_bins","cov_bin","reads")
jpeg(paste("coverage_",basename(input_file),".jpg",sep=""),height=1000,width=1000)
plot(read_coverage$cov_bin,read_coverage$reads,type="l",main=paste("Coverage ",basename(input_file),sep=" "),xlab="Coverage bin",ylab="Reads number")
grid(5,5)
dev.off()


#write.table(summary(cov_info),file=paste(input_file,'summary.txt',sep='_'), sep="\t", row.names=FALSE, col.names=TRUE, quote=F)



#q(save="yes")
