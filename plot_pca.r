
rm(list=ls())
args=commandArgs(trailing=TRUE)
#args[[1]] = cohort
#args[[2]] = chr
#args[[3]] = mds file
#args[[4]] = filtering
twins <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples_pheno.centre", sep=" ",header=F)
alspac <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_ALSPAC_samples_pheno.centre", sep=" ",header=F)
joint <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples_pheno.centre_pop", sep=" ",header=F)

args

args[[1]]
args[[2]]
args[[3]]
args[[4]]

mds_file <- args[[3]]
mds <- read.table(mds_file,sep="\t",header=T)

#mds <- read.table("chr20.mds.mod",sep="\t",header=T)
colnames(twins) <- c("FID","IID","SEQ")
colnames(alspac) <- c("FID","IID","SEQ")
colnames(joint) <- c("FID","IID","SEQ","POP")
cohort <- args[[1]]

if ( cohort == "ALSPAC" ){
merged_new <- merge(mds,alspac,by.x="IID",by.y="IID",all.x)
}

if ( cohort == "TWINSUK" ){
merged_new <- merge(mds,twins,by.x="IID",by.y="IID",all.x)
}
if (cohort == "JOINT" ){
merged_new <- merge(mds,joint,by.x="IID",by.y="IID",all.x)
}

filter <- args[[4]]
#filter <- "5e-6"
#filter <- "1e-5"
#filter <- "UNFILTERED"
#filter <- "5e-8"
#filter <- "MAF_ge_01"
if (filter == "UNFILTERED"){
type <- "UNFILTERED"
}else{
type <- "FILTERED"
}
chr <- args[[2]]
basefolder <- "/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/"

###########SINGLE CHR###############
xmin <- min(merged_new$C1)-0.001
ymin <- min(merged_new$C2)+0.001

pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_",chr,".pdf",sep=""),width=11.7, height=8.3)
par(lab=c(4,4,6),mfrow=c(2,2))
plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C2,col="red")
legend(xmin,ymin,legend=c("BGI","SC"),fill=c("black","red"))

plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C3,xlim=c(min(merged_new$C1),max(merged_new$C1)),ylim=c(min(merged_new$C3),max(merged_new$C3)),main=paste("PCA C1 vs C3, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C3,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C4,xlim=c(min(merged_new$C1),max(merged_new$C1)),ylim=c(min(merged_new$C4),max(merged_new$C4)),main=paste("PCA C1 vs C4, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C4,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C2,merged_new[which(merged_new$SEQ == 1),]$C3,xlim=c(min(merged_new$C2),max(merged_new$C2)),ylim=c(min(merged_new$C3),max(merged_new$C3)),main=paste("PCA C2 vs C3, pruned ",cohort,"  - ",type,sep=""),xlab="C2",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2),]$C2,merged_new[which(merged_new$SEQ == 2),]$C3,col="red")
dev.off()
