#script to read and plot batch effects for sequencing centre
rm(list=ls())
alspac <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_ALSPAC_samples_pheno.centre", sep=" ",header=F)
joint <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples_pheno.centre_pop", sep=" ",header=F)
twins <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples_pheno.centre", sep=" ",header=F)
mds <- read.table("all_chr_merged.mds",header=T)
#mds <- read.table("chr20.mds.mod",sep="\t",header=T)
colnames(twins) <- c("FID","IID","SEQ")
colnames(alspac) <- c("FID","IID","SEQ")
colnames(joint) <- c("FID","IID","SEQ","POP")
merged_new <- merge(mds,alspac,by.x="IID",by.y="IID",all.x)
merged_new <- merge(mds,twins,by.x="IID",by.y="IID",all.x)
merged_new <- merge(mds,joint,by.x="IID",by.y="IID",all.x)

#set a formula for regression association
# formula_glm <- "SEQ ~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10"
# res <- glm(data=merged_new,formula_glm)

# complete_samples <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples.centre",sep="\t",header=F)
# old_mds <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/old_overall.mds",sep=" ",header=T)

#filter
filter <- "5e-6"
filter <- "1e-5"
filter <- "1e-5_RR"
filter <- "1e-2_WR"
filter <- "1e-2"
filter <- "1e-2_RR"
filter <- "1e-2_gwa"
filter <- "1e-2_gwa_5"
filter <- "1e-2_gwa_1_5"
filter <- "UNFILTERED"
filter <- "5e-8"
filter <- "MAF_ge_01"
filter <- "KLAUDIA_MAF_ge_01"
filter <- "KLAUDIA"
#cohort
cohort <- "TWINS"
cohort <- "JOINT"
cohort <- "ALSPAC"
#type
type <- "FILTERED"
type <- "UNFILTERED"

#geno set
geno <- "PRUNED"
geno <- "PRUNED_maf_ge01"
geno <- "UNPRUNED"
geno <- "UNPRUNED_maf_ge01"
chr <- "20"
#basefolder <- "/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/"
basefolder <- "/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/PLOTS/UNPRUNED/"
basefolder <- "/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/PLOTS/PRUNED/"


##### SINGLE POPULATION ###########
###########SINGLE CHR###############
xmin <- min(merged_new$C1)-0.001
ymin <- min(merged_new$C2)+0.0015
# pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_",chr,".pdf",sep=""),width=11.7, height=8.3)
jpeg(paste(basefolder,"PCA_",cohort,"_batch_",geno,"_",type,"_",filter,"_",chr,".jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(4,4,6),mfrow=c(2,3))
plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, ",geno," ",cohort," chr",chr,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C2,col="red")
legend(xmin,ymin,legend=c("BGI","SC"),fill=c("black","red"))

plot(merged_new[which(merged_new$SEQ == 1),]$C3,merged_new[which(merged_new$SEQ == 1),]$C4,xlim=c(min(merged_new$C3),max(merged_new$C3)),ylim=c(min(merged_new$C4),max(merged_new$C4)),main=paste("PCA C3 vs C4, ",geno," ",cohort," chr",chr,"  - ",type,sep=""),xlab="C3",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2),]$C3,merged_new[which(merged_new$SEQ == 2),]$C4,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C5,merged_new[which(merged_new$SEQ == 1),]$C6,xlim=c(min(merged_new$C5),max(merged_new$C5)),ylim=c(min(merged_new$C6),max(merged_new$C6)),main=paste("PCA C5 vs C6, ",geno," ",cohort," chr",chr,"  - ",type,sep=""),xlab="C5",ylab="C6")
points(merged_new[which(merged_new$SEQ == 2),]$C5,merged_new[which(merged_new$SEQ == 2),]$C6,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C7,merged_new[which(merged_new$SEQ == 1),]$C8,xlim=c(min(merged_new$C7),max(merged_new$C7)),ylim=c(min(merged_new$C8),max(merged_new$C8)),main=paste("PCA C7 vs C8, ",geno," ",cohort," chr",chr,"  - ",type,sep=""),xlab="C7",ylab="C8")
points(merged_new[which(merged_new$SEQ == 2),]$C7,merged_new[which(merged_new$SEQ == 2),]$C8,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C9,merged_new[which(merged_new$SEQ == 1),]$C10,xlim=c(min(merged_new$C9),max(merged_new$C9)),ylim=c(min(merged_new$C10),max(merged_new$C10)),main=paste("PCA C9 vs C10, ",geno," ",cohort," chr",chr,"  - ",type,sep=""),xlab="C9",ylab="C10")
points(merged_new[which(merged_new$SEQ == 2),]$C9,merged_new[which(merged_new$SEQ == 2),]$C10,col="red")

dev.off()

########### GWA ###############
#legend position:
xmin <- min(merged_new$C1)-0.001
ymin <- min(merged_new$C2)+0.001

pdf(paste(basefolder,"PCA_",cohort,"_batch_",geno,"_",type,"_",filter,".pdf",sep=""),width=11.7, height=8.3)
par(lab=c(4,4,6),mfrow=c(2,2))
plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, ",geno," ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C2,col="red")
legend(xmin,ymin,legend=c("BGI","SC"),fill=c("black","red"))

plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C3,xlim=c(min(merged_new$C1),max(merged_new$C1)),ylim=c(min(merged_new$C3),max(merged_new$C3)),main=paste("PCA C1 vs C3, ",geno," ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C3,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C4,xlim=c(min(merged_new$C1),max(merged_new$C1)),ylim=c(min(merged_new$C4),max(merged_new$C4)),main=paste("PCA C1 vs C4, ",geno," ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C4,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C2,merged_new[which(merged_new$SEQ == 1),]$C3,xlim=c(min(merged_new$C2),max(merged_new$C2)),ylim=c(min(merged_new$C3),max(merged_new$C3)),main=paste("PCA C2 vs C3, ",geno," ",cohort,"  - ",type,sep=""),xlab="C2",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2),]$C2,merged_new[which(merged_new$SEQ == 2),]$C3,col="red")
dev.off()


###########################################################################
# consecutive components
xmax <- max(merged_new$C1)-0.001
ymax <- max(merged_new$C2)+0.001

# pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_gwa_cons.pdf",sep=""),width=11.7, height=8.3)
jpeg(paste(basefolder,"PCA_",cohort,"_batch_",geno,"_",type,"_",filter,"_cons.jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(4,4,6),mfrow=c(2,3))
plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, ",geno," ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C2,col="red")
legend(xmax,ymax,legend=c("BGI","SC"),fill=c("black","red"))

plot(merged_new[which(merged_new$SEQ == 1),]$C3,merged_new[which(merged_new$SEQ == 1),]$C4,xlim=c(min(merged_new$C3),max(merged_new$C3)),ylim=c(min(merged_new$C4),max(merged_new$C4)),main=paste("PCA C3 vs C4, ",geno," ",cohort,"  - ",type,sep=""),xlab="C3",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2),]$C3,merged_new[which(merged_new$SEQ == 2),]$C4,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C5,merged_new[which(merged_new$SEQ == 1),]$C6,xlim=c(min(merged_new$C5),max(merged_new$C5)),ylim=c(min(merged_new$C6),max(merged_new$C6)),main=paste("PCA C5 vs C6, ",geno," ",cohort,"  - ",type,sep=""),xlab="C5",ylab="C6")
points(merged_new[which(merged_new$SEQ == 2),]$C5,merged_new[which(merged_new$SEQ == 2),]$C6,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C7,merged_new[which(merged_new$SEQ == 1),]$C8,xlim=c(min(merged_new$C7),max(merged_new$C7)),ylim=c(min(merged_new$C8),max(merged_new$C8)),main=paste("PCA C7 vs C8, ",geno," ",cohort,"  - ",type,sep=""),xlab="C7",ylab="C8")
points(merged_new[which(merged_new$SEQ == 2),]$C7,merged_new[which(merged_new$SEQ == 2),]$C8,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C9,merged_new[which(merged_new$SEQ == 1),]$C10,xlim=c(min(merged_new$C9),max(merged_new$C9)),ylim=c(min(merged_new$C10),max(merged_new$C10)),main=paste("PCA C9 vs C10, ",geno," ",cohort,"  - ",type,sep=""),xlab="C9",ylab="C10")
points(merged_new[which(merged_new$SEQ == 2),]$C9,merged_new[which(merged_new$SEQ == 2),]$C10,col="red")

dev.off()
######################################################################################

#SEQ=1 -> BGI
#SEQ=2 -> SANGER
#POP=1 -> ALSPAC
#POP=2 -> TWINS
########### GWA JOINT 3621 samples ###############
##PCA
xmin <- min(merged_new$C1)-0.0025
ymax <- max(merged_new$C2)+0.0025

# pdf(paste(basefolder,"PCA_",cohort,"_batch_",geno,"_",type,"_",filter,".pdf",sep=""),width=11.7, height=8.3)
jpeg(paste(basefolder,"PCA_",cohort,"_batch_",geno,"_",type,"_",filter,".jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(2,2,4),mfrow=c(2,2))
plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,xlim=c(min(merged_new$C1)-0.005,max(merged_new$C1)+0.005),ylim=c(min(merged_new$C2)-0.005,max(merged_new$C2)+0.005),main=paste("PCA C1 vs C2, ",geno," ",cohort," - ",filter,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,col="green")
legend(xmin,ymax,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,xlim=c(min(merged_new$C1)-0.003,max(merged_new$C1)+0.003),ylim=c(min(merged_new$C3)-0.003,max(merged_new$C3)+0.003),main=paste("PCA C1 vs C3, ",geno," ",cohort," - ",filter,sep=""),xlab="C1",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C4,xlim=c(min(merged_new$C1)-0.003,max(merged_new$C1)+0.003),ylim=c(min(merged_new$C4)-0.003,max(merged_new$C4)+0.003),main=paste("PCA C1 vs C4, ",geno," ",cohort," - ",filter,sep=""),xlab="C1",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C4,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C4,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C4,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,xlim=c(min(merged_new$C2)-0.003,max(merged_new$C2)+0.003),ylim=c(min(merged_new$C3)-0.003,max(merged_new$C3)+0.003),main=paste("PCA C2 vs C3, ",geno," ",cohort," - ",filter,sep=""),xlab="C2",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,col="green")
dev.off()

###################################################
##Consecutive PCA
xmin <- min(merged_new$C1)-0.0025
ymax <- max(merged_new$C2)+0.0025

# pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_gwa_cons.pdf",sep=""),width=11.7, height=8.3)
jpeg(paste(basefolder,"PCA_",cohort,"_batch_",geno,"_",type,"_",filter,"_cons.jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(2,2,4),mfrow=c(2,3))
plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,xlim=c(min(merged_new$C1)-0.0025,max(merged_new$C1)+0.0025),ylim=c(min(merged_new$C2)-0.0025,max(merged_new$C2)+0.0025),main=paste("PCA C1 vs C2, ",geno," ",cohort," - ",filter,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,col="green")
legend(xmin,ymax,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C4,xlim=c(min(merged_new$C3)-0.0025,max(merged_new$C3)+0.0025),ylim=c(min(merged_new$C4)-0.0025,max(merged_new$C4)+0.0025),main=paste("PCA C3 vs C4, ",geno," ",cohort," - ",filter,sep=""),xlab="C3",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C4,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C4,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C4,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C5,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C6,xlim=c(min(merged_new$C5)-0.0025,max(merged_new$C5)+0.0025),ylim=c(min(merged_new$C6)-0.0025,max(merged_new$C6)+0.0025),main=paste("PCA C5 vs C6, ",geno," ",cohort," - ",filter,sep=""),xlab="C5",ylab="C6")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C5,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C6,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C5,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C6,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C5,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C6,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C7,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C8,xlim=c(min(merged_new$C7)-0.0025,max(merged_new$C7)+0.0025),ylim=c(min(merged_new$C8)-0.0025,max(merged_new$C8)+0.0025),main=paste("PCA C7 vs C8, ",geno," ",cohort," - ",filter,sep=""),xlab="C7",ylab="C8")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C7,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C8,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C7,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C8,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C7,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C8,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C9,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C10,xlim=c(min(merged_new$C9)-0.0025,max(merged_new$C9)+0.0025),ylim=c(min(merged_new$C10)-0.0025,max(merged_new$C10)+0.0025),main=paste("PCA C9 vs C10, ",geno," ",cohort," - ",filter,sep=""),xlab="C9",ylab="C10")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C9,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C10,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C9,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C10,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C9,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C10,col="green")

dev.off()


############################ JOINT by CHR
xmin <- min(merged_new$C1)-0.0055
ymax <- max(merged_new$C2)+0.003

jpeg(paste(basefolder,"PCA_",cohort,"_batch_",geno,"_",type,"_",filter,"_",chr,"_cons.jpg",sep=""),width=1754, height=1024,pointsize = 20)
par(lab=c(2,2,4),mfrow=c(2,3))
plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,xlim=c(min(merged_new$C1)-0.005,max(merged_new$C1)+0.005),ylim=c(min(merged_new$C2)-0.005,max(merged_new$C2)+0.005),main=paste("PCA C1 vs C2, ",geno," ",cohort," - ",filter,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,col="green")
legend(xmin,ymax,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C4,xlim=c(min(merged_new$C3)-0.003,max(merged_new$C3)+0.003),ylim=c(min(merged_new$C4)-0.003,max(merged_new$C4)+0.003),main=paste("PCA C3 vs C4, ",geno," ",cohort," - ",filter,sep=""),xlab="C3",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C4,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C4,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C4,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C5,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C6,xlim=c(min(merged_new$C5)-0.003,max(merged_new$C5)+0.003),ylim=c(min(merged_new$C6)-0.003,max(merged_new$C6)+0.003),main=paste("PCA C5 vs C6, ",geno," ",cohort," - ",filter,sep=""),xlab="C5",ylab="C6")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C5,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C6,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C5,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C6,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C5,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C6,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C7,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C8,xlim=c(min(merged_new$C7)-0.003,max(merged_new$C7)+0.003),ylim=c(min(merged_new$C8)-0.003,max(merged_new$C8)+0.003),main=paste("PCA C7 vs C8, ",geno," ",cohort," - ",filter,sep=""),xlab="C7",ylab="C8")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C7,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C8,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C7,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C8,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C7,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C8,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C9,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C10,xlim=c(min(merged_new$C9)-0.003,max(merged_new$C9)+0.003),ylim=c(min(merged_new$C10)-0.003,max(merged_new$C10)+0.003),main=paste("PCA C9 vs C10, ",geno," ",cohort," - ",filter,sep=""),xlab="C9",ylab="C10")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C9,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C10,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C9,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C10,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C9,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C10,col="green")
dev.off()

#correlation between components after filtering
corr_matrix <- NULL
for(i in 4:length(colnames(mds_unfiltered))){
  corr_p <-NULL
  for(j in 4:length(colnames(mds_unfiltered))){
    cort <- cor.test(mds_unfiltered[,i],mds_filtered_5e8[,j])
    corr_p <- c(corr_p,cort$p.value)
  }
  assign(paste("corr_p",i,sep=""),corr_p)
  corr_matrix <- as.matrix(rbind(corr_matrix,get(paste("corr_p",i,sep=""))))
}

corr_matrix <- NULL
for(i in 4:length(colnames(mds_unfiltered))){
  corr_p <-NULL
  for(j in 4:length(colnames(mds_unfiltered))){
    cort <- cor.test(mds_unfiltered[,i],mds_filtered_5e6[,j])
    corr_p <- c(corr_p,cort$p.value)
  }
  assign(paste("corr_p",i,sep=""),corr_p)
  corr_matrix <- as.matrix(rbind(corr_matrix,get(paste("corr_p",i,sep=""))))
}

corr_matrix <- NULL
for(i in 4:length(colnames(mds_unfiltered))){
  corr_p <-NULL
  for(j in 4:length(colnames(mds_unfiltered))){
    cort <- cor.test(mds_unfiltered[,i],mds_filtered_1e5[,j])
    corr_p <- c(corr_p,cort$p.value)
  }
  assign(paste("corr_p",i,sep=""),corr_p)
  corr_matrix <- as.matrix(rbind(corr_matrix,get(paste("corr_p",i,sep=""))))
}

#
corr_p <-NULL
for(i in 4:13){
    cort <- cor.test(merged_new$SEQ,merged_new[,i])
  assign(paste("corr_p",i,sep=""),cort$p.value)
    print(cort$p.value)
    corr_p <- c(corr_p,cort$p.value)
}

corr_p <-NULL
for(i in 4:13){
    cort <- cor.test(merged_new$POP,merged_new[,i])
  assign(paste("corr_p",i,sep=""),cort$p.value)
    print(cort$p.value)
    corr_p <- c(corr_p,cort$p.value)
}


#calculate the variance explained from the pca components
#Klaudia way:
f2 <- read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/FILTERED_1e2_gwa/all_chr_merged.mds", header=T, stringsAsFactors=F)
points2 <- f2[,4:13]
## Eigen values
eig2 <- rep(0,10)
for (i in 1:10) {
  eig2[i] <- sum(points2[,i]^2)
}
p <- eig2/(sum(eig2))

#Max way:
mzzzzzz_f2 <- read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/UNFILTERED_1/all_chr_merged.mds", header=T, stringsAsFactors=F)

m_f2 <- read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/FILTERED_1e2_gwa/all_chr_merged.mds", header=T, stringsAsFactors=F)
m_points2 <- m_f2[,4:13]
m_cov_mat_pca <- cov(m_points2)
m_eig_pca <- eigen(m_cov_mat_pca)
tot_var <- (sum(m_eig_pca$values))
m_p_pca <- m_eig_pca$values/(sum(m_eig_pca$values))


var.mds <- data.frame(apply(m_f2[, 4:13], 2,function(x){c(variance = var(x, na.rm = TRUE))}))
names(var.mds) <- "varC"

###################################################################################################################################################################################
##### Distance matrix from Plink #####

f <- read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/UNFILTERED_1/data.genome", header=T, stringsAsFactors=F)
f <- read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/FILTERED_1e2_gwa/data.genome", stringsAsFactors=F, header=T)   ## 1,537,551
## remove header lines
ind <- which(f[,1]=="FID1")
f <- f[-ind,]                 ## 1,537,381 
## IDs
id1 <- unique(f$IID1)
id2 <- unique(f$IID2)
setdiff(id1,id2)   ## QTL190044
setdiff(id2,id1)   ## UK10K_TW5120489
which(id1=="QTL190044")
which(id2=="UK10K_TW5120489")
id <- sort(unique(c(id1,id2)))
## find indices
indr <- sapply(f$IID1, function(v) which(v==id))
indc <- sapply(f$IID2, function(v) which(v==id))
dat <- cbind(as.numeric(indr), as.numeric(indc), as.numeric(f$DST))
m <- matrix(0, nrow=length(id), ncol=length(id))
for (i in 1:nrow(dat))  {
  print(i)
  x <- dat[i,1]
  y <- dat[i,2]
  m[x,y] <- dat[i,3]
}
fdat <- data.frame(m, row.names=id)
names(fdat) <- id
## change into full symmetric matrix
ffdat <- fdat + t(fdat)
#write.table(ffdat, "/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/analysis/dist_matrix_twins.txt", quote=F, sep="\t")

## alternatively
m <- matrix(0, nrow=length(id), ncol=length(id))
rownames(m) <- colnames(m) <- id
for (i in 1:nrow(f))  {
  if (i%%1000==0)  print(i)
  m[f[i,1],f[i,3]] <- f$DST[i]
}
nm <- apply(m, 2, as.numeric)
fdat <- nm + t(nm)
rownames(fdat) <- colnames(fdat)
#write.table(fdat, "/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/analysis/dist_matrix_twins_v2.txt", quote=F, sep="\t")

## replaces double for-loop
d <- matrix(0,length(id),length(id))
rownames(d) <- colnames(d) <- id
d[cbind(f$ID1,f$ID2)] <- f$DIST  ## NEW: assign all values in one go
d <- d + t(d)                    ## symmetrize                              

## test ----------
(m <- matrix(1:25,5,5))
rownames(m) <- colnames(m) <- c('A','B','C','D','E')
i <- c('A','E')
j <- c('C','D')
(inds <- cbind(i,j))
m[inds] <- c(-1,-2)
m
## end test -----------



##### Run MDS #####

f <- read.delim("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/analysis/dist_matrix_twins.txt")         ## 1754 x 1754
## multi-dimensional scaling
mds <- cmdscale(d=ffdat, k=10, eig=TRUE, x.ret=T)
sum(mds$points[,1]^2)   ## == eig[1]
sum(mds$points[,1]*mds$points[,2])    ## 0 because orthogonal
sum(mds$points[,2]*mds$points[,3])    ## 0 because orthogonal

## what happens is the following
## obtain point correlation matrix from distance matrix:
A <- -0.5 * f^2
## centering matrix H = I - 1/n * J
## can be used to simulate the effect of column centering 
## (shift each coordinate to mean 0) on covariance matrix;
## ie no need to go back to original centered coordinates for covar
n <- nrow(f)
H <- diag(1,n) - matrix(1,n,n)/n
B <- H %*% A %*% H

## all these result in the same eigenvalues
eig <- eigen(B)$val

## variance explained
p <- eig/sum(eig)
mds.eig <- mds$eig/sum(mds$eig)    ## same as eig/sum(eig)

#pdf("/lustre/scratch113/projects/uk10k/users/kw8/plots3/barplot_twins_var_explained.pdf", width=5.5, height=5, pointsize=11)
barplot(mds.eig[1:10], col="red", ylab="variance explained", names=paste("C",1:10,sep=""), xlab="components", cex.lab=1.3)
dev.off()
 
#pdf("/lustre/scratch113/projects/uk10k/users/kw8/plots3/plot_twins_var_explained.pdf", width=5, height=5, pointsize=11)
plot(mds.eig, col="red", ty="b", pch=18, ylab="variance explained", xlab="components", cex.lab=1.3)
dev.off()


##### MDS output from Plink #####

f2 <- read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/GW/FILTERED_1e2_gwa/all_chr_merged.mds", header=T, stringsAsFactors=F)
points2 <- f2[,4:13]
## Eigen values
eig2 <- rep(0,10)
for (i in 1:10) {
  eig2[i] <- sum(points2[,i]^2)
}
p <- eig2/(sum(eig2))

#pdf("/lustre/scratch113/projects/uk10k/users/kw8/plots3/barplot_twins_var_explained_plink.pdf", width=5.5, height=5, pointsize=11)
barplot(p, col="red", ylab="variance explained", names=paste("C",1:10,sep=""), cex.lab=1.3, main="Variance explained (MDS in Plink)")
dev.off()