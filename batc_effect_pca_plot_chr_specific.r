#script to read and plot batch effects for sequencing centre
twins <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples_pheno.centre", sep=" ",header=F)
rm(list=ls())
alspac <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_ALSPAC_samples_pheno.centre", sep=" ",header=F)
joint <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples_pheno.centre_pop", sep=" ",header=F)
mds <- read.table("all.mds.mod",sep="\t",header=T)
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

filter <- "5e-6"
filter <- "1e-5"
filter <- "UNFILTERED"
filter <- "5e-8"

cohort <- "JOINT"
filter <- "MAF_ge_01"
cohort <- "TWINS"
cohort <- "ALSPAC"
chr <- "8"
type <- "UNFILTERED"
type <- "FILTERED"
basefolder <- "/nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/"

#Code to generalize: now only for CHR4
# library("kernlab")
library("fpc")

#sc <- kkmeans(as.matrix(merged_new[,c(4,5)]), centers=4)
#merged_new$CLUSTER_1 <- sc@.Data

#sc <- kmeans(as.matrix(merged_new[,c(4,5)]), centers=4)
sc <- pam(as.matrix(merged_new[,c(6,7)]), k=2)
sc2 <- pam(as.matrix(merged_new[,c(8,9)]), k=4)
sc <- pam(as.matrix(merged_new[,c(4,5)]), k=6)
merged_new$CLUSTER_1 <- sc$clustering
merged_new$CLUSTER_2 <- sc2$clustering

x11()
plot(merged_new[which(merged_new$CLUSTER_1 == 1),]$C1,merged_new[which(merged_new$CLUSTER_1 == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$CLUSTER_1 == 2),]$C1,merged_new[which(merged_new$CLUSTER_1 == 2),]$C2,col="red")
points(merged_new[which(merged_new$CLUSTER_1 == 3),]$C1,merged_new[which(merged_new$CLUSTER_1 == 3),]$C2,col="green")
points(merged_new[which(merged_new$CLUSTER_1 == 4),]$C1,merged_new[which(merged_new$CLUSTER_1 == 4),]$C2,col="blue")
points(merged_new[which(merged_new$CLUSTER_1 == 5),]$C1,merged_new[which(merged_new$CLUSTER_1 == 5),]$C2,col="purple")
points(merged_new[which(merged_new$CLUSTER_1 == 6),]$C1,merged_new[which(merged_new$CLUSTER_1 == 6),]$C2,col="brown")
abline(v=cluster_thr, col="green")


plot(merged_new[which(merged_new$CLUSTER_1 == 1),]$C3,merged_new[which(merged_new$CLUSTER_1 == 1),]$C4,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$CLUSTER_1 == 2),]$C3,merged_new[which(merged_new$CLUSTER_1 == 2),]$C4,col="red")
x11()
plot(merged_new[which(merged_new$CLUSTER_2 == 1),]$C5,merged_new[which(merged_new$CLUSTER_2 == 1),]$C6,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$CLUSTER_2 == 2),]$C5,merged_new[which(merged_new$CLUSTER_2 == 2),]$C6,col="red")
points(merged_new[which(merged_new$CLUSTER_2 == 3),]$C5,merged_new[which(merged_new$CLUSTER_2 == 3),]$C6,col="green")
points(merged_new[which(merged_new$CLUSTER_2 == 4),]$C5,merged_new[which(merged_new$CLUSTER_2 == 4),]$C6,col="blue")


#to work on chr8 and use only extreme clustering
merged_extreme <- merged_new[-which(merged_new$CLUSTER_1 == 2),]
merged_extreme <- merged_extreme[-which(merged_extreme$CLUSTER_1 == 4),]
merged_old <- merged_new
merged_new <- merged_extreme
dim(merged_old)
dim(merged_new)

# x11()
# plot(merged_new[which(merged_new$village_cod == 1),]$C1,merged_new[which(merged_new$village_cod == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main="PCA C1 vs C2, FVG",xlab="C1",ylab="C2")
# points(merged_new[which(merged_new$village_cod == 2),]$C1,merged_new[which(merged_new$village_cod == 2),]$C2,col="red")
# points(merged_new[which(merged_new$village_cod == 3),]$C1,merged_new[which(merged_new$village_cod == 3),]$C2,col="green")
# points(merged_new[which(merged_new$village_cod == 4),]$C1,merged_new[which(merged_new$village_cod == 4),]$C2,col="blue")
# points(merged_new[which(merged_new$village_cod == 5),]$C1,merged_new[which(merged_new$village_cod == 5),]$C2,col="yellow")
# points(merged_new[which(merged_new$village_cod == 6),]$C1,merged_new[which(merged_new$village_cod == 6),]$C2,col="purple")
# legend(0.02,-0.02,legend=c("Clauzetto","Erto","Illegio","Resia","Sauris","SMC"),fill=c("black","red","green","blue","yellow","purple"))

# twins_cluster_thr <- 0.01 
# alspac_cluster_thr <- -0.01
right_cluster <- merged_new[which(merged_new$CLUSTER_1 == 1 | merged_new$CLUSTER_1 == 6),]
left_cluster <- merged_new[which(merged_new$CLUSTER_1 == 5 | merged_new$CLUSTER_1 == 3),]
# right_cluster <- merged_new[which(merged_new$C2 > alspac_cluster_thr),]
# left_cluster <- merged_new[which(merged_new$C2 <= alspac_cluster_thr),]
# right_cluster <- merged_new[which(merged_new$C3 > twins_cluster_thr ),]
# left_cluster <- merged_new[which(merged_new$C3 <= twins_cluster_thr ),]
right_cluster$CLUSTER_2 <- 1
left_cluster$CLUSTER_2 <- 2
dim(right_cluster)
dim(left_cluster)

samples_cluster <- rbind(right_cluster[,c(1,2,17)],left_cluster[,c(1,2,17)])
# samples_cluster <- rbind(right_cluster[,c(1,2,16)],left_cluster[,c(1,2,16)])
# write.table(samples_cluster,file=paste(cohort,chr,"cluster.pheno",sep="_"),sep=" ",col.names=F,quote=F,row.names=F)
# write.table(samples_cluster,file=paste(cohort,chr,"cluster_extreme.pheno",sep="_"),sep=" ",col.names=F,quote=F,row.names=F)
write.table(samples_cluster,file=paste(cohort,chr,"cluster_middle_vs_low.pheno",sep="_"),sep=" ",col.names=F,quote=F,row.names=F)

#check if the clustering is ok
x11()
plot(right_cluster[which(right_cluster$CLUSTER_2 == 1),]$C1,right_cluster[which(right_cluster$CLUSTER_2 == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(left_cluster[which(left_cluster$CLUSTER_2 == 2),]$C1,left_cluster[which(left_cluster$CLUSTER_2 == 2),]$C2,col="red")

x11()
plot(right_cluster[which(right_cluster$CLUSTER_2 == 1),]$C1,right_cluster[which(right_cluster$CLUSTER_2 == 1),]$C3,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C3)-0.001,max(merged_new$C3)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(left_cluster[which(left_cluster$CLUSTER_2 == 2),]$C1,left_cluster[which(left_cluster$CLUSTER_2 == 2),]$C3,col="red")
abline(v=cluster_thr, col="green")

#for CHR6
sc <- pam(as.matrix(merged_new[,c(4,5)]), k=4)
x11()
plot(merged_new[which(merged_new$CLUSTER_1 == 1),]$C1,merged_new[which(merged_new$CLUSTER_1 == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$CLUSTER_1 == 2),]$C1,merged_new[which(merged_new$CLUSTER_1 == 2),]$C2,col="red")
points(merged_new[which(merged_new$CLUSTER_1 == 3),]$C1,merged_new[which(merged_new$CLUSTER_1 == 3),]$C2,col="green")
points(merged_new[which(merged_new$CLUSTER_1 == 4),]$C1,merged_new[which(merged_new$CLUSTER_1 == 4),]$C2,col="blue")


samples_cluster <- rbind(right_cluster[,c(1,2,17)],left_cluster[,c(1,2,17)])
write.table(merged_new[,c(1,2,17)],file=paste(cohort,chr,"cluster.pheno",sep="_"),sep=" ",col.names=F,quote=F,row.names=F)
write.table(merged_new[,c(1,2,16)],file=paste(cohort,chr,"cluster.pheno",sep="_"),sep=" ",col.names=F,quote=F,row.names=F)

#then we perform a logistic analysis on cluster as phenotype and center as covariate

# plink --noweb \
# --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr4.pruned.bed \
# --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr4.pruned.bim.COPY \
# --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr4.pruned.fam \
# --allow-no-sex \
# --maf 0.01 \
# --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples.keeplist \
# --pheno /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/UNFILTERED_1/CHR4/TWINS_4_cluster.pheno \
# --covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples_pheno.centre \
# --logistic \
# --adjust \
# --out chr4_cluster_center.txt

#now plot association results
source("/nfs/users/nfs_m/mc14/Work/r_scripts/qqman.R")
chr <- 12
cohort <- "ALSPAC"
cohort <- "TWINS"

# to_plot_unfiltered <- read.table("chr12_cluster_c3_center.txt.assoc.logistic", header=T)
to_plot_unfiltered <- read.table(paste("chr",chr,"_cluster_center.txt.assoc.logistic",sep=''), header=T)

to_plot_unfiltered <- to_plot_unfiltered[which(to_plot_unfiltered$TEST == "ADD"),]
# to_plot_unfiltered = read.csv("chr6.txt.assoc.logistic.mod", sep="\t", header=T)

to_plot_mafs <- read.table("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/joint_NEW/MDS/GWA/gwa_data_freq.frq.mod",header=T,sep="\t")

# to add maf info
merged_new <- merge(to_plot_unfiltered,to_plot_mafs,by.x="SNP",by.y="SNP",all.x)
merged_new$SNP <- as.character(merged_new$SNP)
merged_new$CHR.y <- NULL
merged_new$A1.y <- NULL
names(merged_new)[2] <- "CHR"
names(merged_new)[4] <- "A1"

dim(to_plot_unfiltered)
dim(merged_new)

#colnames(to_plot_unfiltered)=c("SNP","CHR","BP","P")
# colnames(to_plot_unfiltered)[1] <- "CHR"
# colnames(to_plot_unfiltered)[2] <- "SNP"
# colnames(to_plot_unfiltered)[3] <- "BP"
# colnames(to_plot_unfiltered)[9] <- "P"

dim(to_plot_unfiltered[which(to_plot_unfiltered$P < 1 & to_plot_unfiltered$P >= 0.1),])
dim(to_plot_unfiltered[which(to_plot_unfiltered$P < 0.1 & to_plot_unfiltered$P >= 0.01),])
dim(to_plot_unfiltered[which(to_plot_unfiltered$P < 0.01),])


merged_new_reg <- merged_new[which(merged_new$BP >= 1295197 & merged_new$BP <= 1820529),]
low_to_comm_reg <- merged_new_reg[which(merged_new_reg$MAF >= 0.01 & merged_new_reg$MAF < 0.05),]$SNP
manhattan(merged_new_reg, annotate=low_to_comm_reg, pch=16, main=paste("chr",chr," unfiltered",sep=""))
write.table(merged_new_reg,file=paste(cohort,chr,"detail.region",sep="_"),sep=" ",col.names=T,quote=F,row.names=F)

#extract snps to color
low_to_comm <- merged_new[which(merged_new$MAF >= 0.01 & merged_new$MAF < 0.05),]$SNP
# jpeg(paste(args[[1]], "_", args[[2]], ".manhattan.unfiltered.jpg", sep=""), width = 1300, height = 600, pointsize = 16)
jpeg(paste(cohort,"_test_logistic_chr",chr,".manhattan.unfiltered.jpg",sep=""), width = 1300, height = 600, pointsize = 16)
# manhattan(merged_new, annotate=low_to_comm, pch=16, main=paste(args[[1]], "_", args[[2]], " unfiltered", sep=""))
manhattan(merged_new, annotate=low_to_comm, pch=16, main=paste("chr",chr," unfiltered",sep=""))
# jpeg(paste("TC", "_", "22", ".manhattan.unfiltered.jpg", sep=""), width = 1300, height = 600, pointsize = 16)
# manhattan(to_plot_unfiltered, pch=16, main=paste("TC", "_", "22", " unfiltered", sep=""))
dev.off()

# jpeg(paste(args[[1]], "_", args[[2]], ".qq.unfiltered.jpg", sep=""), width = 1300, height = 600, pointsize = 16)
jpeg(paste(cohort,"_test_logistic_chr",chr,".qq.unfiltered.jpg",sep=""), width = 1300, height = 600, pointsize = 16)
# qq(to_plot_unfiltered$P, main=paste(args[[1]], "_", args[[2]], " unfiltered", sep=""))
qq(to_plot_unfiltered$P, main=paste("chr",chr," unfiltered",sep=""))
# jpeg(paste("TC", "_", "22", ".qq.unfiltered.jpg", sep=""), width = 1300, height = 600, pointsize = 16)
# qq(to_plot_unfiltered$P, pch=16, main=paste("TC", "_", "22", " unfiltered", sep=""))
dev.off()

###########SINGLE CHR###############
xmin <- min(merged_new$C1)-0.001
ymin <- min(merged_new$C2)+0.001

pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_",chr,".pdf",sep=""),width=11.7, height=8.3)
par(lab=c(4,4,6),mfrow=c(2,2))
plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C2,xlim=c(min(merged_new$C1)-0.001,max(merged_new$C1)+0.001),ylim=c(min(merged_new$C2)-0.001,max(merged_new$C2)+0.001),main=paste("PCA C1 vs C2, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C2,col="red")
legend(-0.0035,0.0045,legend=c("BGI","SC"),fill=c("black","red"))

plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C3,xlim=c(min(merged_new$C1),max(merged_new$C1)),ylim=c(min(merged_new$C3),max(merged_new$C3)),main=paste("PCA C1 vs C3, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C3,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C1,merged_new[which(merged_new$SEQ == 1),]$C4,xlim=c(min(merged_new$C1),max(merged_new$C1)),ylim=c(min(merged_new$C4),max(merged_new$C4)),main=paste("PCA C1 vs C4, pruned ",cohort,"  - ",type,sep=""),xlab="C1",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2),]$C1,merged_new[which(merged_new$SEQ == 2),]$C4,col="red")

plot(merged_new[which(merged_new$SEQ == 1),]$C2,merged_new[which(merged_new$SEQ == 1),]$C3,xlim=c(min(merged_new$C2),max(merged_new$C2)),ylim=c(min(merged_new$C3),max(merged_new$C3)),main=paste("PCA C2 vs C3, pruned ",cohort,"  - ",type,sep=""),xlab="C2",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2),]$C2,merged_new[which(merged_new$SEQ == 2),]$C3,col="red")
dev.off()

########### GWA ###############
#legend position:
xmin <- min(merged_new$C1)-0.001
ymin <- min(merged_new$C2)+0.001

pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_gwa.pdf",sep=""),width=11.7, height=8.3)
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

#SEQ=1 -> BGI
#SEQ=2 -> SANGER
#POP=1 -> ALSPAC
#POP=2 -> TWINS
########### GWA JOINT 3621 samples ###############
pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_gwa.pdf",sep=""),width=11.7, height=8.3)
par(lab=c(2,2,4),mfrow=c(2,2))
plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,xlim=c(min(merged_new$C1)-0.005,max(merged_new$C1)+0.005),ylim=c(min(merged_new$C2)-0.005,max(merged_new$C2)+0.005),main=paste("PCA C1 vs C2, pruned ",cohort," - ",filter,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,col="green")
legend(-0.015,0.01,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,xlim=c(min(merged_new$C1)-0.003,max(merged_new$C1)+0.003),ylim=c(min(merged_new$C3)-0.003,max(merged_new$C3)+0.003),main=paste("PCA C1 vs C3, pruned ",cohort," - ",filter,sep=""),xlab="C1",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,col="green")
# legend(0.003,0.006,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C4,xlim=c(min(merged_new$C1)-0.003,max(merged_new$C1)+0.003),ylim=c(min(merged_new$C4)-0.003,max(merged_new$C4)+0.003),main=paste("PCA C1 vs C4, pruned ",cohort," - ",filter,sep=""),xlab="C1",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C4,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C4,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C4,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,xlim=c(min(merged_new$C2)-0.003,max(merged_new$C2)+0.003),ylim=c(min(merged_new$C3)-0.003,max(merged_new$C3)+0.003),main=paste("PCA C2 vs C3, pruned ",cohort," - ",filter,sep=""),xlab="C2",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,col="green")
# legend(0.003,0.006,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))
dev.off()

############################ JOINT by CHR
xmin <- min(merged_new$C1)-0.001
ymin <- min(merged_new$C2)+0.001

pdf(paste(basefolder,"PCA_",cohort,"_batch_PRUNED_",type,"_",filter,"_",chr,".pdf",sep=""),width=11.7, height=8.3)
par(lab=c(2,2,4),mfrow=c(2,2))
plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,xlim=c(min(merged_new$C1)-0.005,max(merged_new$C1)+0.005),ylim=c(min(merged_new$C2)-0.005,max(merged_new$C2)+0.005),main=paste("PCA C1 vs C2, pruned ",cohort," - ",filter,sep=""),xlab="C1",ylab="C2")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,col="green")
legend(xmin,ymin,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,xlim=c(min(merged_new$C1)-0.003,max(merged_new$C1)+0.003),ylim=c(min(merged_new$C3)-0.003,max(merged_new$C3)+0.003),main=paste("PCA C1 vs C3, pruned ",cohort," - ",filter,sep=""),xlab="C1",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,col="green")
# legend(0.003,0.006,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C4,xlim=c(min(merged_new$C1)-0.003,max(merged_new$C1)+0.003),ylim=c(min(merged_new$C4)-0.003,max(merged_new$C4)+0.003),main=paste("PCA C1 vs C4, pruned ",cohort," - ",filter,sep=""),xlab="C1",ylab="C4")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C4,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C4,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C1,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C4,col="green")

plot(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C2,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 1),]$C3,xlim=c(min(merged_new$C2)-0.003,max(merged_new$C2)+0.003),ylim=c(min(merged_new$C3)-0.003,max(merged_new$C3)+0.003),main=paste("PCA C2 vs C3, pruned ",cohort," - ",filter,sep=""),xlab="C2",ylab="C3")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C2,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 1),]$C3,col="blue")
points(merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C2,merged_new[which(merged_new$SEQ == 1 & merged_new$POP == 2),]$C3,col="red")
points(merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C2,merged_new[which(merged_new$SEQ == 2 & merged_new$POP == 2),]$C3,col="green")
# legend(0.003,0.006,legend=c("ALSPAC BGI","ALSPAC SC","TWINS BGI","TWINS SC"),fill=c("black","blue","red","green"))
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
