#script to plot all chr files together
rm(list=ls())
require(gplots)

#Set arguments to use from command line
# 
#read commandline args
args <- commandArgs(trailing=TRUE)

#set population files names and paths

pop_1_name <- args[[1]]
pop_2_name <- args[[2]]
pop_overlap_filepath <- args[[3]]

pop_overlap <- read.table(pop_overlap_filepath, sep='\t',header=F,comment.char="",colClasses=c("character","numeric","character","character","numeric","numeric","numeric","factor"))

colnames(pop_overlap) <- c("CHROM","POS",paste(pop_1_name,"ID",sep="_"),paste(pop_2_name,"ID",sep="_"),paste(pop_1_name,"MAF",sep="_"),paste(pop_2_name,"MAF",sep="_"),"D_MAF","SIGN")

#split increased maf and decreased maf
increase_diff <- pop_overlap[which(pop_overlap$D_MAF > 0 & as.character(pop_overlap$SIGN) == "+"),]
decrease_diff <- pop_overlap[which(pop_overlap$D_MAF > 0 & as.character(pop_overlap$SIGN) == "-"),]
same_diff <- pop_overlap[which(as.character(pop_overlap$SIGN) == "="),]


#write tables
write.table(increase_diff,file=paste(pop_1_name,"vs",pop_2_name,"increase.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)
write.table(decrease_diff,file=paste(pop_1_name,"vs",pop_2_name,"decrease.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)
write.table(same_diff,file=paste(pop_1_name,"vs",pop_2_name,"same.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)

#write summaries
write.table(summary(increase_diff),file=paste("summary",pop_1_name,"vs",pop_2_name,"increase.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)
write.table(summary(decrease_diff),file=paste("summary",pop_1_name,"vs",pop_2_name,"decrease.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)
write.table(summary(same_diff),file=paste("summary",pop_1_name,"vs",pop_2_name,"same.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)
write.table(summary(pop_overlap),file=paste("summary",pop_1_name,"vs",pop_2_name,"all_overlap.txt",sep="_"),sep="\t",col.names=T,quote=F,row.names=F)

#plot density for maf and dmaf for all pop
plot_name <- paste(pop_1_name,pop_2_name,"DENSITY_MAF_DELTA_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
	par(lab=c(15,15,17),cex=2)
	plot(density(pop_overlap$D_MAF), main=paste("Distribution of Delta MAF",pop_1_name,"vs",pop_2_name,sep=" "),xlab="Delta MAF",ylab="Delta MAF density")
	lines(density(pop_overlap[,5]),col=c("red"))
	lines(density(pop_overlap[,6]),col=c("blue"))
	smartlegend(x="right",y="center", inset = 0,c("DMAF", "VBI MAF","EUR MAF"),fill = c("black", "red","blue"))
dev.off()

#plot density for maf and dmaf for all pop
plot_name <- paste(pop_1_name,pop_2_name,"DENSITY_DELTA_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
	par(lab=c(15,15,17),cex=2)
	plot(density(pop_overlap$D_MAF), main=paste("Distribution of Delta MAF",pop_1_name,"vs",pop_2_name,sep=" "),xlab="Delta MAF",ylab="Delta MAF density")
dev.off()

#plot maf density for each population
plot_name <- paste(pop_1_name,"DENSITY_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
	par(lab=c(15,15,17),cex=2)
	plot(density(pop_overlap[,5]), main=paste("Distribution of MAF",pop_1_name,sep=" "),xlab="MAF",ylab="MAF density distribution")
dev.off()

plot_name <- paste(pop_2_name,"DENSITY_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
	par(lab=c(15,15,17),cex=2)
	plot(density(pop_overlap[,6]), main=paste("Distribution of MAF",pop_2_name,sep=" "),xlab="MAF",ylab="MAF density distribution")
dev.off()

#now make a boxplot different populations
box_plot_name <- paste(pop_1_name,pop_2_name,"BOX_MAF.jpg",sep="_")
jpeg(box_plot_name, width=1000, height=1000)
	par(cex=2)
	boxplot(pop_overlap[,5],pop_overlap[,6],names=c("VBI","EUR"),main="Minor allele frequencies", xlab="Populations",ylab="MAF")
dev.off()

#now make a boxplot different populations
#increase:
box_plot_name <- paste(pop_1_name,pop_2_name,"BOX_INCREASE_MAF.jpg",sep="_")
jpeg(box_plot_name, width=1000, height=1000)
	par(cex=2)
	boxplot(increase_diff[,5],increase_diff[,6],names=c("VBI","EUR"),main="Minor allele frequencies", xlab="Populations",ylab="MAF")
dev.off()

#decrease:
box_plot_name <- paste(pop_1_name,pop_2_name,"BOX_DECREASE_MAF.jpg",sep="_")
jpeg(box_plot_name, width=1000, height=1000)
	par(cex=2)
	boxplot(decrease_diff[,5],decrease_diff[,6],names=c("VBI","EUR"),main="Minor allele frequencies", xlab="Populations",ylab="MAF")
dev.off()

#same freq
box_plot_name <- paste(pop_1_name,pop_2_name,"BOX_SAME_MAF.jpg",sep="_")
jpeg(box_plot_name, width=1000, height=1000)
	par(cex=2)
	boxplot(same_diff[,5],same_diff[,6],names=c("VBI","EUR"),main="Minor allele frequencies", xlab="Populations",ylab="MAF")
dev.off()

#Now plot density of dmaf for increase and decrease, and maf for =
increase_plot_name <- paste(pop_1_name,pop_2_name,"INCREASE_DENSITY_DMAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
        par(lab=c(15,15,17),cex=2)
        plot(density(increase_diff$D_MAF), main=paste("Distribution of Delta MAF",pop_1_name,"vs",pop_2_name,sep=" "),xlab="Delta MAF",ylab="Delta MAF density distribution")
dev.off()

decrease_plot_name <- paste(pop_1_name,pop_2_name,"DECREASE_DENSITY_DMAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
        par(lab=c(15,15,17),cex=2)
        plot(density(decrease_diff$D_MAF), main=paste("Distribution of Delta MAF",pop_1_name,"vs",pop_2_name,sep=" "),xlab="Delta MAF",ylab="Delta MAF density distribution")
dev.off()

same_plot_name <- paste(pop_1_name,pop_2_name,"SAME_DENSITY_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
        par(lab=c(15,15,17),cex=2)
        plot(density(same_diff[,5]), main=paste("Distribution of MAF",pop_1_name,"vs",pop_2_name,sep=" "),xlab="MAF",ylab="MAF density distribution")
dev.off()

#plot increase and decrease (for both pop) and same maf: density (and manhattan plot?)
#increase:
plot_name <- paste(pop_1_name,pop_2_name,"DENSITY_INCREASED_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
        par(lab=c(15,15,17),cex=2,lwd=2,xaxp=c(min(c(min(increase_diff[,5]),min(increase_diff[,6]))),max(c(max(increase_diff[,5]),max(increase_diff[,6]))),1),yaxp=c(min(c(min(density(increase_diff[,5])$x),min(density(increase_diff[,6])$x))),max(c(max(density(increase_diff[,5])$x),max(density(increase_diff[,6])$x))),1)
        plot(density(increase_diff[,5]), main=paste("Distribution of Increased MAF in",pop_1_name,"vs",pop_2_name,sep=" "),xlab="MAF",ylab="MAF density", col=c("red"))
        lines(density(increase_diff[,6]),col=c("blue"))
        smartlegend(x="right",y="center", inset = 0,c("VBI MAF","EUR MAF"),fill = c("red","blue"))
dev.off()

#decrease:
plot_name <- paste(pop_1_name,pop_2_name,"CHR",chr,"DENSITY_decreased_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
        par(lab=c(15,15,17),cex=2,lwd=2,xaxp=c(min(c(min(decrease_diff[,5]),min(decrease_diff[,6]))),max(c(max(decrease_diff[,5]),max(decrease_diff[,6]))),1),yaxp=c(min(c(min(density(decrease_diff[,5])$x),min(density(decrease_diff[,6])$x))),max(c(max(density(decrease_diff[,5])$x),max(density(decrease_diff[,6])$x))),1)
        plot(density(decrease_diff[,5]), main=paste("Distribution of Decreased MAF in",pop_1_name,"vs",pop_2_name,sep=" "),xlab="MAF",ylab="MAF density", col=c("red"))
        lines(density(decrease_diff[,6]),col=c("blue"))
        smartlegend(x="right",y="center", inset = 0,c("VBI MAF","EUR MAF"),fill = c("red","blue"))
dev.off()

#same:
plot_name <- paste(pop_1_name,pop_2_name,"CHR",chr,"DENSITY_same_MAF.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
        par(lab=c(15,15,17),cex=2,lwd=2,xaxp=c(min(c(min(same_diff[,5]),min(same_diff[,6]))),max(c(max(same_diff[,5]),max(same_diff[,6]))),1),yaxp=c(min(c(min(density(same_diff[,5])$x),min(density(same_diff[,6])$x))),max(c(max(density(same_diff[,5])$x),max(density(same_diff[,6])$x))),1)
        plot(density(same_diff[,5]), main=paste("Distribution of constant MAF in",pop_1_name,"vs",pop_2_name,sep=" "),xlab="MAF",ylab="MAF density", col=c("red"))
        lines(density(same_diff[,6]),col=c("blue"))
        smartlegend(x="right",y="center", inset = 0,c("VBI MAF","EUR MAF"),fill = c("red","blue"))
dev.off()


#make a test to know if differences are significatives
#for the whole set
wilcox_test <- wilcox.test(pop_overlap[,5],pop_overlap[,6], alternative="greater")
#write(wilcox_test,file=paste(pop_1_name,"vs",pop_2_name,"wilcox.txt",sep="_"),sep="\t")
sink(file=paste(pop_1_name,"vs",pop_2_name,"wilcox.txt",sep="_"))
print(wilcox_test)
sink()

#only for increase
increase_wilcox_test <- wilcox.test(increase_diff[,5],increase_diff[,6], alternative="greater")
#write(increase_wilcox_test,file=paste(pop_1_name,"vs",pop_2_name,"increase_wilcox.txt",sep="_"),sep="\t")
sink(file=paste(pop_1_name,"vs",pop_2_name,"increase_wilcox.txt",sep="_"))
print(increase_wilcox_test)
sink()

#only for decrease
decrease_wilcox_test <- wilcox.test(decrease_diff[,5],decrease_diff[,6], alternative="less")
#write(decrease_wilcox_test,file=paste(pop_1_name,"vs",pop_2_name,"decrease_wilcox.txt",sep="_"),sep="\t")
sink(file=paste(pop_1_name,"vs",pop_2_name,"decrease_wilcox.txt",sep="_"))
print(decrease_wilcox_test)
sink()


#Use a for cicle to plot and extract info for maf classes:
classes <- c("increase","decrease","same")
#now split all maf in frequencies bins
#this is a function that splits in bins of maf and writes some summaries
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
populations <- c(pop_1_name,pop_2_name)

for ( i in 1:length(classes)){
  current_class <- classes[i]
  current_snp_set <- get(paste(current_class,"diff",sep='_'))

  #this is a function that splits in bins of maf and writes some summaries
  assign(paste(pop_1_name,current_class,"maf_classes",sep="_"), split_bins(maf_classes,as.data.frame(cbind(CHROM=current_snp_set$CHROM,MAF=current_snp_set[,5])),pop_1_name))
  assign(paste(pop_2_name,current_class,"maf_classes",sep="_"), split_bins(maf_classes,as.data.frame(cbind(CHROM=current_snp_set$CHROM,MAF=current_snp_set[,6])),pop_2_name))

  #write a cute output
  for(i in 1:length(populations)){
    sink(paste(populations[i],"CHR",chr,current_class,"maf_classes_bin_resume.txt",sep="_"))
    print(get(paste(populations[i],current_class,"maf_classes",sep="_")))
    sink()  
  }
  
}


q(save="no")
