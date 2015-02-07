#script to plot maf and dmaf for general files 
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
current_class <- args[[4]]

pop_overlap <- read.table(pop_overlap_filepath, sep='\t',header=F,comment.char="",colClasses=c("character","numeric","character","character","numeric","numeric","numeric","factor"))

colnames(pop_overlap) <- c("CHROM","POS",paste(pop_1_name,"ID",sep="_"),paste(pop_2_name,"ID",sep="_"),paste(pop_1_name,"MAF",sep="_"),paste(pop_2_name,"MAF",sep="_"),"D_MAF","SIGN")

#write summaries
write.table(summary(pop_overlap),file=paste("summary_",pop_overlap_filepath,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#plot density for dmaf for all pop
plot_name <- paste(pop_overlap_filepath,current_class,"DELTA_MAF_DENSITY.jpg",sep="_")
jpeg(plot_name, width=1000, height=1000)
        par(lab=c(15,15,17),cex=2)
        plot(density(pop_overlap$D_MAF), main=paste("Distribution of Delta MAF",pop_1_name,"vs",pop_2_name,sep=" "),xlab="Delta MAF",ylab="Delta MAF density")
dev.off()

#plot density for maf and dmaf for all pop
plot_name <- paste(pop_overlap_filepath,"DENSITY_MAF.jpg",sep="_")
minxlim <- min(c(min(density(pop_overlap[,5])$x),min(density(pop_overlap[,6])$x)))
maxxlim <- max(c(max(density(pop_overlap[,5])$x),max(density(pop_overlap[,6])$x)))
minylim <- min(c(min(density(pop_overlap[,5])$y),min(density(pop_overlap[,6])$y)))
maxylim <- max(c(max(density(pop_overlap[,5])$y),max(density(pop_overlap[,6])$y)))

jpeg(plot_name, width=1000, height=1000)
  par(lab=c(15,15,17),cex=2,lwd=2)
  plot(density(pop_overlap[,6]), main=paste("Distribution of MAF",pop_1_name,"vs",pop_2_name,sep=" "),xlab="MAF",ylab="MAF density", xlim=c(minxlim,maxxlim), ylim=c(minylim,maxylim),col=c("blue"))
  lines(density(pop_overlap[,5]),col=c("red"))
  #lines(density(pop_overlap[,6]),col=c("blue"))
  smartlegend(x="right",y="center", inset = 0,c(paste(pop_1_name,"MAF",sep=" "),paste(pop_2_name,"MAF",sep=" ")),fill = c("red","blue"))
dev.off()

#now split all maf in frequencies bins
#this is a function that splits in bins of maf and writes some summaries
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
populations <- c(pop_1_name,pop_2_name)

#this is a function that splits in bins of maf and writes some summaries
assign(paste(pop_1_name,current_class,"maf_classes",sep="_"), split_bins(maf_classes,as.data.frame(cbind(CHROM=pop_overlap$CHROM,MAF=pop_overlap[,5],D_MAF=pop_overlap$D_MAF)),pop_1_name))
assign(paste(pop_2_name,current_class,"maf_classes",sep="_"), split_bins(maf_classes,as.data.frame(cbind(CHROM=pop_overlap$CHROM,MAF=pop_overlap[,6],D_MAF=pop_overlap$D_MAF)),pop_2_name))

#write a cute output
for(i in 1:length(populations)){
  sink(paste(populations[i],pop_overlap_filepath,"maf_classes_bin_resume.txt",sep="_"))
  print(get(paste(populations[i],current_class,"maf_classes",sep="_")))
  sink()  


	maf_and_dmaf_plot_name <- paste(populations[i],"ALL_MAF.jpg",sep="_")

	#we need to plot all dmaf in the same figure
	jpeg(maf_and_dmaf_plot_name, width=2000, height=1000)
	par(lab=c(15,15,17),cex=2,lwd=2)
	minxlim <- 0
	maxxlim <- 0.1
        minylim <- 0
        maxylim <- 800

	#Now plot for each maf cass and for each population,maf and dmaf, and maf vs dmaf
	legend_string <- NULL
	legend_colors <- NULL

	for(j in 2:length(maf_classes)){
		current_maf_class <- maf_classes[j]
		
		#read file:
		filename <- paste(populations[i],"maf_lte",current_maf_class,"table.txt",sep="_")
		current_set <- read.table(filename,sep='\t',header=T,comment.char="",colClasses=c("character","numeric","numeric"))
		colnames(current_set) <- c("CHROM","MAF","D_MAF")

		if( j ==2 ){
			plot(density(current_set$D_MAF), main=paste("Distribution of Delta MAF for",populations[i],sep=" "),xlab="Delta MAF",ylab="Density distribution", xlim=c(minxlim,maxxlim), ylim=c(minylim,maxylim),col=c(colors()[j*11]))
		}else{
			lines(density(current_set$D_MAF),col=c(colors()[j*11]))
		}
	
	legend <- paste("MAF bin",current_maf_class, sep=" ")
	legend_string <- c(legend_string,legend)
	legend_colors <- c(legend_colors,colors()[j*11])
	
	}
	
	smartlegend(x="right",y="center", inset = 0,legend_string,fill = legend_colors)
	dev.off()


#Now plot for each maf cass and for each population,maf and dmaf, and maf vs dmaf
	for(j in 2:length(maf_classes)){
		current_maf_class <- maf_classes[j]
		maf_and_dmaf_plot_name <- paste(populations[i],"maf_lte",current_maf_class,"MAF_AND_DMAF.jpg",sep="_")
		#read file:
		filename <- paste(populations[i],"maf_lte",current_maf_class,"table.txt",sep="_")
		current_set <- read.table(filename,sep='\t',header=T,comment.char="",colClasses=c("character","numeric","numeric"))
		colnames(current_set) <- c("CHROM","MAF","D_MAF")
		if (populations[i] == "VBI"){
			plot_col = c("red")
		}
		if (populations[i] == "EUR"){
                        plot_col = c("blue")
                }

		#plot only delta maf
		dmaf_plot_name <-paste(populations[i],"maf_lte",current_maf_class,"DMAF.jpg",sep="_")
		jpeg(dmaf_plot_name, width=1000, height=1000)
			par(lab=c(15,15,17),cex=2,lwd=2)
			plot(density(current_set$D_MAF), main=paste("Distribution of Delta MAF for",populations[i],"in",current_maf_class,sep=" "),xlab="Delta MAF",ylab="Density distribution",col=plot_col)
#		polygon(density(current_set$D_MAF),col=colors()[j*11],border=plot_col)
		dev.off()

		#plot maf vs delta maf
		maf_vs_dmaf_plot_name <-paste(populations[i],"maf_lte",current_maf_class,"MAF_vs_DMAF.jpg",sep="_")
		jpeg(maf_vs_dmaf_plot_name, width=1000, height=1000)
			par(lab=c(15,15,17),cex=2,lwd=2,pch=20)
			plot(current_set$MAF,current_set$D_MAF, main=paste("MAF vs Delta MAF for",populations[i],"in",current_maf_class,sep=" "),xlab="MAF",ylab="Delta MAF",col=plot_col)
		dev.off()

	}
}

#now make a boxplot different populations
box_plot_name <- paste(pop_overlap_filepath,"BOX_MAF.jpg",sep="_")
jpeg(box_plot_name, width=1000, height=1000)
  par(cex=2)
  boxplot(pop_overlap[,5],pop_overlap[,6],names=c(pop_1_name,pop_2_name),main="Minor allele frequencies", xlab="Populations",ylab="MAF")
dev.off()

#make a test to know if differences are significatives
#for the whole set
wilcox_test <- wilcox.test(pop_overlap[,5],pop_overlap[,6], alternative="greater")
#write(wilcox_test,file=paste(pop_1_name,"vs",pop_2_name,"wilcox.txt",sep="_"),sep="\t")
sink(file=paste(pop_overlap_filepath,"wilcox.txt",sep="_"))
print(wilcox_test)
sink()


q(save="no")
