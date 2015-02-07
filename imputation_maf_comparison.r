#Script to compare maf before and after imputation 
require(gplots)

#clean workspace
rm(list=ls())

#open low freq table
table_path <- '/lustre/scratch110/sanger/mc14/GENOTIPI/SEQUENCED/SEQUENCES/new_rel/RSIDFIX/NEW_ANNOTATED_RELEASE/CHR_SPLITTED/MAF/MAF_stats_7BINS/VBI_seq_maf_lte_0.02_table.txt'
impute_path <- '/lustre/scratch110/sanger/mc14/GENOTIPI/IMP_VBI_1KGP/LOW_FREQ_COMPARISON/MAF_lte_002_tot/'

low_freq <- read.table(table_path,sep="\t", header=T, stringsAsFactors=F)

#create a dataframe to collect all chr
all_chr_merged <- NULL

#now for each chr, we need to upload the file obtained before
for (i in 1:22){
	current_set_name <- paste("chr",i,"_imputed",sep="")
	current_chr <- read.table(paste(impute_path,'chr',i,'_low_freq.txt',sep=""),sep=" ",head=F,stringsAsFactors=F)
	colnames(current_chr) <- c('CHROM','ID','POS','exp_freq_a1','info','certainty','type','info_type1','concord_type1','r2_type1','info_type0','concord_type0','r2_type0','MAF')
	current_chr$info_type1 <- NULL
	current_chr$info_type0 <- NULL
	current_chr$concord_type1 <- NULL
	current_chr$concord_type0 <- NULL
	current_chr$r2_type1<- NULL
	current_chr$r2_type0<- NULL
	current_chr_merged <- merge(current_chr,low_freq[which(low_freq$CHROM == i),],by="POS",all.x=F)
	colnames(current_chr_merged) <- c("POS","CHROM","ID","exp_freq_a1","info","certainty","type","MAF","CHROM_bef","ID_bef","AC","AN","AF","MAF_bef")
	current_chr_merged$exp_maf_1664 <- ((1664*current_chr_merged$MAF_bef)/(110))
	
	#plots and summaries
	jpeg(paste("imp_maf_density",current_set_name,"jpg",sep="."),width=1000,height=1000)
	par(lab=c(15,15,17),cex=2)
		plot(density(current_chr_merged$info),main=paste('Chromosome',i,'info score density',sep=' '))
	dev.off()

	jpeg(paste("qqplot",current_set_name,"jpg",sep="."),width=1000,height=1000)
	par(lab=c(15,15,17),cex=2)
		qqplot(current_chr_merged$MAF,current_chr_merged$info, main=paste('Chromosome',i,'MAF vs info score',sep=' '))
	dev.off()

	jpeg(paste("maf_bef",current_set_name,"jpg",sep="."),width=1000,height=1000)
	par(lab=c(15,15,17),cex=2)
		plot(current_chr_merged$POS,current_chr_merged$MAF_bef, main=paste('Chromosome',i,'POS vs MAF before imputation',sep=' '),xlab="POS",ylab="VBI seq MAF")
	dev.off()

	jpeg(paste("maf_imp",current_set_name,"jpg",sep="."),width=1000,height=1000)
	par(lab=c(15,15,17),cex=2)
		plot(current_chr_merged$POS,current_chr_merged$MAF,main=paste('Chromosome',i,'POS vs MAF after imputation',sep=' '),xlab="POS",ylab="VBI MAF post imputation")
	dev.off()
	
	all_chr_merged <- rbind(all_chr_merged,current_chr_merged)
	rm(current_chr_merged)
	gc()
}

#now split all maf in frequencies bins
#this is a function that splits in bins of maf and writes some summaries
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
#maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
maf_classes <- c(0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)

pop_name <- "002_MAF_class"

all_maf_classes <- split_bins(maf_classes,all_chr_merged,pop_name)
gc()
#write a cute output
sink('maf_bin_resume.txt')
print(all_maf_classes)
sink()

#print for all chr, maf density before and after imputation
jpeg(paste("all_mb_mi_density","jpg",sep="."),width=1000,height=1000)
	par(lab=c(15,15,17),cex=2)
	plot(density(all_chr_merged$MAF),main='Overall maf density distribution',ylim=c(0,max(density(all_chr_merged$MAF_bef)$y)))
	lines(density(all_chr_merged$MAF_bef),col=c('red'))
	smartlegend(x="right",y="center", inset = 0,c("MAF of imputed sites", "MAF from seq sites"),fill = c("black", "red"))
dev.off()

for (i in 2:length(maf_classes)){ 
        current_class_path <- paste(pop_name,'maf_lte',maf_classes[i],'table.txt',sep="_") 
        current_class_table <- read.table(current_class_path,sep="\t",header=T)
        
        jpeg(paste(current_class_path,"info_density.jpg",sep="_"),width=1000,height=1000)
                par(lab=c(15,15,17),cex=2)
                plot(density(current_class_table$info),main='Info score density distribution')
        dev.off()
        jpeg(paste(current_class_path,"info_MAF.jpg",sep="_"),width=1000,height=1000)
                par(lab=c(15,15,17),cex=2)
                qqplot(current_class_table$MAF,current_class_table$info, main='MAF vs info score', xlab='MAF',ylab='Impute info score')
        dev.off()
}


