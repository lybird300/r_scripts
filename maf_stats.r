#script for statistics on maf...
###### REPLOT with ggplot
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
require(plotrix)
require(ggplot2)
require(reshape2)
base_folder <- getwd()

#Set arguments to use from command line
# 
#read commandline args
args <- commandArgs(trailing=TRUE)

filepath <- args[[1]]
pop_name <- args[[2]]

#fixed parameters for pops
#FVG
# filepath <- "/lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140410_ANNOTATE/MAF/ALL.maf_table.snp.tab"
# pop_name <- "FVG"

#VBI
# filepath <- "/lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140410_ANNOTATE/MAF/ALL.maf_table.snp.tab"
# pop_name <- "VBI"
#chr <- args[[2]]

pops <- c("FVG")
pops <- c("VBI","FVG","CARL")
all_pop_MAF_count <- NULL
all_pop_MAF_freq <- NULL
for (pop in pops) {
	if( pop == "CARL") {
	# filepath <- "/lustre/scratch113/projects/carl_seq/variant_refinement/13102015_RELEASE/ALL_CARL_20151013_freq_only.frq"
	filepath <- "/lustre/scratch113/projects/carl_seq/variant_refinement/13102015_RELEASE/ALL_CARL_20151013_freq_only_nomono.frq"
	}
	if( pop == "VBI") {
	# filepath <- "/lustre/scratch113/projects/esgi-vbseq/08092015/13102015_RELEASE/ALL_VBI_20151013_freq_only.frq"
	filepath <- "/lustre/scratch113/projects/esgi-vbseq/08092015/13102015_RELEASE/ALL_VBI_20151013_freq_only_nomono.frq"
	}
	if( pop == "FVG") {
	# filepath <- "/lustre/scratch113/projects/fvg_seq/16092015/13102015_RELEASE/ALL_FVG_20151013_freq_only.frq"
	filepath <- "/lustre/scratch113/projects/fvg_seq/16092015/13102015_RELEASE/ALL_FVG_20151013_freq_only_nomono.frq"
	}

	current_pop_maf <- read.table(filepath,header=T,stringsAsFactors=F, comment.char="")
	colnames(current_pop_maf) <- c(pop)
	#create hist for counts and freqs
	current_pop_hist <- multhist(current_pop_maf,freq=FALSE,breaks=20,plot=F)
	current_pop_hist_count <- as.data.frame(cbind(current_pop_hist[[1]]$mids,current_pop_hist[[1]]$counts))
	current_pop_hist_count$pop <- pop
	colnames(current_pop_hist_count) <- c("breaks","counts","pop")
	current_pop_hist_freq <- as.data.frame(cbind(current_pop_hist[[1]]$mids,current_pop_hist[[1]]$density))
	current_pop_hist_freq$pop <- pop
	colnames(current_pop_hist_freq) <- c("breaks","freq","pop")

	all_pop_MAF_count <- rbind(all_pop_MAF_count,current_pop_hist_count)
	all_pop_MAF_freq <- rbind(all_pop_MAF_freq,current_pop_hist_freq)
	
}


save(all_pop_MAF_count,file="all_pop_MAF_count.RData")
save(all_pop_MAF_freq,file="all_pop_MAF_freq.RData")

#Upload data:
#all_maf <- read.table(filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="",colClasses=c("character","character","numeric","numeric","numeric","numeric","numeric","numeric"))
#all_maf <- read.table(filepath,header=F,sep=' ', stringsAsFactors=F, comment.char="")
# all_maf <- read.table(filepath,header=F,sep='\t', stringsAsFactors=F, comment.char="")

#FIXME: the best would be to have an headed file...
#colnames(all_maf) <- c('CHROM','POS','ID','AN','AC','DP','AF','MAF')
# colnames(all_maf) <- c('CHROM','POS','ID','AC','AN','AF','MAF')
# colnames(all_maf) <- c('CHROM','POS','REF','ALT','MAF')
#FIXME:different header for different format...maybe we could set up so sort of recognition for columns we want to use...
##################################################
#BEST FIX: provide the header as a parameter!!!###
##################################################
# colnames(all_maf) <- c('rsid','pos','allele_A','allele_B','index','average_maximum_posterior_call','info','cohort_1_AA','cohort_1_AB','all_BB','all_NULL','MAF','missing_data_proportion','cohort_1_hwe','all_A_freq','all_B_freq','minor_ALL')
#colnames(all_maf) <- c('CHROM','ID','POS','exp_freq_a1','info','certainty','type','info_type1','concord_type1','r2_type1','info_type0','concord_type0','r2_type0','MAF')
#colnames(all_maf) <- c('CHROM','ID','POS','exp_freq_a1','info','certainty','type','MAF')

#if you want to extract only SNPs after uploading the table:
all_maf_snps <- all_maf[which(nchar(all_maf$ALT) == nchar(all_maf$REF)),]

#first write a summary (for snps only)
# write.table(summary(all_maf),file="complete_summary.txt",sep="\t",col.names=T,quote=F,row.names=F)
write.table(summary(all_maf_snps),file="complete_snp_summary.txt",sep="\t",col.names=T,quote=F,row.names=F)

#now plot maf density
jpeg(paste(pop_name,"maf_density_plot.jpg",sep="_"),width=800, height=800)
par(lab=c(10,10,12),cex=2)
  plot(density(all_maf_snps$MAF), main="Maf density distribution", xlab="Maf",ylab="Maf density")
dev.off()

#count monomorphic sites
all_mono_snps <- all_maf_snps[which(all_maf_snps$MAF == 0),]
dim(all_mono_snps)
#count recursively until AC = 5
all_count_snps <- NULL

for (i in 1:50) {
  all_current_ac <- all_maf_snps[which(all_maf_snps$AC == i | all_maf_snps$AC == (all_maf_snps$AN - i)),]
  # print(paste("Current allele count:",i,sep=""))
  # print(dim(all_current_ac))
  # assign(paste("all_ac_",i,sep=""),all_current_ac)
  all_count_current <- c(i,length(all_current_ac$CHROM))
  all_count_snps <- rbind(all_count_snps,all_count_current)
}

#set the class to data.frame
all_count_snps <- as.data.frame(all_count_snps)
rownames(all_count_snps) <- NULL
colnames(all_count_snps) <- c("AC","count")

#plot those numbers
jpeg(paste(pop_name,"site_allele_count_plot.jpg",sep="_"),width=800, height=800)
  barplot(all_count_snps$count,names.arg=all_count_snps$AC,col=colors()[72])
dev.off()

#now split all maf in frequencies bins
#this is a function that splits in bins of maf and writes some summaries
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
# maf_classes <- c(0,0.005,0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)
maf_classes <- c(0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

# maf_classes <- c(0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
# maf_classes <- c(0,0.02,0.05,0.1,0.2,0.3,0.4,0.5)

all_maf_classes <- split_bins(maf_classes,all_maf_snps[which(all_maf_snps$MAF > 0),],pop_name)
# all_maf_classes <- split_bins(maf_classes,maf_table,pop_name)
gc()
#write a cute output
sink(paste(pop_name,'maf_bin_resume.txt',sep="_"))
print(all_maf_classes)
sink()

jpeg(paste(pop_name,"site_count_plot.jpg",sep="_"),width=800, height=800)
	barplot(as.matrix(all_maf_classes),names.arg=maf_classes[2:length(maf_classes)],col=colors()[72])
dev.off()

q(save='no')


##############################################################################
###13/10/2015

all_pop_MAF_freq$pop <- factor(all_pop_MAF_freq$pop, levels = pops)
colnames(all_pop_MAF_freq) <- c("breaks","value","pop")

all_pop_MAF_count$pop <- factor(all_pop_MAF_count$pop, levels = pops)
colnames(all_pop_MAF_count) <- c("breaks","value","pop")
all_cols <- col_pop(pops)
maf_sets_1 <- c("all_pop_MAF_freq")
maf_sets_2 <- c("all_pop_MAF_count")
# for(set in maf_sets_1){
for(set in maf_sets_2){
	# set <- "all_pop_MAF_freq"
  current_set <- get(set)
  pl <- ggplot(current_set)

  pl <- pl + geom_bar(stat="identity",width=0.5, position = position_dodge(width=0.8),colour="black")
  pl <- pl + aes(x = factor(breaks), y = value, fill=pop)
  pl <- pl + xlab("MAF")
  # pl <- pl + ylab("Proportion of sites (%)")
  pl <- pl + ylab("Site count")
  pl <- pl + scale_fill_manual("", values=all_cols)
  # pl <- pl + facet_wrap( ~ cat, ncol=1)
  pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
  pl <- pl + theme_bw()

  pl <- pl + theme(strip.text.x = element_text(size = 20))
  pl <- pl + theme(axis.text.x=element_text(size = rel(1.2),angle=45))
  pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
  pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
  pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
  
  # png(paste(base_folder,"/1_",set,"_20150527.png",sep=""),width=1400, height=500)
  # print(pl)
  # dev.off()
  ggsave(filename=paste(base_folder,"/1_",set,"_20150525.jpeg",sep=""),width=24, height=14,units="cm",dpi=600,plot=pl)
}

###############################PLOT OVERLAP BY MAF
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
require(plotrix)
require(ggplot2)
require(reshape2)
base_folder <- getwd()

#COMPARISON BASED ON INGI FREQS
# pops <- c("EUR")
# ipops <- c("fvg","vbi")
pops <- c("EUR","TSI","UK10K")
pops <- c("EUR","TSI")
ipops <- c("carl","fvg","vbi")
for (pop in pops){
  print(pop)
  for(ipop in ipops){
    print(ipop)
    all_current_pop <- NULL
    for (chr in seq(1,22,1)){
      print(chr)
      filepath <- paste("/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/",pop,"/UNION/",chr,"/sites_",pop,"_",ipop,"_only.txt",sep="")
      current_chr <- read.table(filepath,header=F,colClasses=c("numeric","numeric","character","character","factor","numeric","numeric"),comment.char="")
      colnames(current_chr) <- c("CHR","POS","REF","ALT","OVERLAP",pop,ipop)
      
      current_chr_ipop <- current_chr[,c(1,5,6,7)]
      all_current_pop <- rbind(all_current_pop,current_chr_ipop)
    }
    assign(paste("all_pop_",pop,"_",ipop,"_only",sep=""),all_current_pop)
    save(file=paste("all_pop_",pop,"_",ipop,"_only",sep=""),all_current_pop)
  }
 
}

#select bins of frequency and check overlapping sites
epops <- c("EUR","TSI")
epops <- c("TSI")
for (pop in epops){
  ingi_sets <- c(paste("all_pop_",pop,"_carl_only",sep=""),paste("all_pop_",pop,"_vbi_only",sep=""),paste("all_pop_",pop,"_fvg_only",sep=""))
  for (ingi in ingi_sets){
    print(ingi)
    # ingi <- ingi_sets[2]
    current_ingi <- get(ingi)
    freqs <- c(0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.30)
    current_ingi_all_bin <- NULL
    for(i in 1:length(freqs)){
      if(i==1){
        current_ingi_current_bin <- current_ingi[which(current_ingi[,4] <= freqs[i]),]
        current_out_current_bin <- current_ingi[which(current_ingi[,3] <= freqs[i]),]
      }else if(i==length(freqs)){
        current_ingi_current_bin <- current_ingi[which(current_ingi[,4] > freqs[i]),]
        current_out_current_bin <- current_ingi[which(current_ingi[,3] > freqs[i]),]
      }else{
        current_ingi_current_bin <- current_ingi[which(current_ingi[,4] > freqs[i-1] & current_ingi[,4] <= freqs[i]),]
        current_out_current_bin <- current_ingi[which(current_ingi[,3] > freqs[i-1] & current_ingi[,3] <= freqs[i]),]
      }
      current_ingi_current_bin_resume <- as.data.frame(cbind(bin=freqs[i],tot_ingi=length(current_ingi_current_bin$CHR),tot_out=length(current_out_current_bin$CHR)))

      current_ingi_all_bin <- rbind(current_ingi_all_bin,current_ingi_current_bin_resume)
    }
      current_ingi_all_bin$diff <- current_ingi_all_bin$tot_ingi - current_ingi_all_bin$tot_out
      assign(paste(ingi,"_all_bin_reverse",sep=""),current_ingi_all_bin)
  }
}



#COMPARISON BASED ON OUTBRED FREQS
for (pop in pops){
  all_chr_pop <- NULL
  # pop="UK10K"
  print(pop)
  all_pop_carl <- NULL
  all_pop_vbi <- NULL
  all_pop_fvg <- NULL
  all_pop_out <- NULL
  for (chr in seq(1,22,1)){
    print(chr)
    filepath <- paste("/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/",pop,"/UNION/",chr,"/sites_",pop,"_carl_vbi_fvg.txt",sep="")
    current_chr <- read.table(filepath,header=F,colClasses=c("numeric","numeric","character","character","factor","numeric","numeric","numeric","numeric"),comment.char="")
    colnames(current_chr) <- c("CHR","POS","REF","ALT","OVERLAP",pop,"CARL","VBI","FVG")

    current_chr_out <- current_chr[,c(1,5,6)]
    
    current_chr_pop_carl <- current_chr[which(current_chr$OVERLAP=="1111" |current_chr$OVERLAP=="1110"|current_chr$OVERLAP=="1101" |current_chr$OVERLAP=="1100" | current_chr$OVERLAP=="1100"),c(1,5,6,7)]
    current_chr_pop_vbi <- current_chr[which(current_chr$OVERLAP=="1111" |current_chr$OVERLAP=="1011" |current_chr$OVERLAP=="1010" |current_chr$OVERLAP=="1110" | current_chr$OVERLAP=="1110" ),c(1,5,6,8)]
    current_chr_pop_fvg <- current_chr[which(current_chr$OVERLAP=="1111" |current_chr$OVERLAP=="1001" |current_chr$OVERLAP=="1101" |current_chr$OVERLAP=="1011" | current_chr$OVERLAP=="1011"),c(1,5,6,9)]
    all_pop_carl <- rbind(all_pop_carl,current_chr_pop_carl)
    all_pop_vbi <- rbind(all_pop_vbi,current_chr_pop_vbi)
    all_pop_fvg <- rbind(all_pop_fvg,current_chr_pop_fvg)
    all_pop_out <- rbind(all_pop_out,current_chr_out)
  }
  assign(paste("all_pop_carl_",pop,sep=""),all_pop_carl)
  assign(paste("all_pop_vbi_",pop,sep=""),all_pop_vbi)
  assign(paste("all_pop_fvg_",pop,sep=""),all_pop_fvg)
  assign(paste("all_pop_",pop,sep=""),all_pop_out)
  
  save(file=paste("all_pop_carl_",pop,".Rdata",sep=""),all_pop_carl)
  save(file=paste("all_pop_vbi_",pop,".Rdata",sep=""),all_pop_vbi)
  save(file=paste("all_pop_fvg_",pop,".Rdata",sep=""),all_pop_fvg)
  save(file=paste("all_pop_",pop,".Rdata",sep=""),all_pop_out)

}


epops <- c("EUR","TSI")
for (pop in epops){
  ingi_sets <- c(paste("all_pop_carl_",pop,sep=""),paste("all_pop_vbi_",pop,sep=""),paste("all_pop_fvg_",pop,sep=""))
  
  for (ingi in ingi_sets){
    # ingi <- ingi_sets[1]
    current_ingi <- get(ingi)
    load(paste("all_pop_",pop,".Rdata",sep=""))
    current_out <- all_pop_out

    freqs <- c(0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.30)
    current_ingi_all_bin <- NULL
    for(i in 1:length(freqs)){
      if(i==1){
        current_ingi_current_bin <- current_ingi[which(current_ingi[,3] <= freqs[i]),]
        current_out_current_bin <- current_out[which(current_out[,3] <= freqs[i]),]
      }else if(i==length(freqs)){
        current_ingi_current_bin <- current_ingi[which(current_ingi[,3] > freqs[i]),]
        current_out_current_bin <- current_out[which(current_out[,3] > freqs[i]),]
      }else{
        current_ingi_current_bin <- current_ingi[which(current_ingi[,3] > freqs[i-1] & current_ingi[,3] <= freqs[i]),]
        current_out_current_bin <- current_out[which(current_out[,3] > freqs[i-1] & current_out[,3] <= freqs[i]),]
      }
      current_ingi_current_bin_resume <- as.data.frame(cbind(bin=freqs[i],tot_out=length(current_out_current_bin$CHR),tot_ingi=length(current_ingi_current_bin$CHR)))

      current_ingi_all_bin <- rbind(current_ingi_all_bin,current_ingi_current_bin_resume)
    }
      current_ingi_all_bin$diff <- current_ingi_all_bin$tot_out - current_ingi_all_bin$tot_ingi
      assign(paste(ingi,"_all_bin",sep=""),current_ingi_all_bin)
  }
}






ipops <- c("CARL","FVG","VBI")
for (ipop in ipops){
  ##### Plots ######
  current_ingi_all_bin <- get(paste("all_pop_",ipop,"_EUR",sep=""))
  current_ingi_all_bin <- get(paste("all_pop_",ipop,"_TSI",sep=""))

  current_ingi_all_bin <- all_pop_carl_EUR_all_bin
  current_ingi_all_bin <- all_pop_carl_TSI_all_bin

  for (i in length(epops)){
  current_ingi_all_bin <- get(paste("all_pop_",ipop,"_",pops[i],"_all_bin",sep=""))
    
  pdf(paste("/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/AF/shared_variants_EUR_",ipop,".pdf",sep=""), width=8, height=4, pointsize=12)
  # pdf(paste("/lustre/scratch113/projects/esgi-vbseq/13102015_SIGU/AF/shared_variants_EUR_","CARL",".pdf",sep=""), width=8, height=4, pointsize=12)
  m <- barplot(current_ingi_all_bin,col=c("white","red","black"),
               las=1,xlab="number of variants (x20000000)",
               xlim=c(0,sxmax+1), border=c("white","red","black"))

  }
  # uk10k <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_REL-2012-06-02_paper.txt", stringsAsFactors=F)
  # g1k <- read.delim("~/uk10k/analysis/overlaps_1000GP_uk10k_REL-2012-06-02_paper.txt", stringsAsFactors=F)

  unuk10k <- current_ingi_all_bin$diff
  # ung1k <- g1k$total - g1k$overlap
  xmax <- 30000000
  whitebar <- xmax - (current_ingi_all_bin$tot_ingi + unuk10k)
  fscale <- 20000000
  sxmax <- xmax/fscale
  # pg1k <- paste(round(100 * g1k$overlap / g1k$total,1),"%", sep="")
  puk10k <- paste(round(100 * current_ingi_all_bin$tot_ingi / current_ingi_all_bin$tot_out,1),"%", sep="")

  # par(mfrow=c(1,2), mai=c(1,1.2,0.1,0.1))
  slab <- c(0.05,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30)
  # labnames=c("singletons    ","doubletons   ","AF<=5%     ","AF<=10%    ","AF>10%     ")
  # m <- barplot(rbind(whitebar, unuk10k, current_ingi_all_bin$tot_ingi)/fscale, horiz=T, col=c("white","red","black"),
               # las=1,
  m <- barplot(rbind(whitebar, unuk10k, current_ingi_all_bin$tot_ingi)/fscale,horiz=T,col=c("white","red","black"),
               las=1,xlab="number of variants (x20000000)",
               xlim=c(0,sxmax+1), border=c("white","red","black"))
  # text(whitebar/fscale-3.5, m, puk10k)
  # axis(1, labels=slab, at=sxmax-slab)
  # barplot(rbind(g1k$overlap, ung1k)/fscale, horiz=T, names.arg=labnames, col=c("black","red"), las=1,
  #         xlab="number of variants (x1000000)             ",
  #         axes=F, xlim=c(0,sxmax+1), border=c("black","red"))
  # text((g1k$overlap+ung1k)/1000000+3.5, m, pg1k)
  # axis(1, labels=slab, at=slab)
  legend(14,5, legend=c("shared","unique"), fill=c("black","red"), cex=0.8)
  dev.off()


  
}


#need to define by hand freq bins, than select and count stuff for each freq bin for each population

  # current_chr_pop_carl_melt <- melt(current_chr_pop_carl, id='OVERLAP')

for(set in maf_sets_3){
  current_set <- get(set)
  pl <- ggplot(current_ingi_all_bin)

  pl <- pl + geom_bar(stat="identity",width=0.5, position = position_dodge(width=0.8),colour="black")
  pl <- pl + aes(x = factor(bin), y = tot_out, fill=tot_ingi)
  pl <- pl + xlab("MAF")
  pl <- pl + ylab("Proportion of sites (%)")
  # pl <- pl + ylab("Site count")
  pl <- pl + scale_fill_manual("", values=all_cols)
  # pl <- pl + facet_wrap( ~ cat, ncol=1)
  pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
  pl <- pl + theme_bw()

  pl <- pl + theme(strip.text.x = element_text(size = 20))
  pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
  pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
  pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
  pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
  
  # png(paste(base_folder,"/1_",set,"_20150912_count.png",sep=""),width=1400, height=500)
  # png(paste(base_folder,"/1_",set,"_20150914.png",sep=""),width=1400, height=500)
  # print(pl)
  # dev.off()
  ggsave(filename=paste(base_folder,"/1_",ingi,"_19102015.jpeg",sep=""),width=12, height=7,units="in",dpi=300,plot=pl)
}


}
}