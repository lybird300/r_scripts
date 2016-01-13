#script to estract statisics fron vcf-stats results
rm(list=ls())

#read commandline args
args <- commandArgs(trailing=TRUE)

#Args list:
#args[[1]]:list with individuals and sex (absolute path) (i.e '/nfs/users/nfs_m/mc14/Work/SANGER/DOCS/ValBorbera/sequenced_id_conversion_sex.txt')
#args[[2]]:the path of vcf-stats files to process (absolute path)
#args[[3]]: prefix used to assign names to vcf-stats output files

#Assign args
#indiv_file <- args[[1]]
#vcfstats_filepath <- args[[2]]
#vcfstats_file_prefix <- args[[3]]
indiv_file <- '/nfs/users/nfs_m/mc14/Work/SANGER/DOCS/ValBorbera/sequenced_id_conversion_sex.txt'
vcfstats_filepath <- 'vcf_stats'
vcfstats_file_prefix <- 'generic_stats'

#upload list with sex
individuals <- read.table(indiv_file, header=F)

#set colnames based on col number
if ( dim(individuals)[2] == 2 ){
colnames(individuals) <- c('SEQ_ID','sex')
}

if ( dim(individuals)[2] == 3 ){
colnames(individuals) <- c('SEQ_ID','PRIVATE_ID','sex')
}
#We assume by default that the last column is the sex column,
# and it is coded as follow:
#MALE:1
#FEMALE:2

#separate male/female
male_inds <- individuals[which(individuals$sex == 1),]
female_inds <- individuals[which(individuals$sex == 2),]

#Open a file stream to write some outputs
#sink('stats_resumes.log')

#INDELS
indels <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'indels',sep="."),sep="/"), header=F, skip=2,col.names=c('count','sample'))
indels$sample <- gsub("samples/","",indels$sample)

indels_summary <- summary(indels)

#cat('INDELS\n')
#cat(indels_summary,labels=names(indels_summary))
#cat('\n\n')

write.table(indels_summary, file='indels_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("indels_density.jpg",width=1000, height=1000)
plot(density(indels$count),main ="Indels density distribution in VBSEQ", xlab="Indels count")
dev.off()

jpeg("indels_boxplot.jpg", width=1000, height=1000)
indels_boxplot <- boxplot(indels$count,main ="Indels in VBSEQ", xlab="All individuals", ylab="Indels count")
dev.off()

#cat('SNPs Outliers\n')
#cat(indels[which(indels$count %in%  indels_boxplot$out),]$sample)
#cat('\n\n')

#to retrieve the outliers from boxplot
indels_outliers <- indels[which(indels$count %in% indels_boxplot$out),]$sample

#male/female boxplot
jpeg("indels_boxplot_male_female.jpg", width=1000, height=1000)
indels_boxplot_mf <- boxplot(indels[which(indels$sample %in% male_inds$SEQ_ID),]$count,indels[which(indels$sample %in% female_inds$SEQ_ID),]$count,main ="Indels in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Indels count")
dev.off()

#to retrieve the outliers from boxplot
#male_indels_outliers <- indels[which(indels$sample %in% male_inds$SEQ_ID && indels$count %in% indels_boxplot$out),]$sample
#female_indels_outliers <- indels[which(indels$count %in% indels_boxplot$out),]$sample

#######PRIVATE SNPS
private_s <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'private',sep="."),sep="/"), header=F, skip=1,col.names=c('count','sample'))

private_s_summary <- summary(private_s)

#cat('PRIVATE SNPS\n')
#cat(private_s_summary)
#cat('\n\n')

write.table(private_s_summary, file='private_snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("private_snps_density.jpg", width=1000, height=1000)
plot(density(private_s$count),main ="Private SNPs density distribution in VBSEQ", xlab="Private SNPs count")
dev.off()

jpeg("private_snps_boxplot.jpg", width=1000, height=1000)
private_s_boxplot <- boxplot(private_s$count,main ="Private SNPs in VBSEQ", xlab="All individuals", ylab="Private SNPs count")
dev.off()

#to retrieve the outliers from boxplot
private_outliers <- private_s[which(private_s$count %in%  private_s_boxplot$out),]$sample

#cat('PRIVATE SNPs Outliers\n')
#cat(private_s[which(private_s$count %in%  private_s_boxplot$out),]$sample)
#cat('\n\n')
#sink()

#male/female boxplot
jpeg("private_snps_boxplot_male_female.jpg", width=1000, height=1000)
private_s_boxplot_mf <- boxplot(private_s[which(private_s$sample %in% male_inds$SEQ_ID),]$count,private_s[which(private_s$sample %in% female_inds$SEQ_ID),]$count,main ="Private SNPs in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Private SNPs count")
dev.off()

################SHARED SNPS
#Number of sites having a non-reference allele in 0,1,2,etc samples

shared_s <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'shared',sep="."),sep="/"), header=F, skip=1,col.names=c('ind_count','freq'))

shared_s_summary <- summary(shared_s$freq)

#write.table(shared_s_summary, file='shared_snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("shared_snps_density.jpg", width=1000, height=1000)
plot(density(shared_s$freq),main ="Density distribution of the number of sites having a non-reference allele in VBSEQ",xlab="SNPs count")
dev.off()

jpeg("shared_snps_hist.jpg", width=1000, height=1000)
barplot(shared_s$freq,names.arg=shared_s$ind_count,main ="Number of sites having a non-reference allele in VBSEQ", xlab="# Individuals", ylab="Number of sites having a non-reference allel (count)")
dev.off()

#SNPS
#Number of positions with SNPs

snps <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'snps',sep="."),sep="/"), header=F, skip=2,col.names=c('count','sample'))
snps$sample <- gsub("samples/","",snps$sample)

snps_summary <- summary(snps$count)

#write.table(snps_summary, file='snps_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("snps_density.jpg", width=1000, height=1000)
plot(density(snps$count),main ="SNPs density distribution in VBSEQ",xlab="SNPs count")
dev.off()

jpeg("snps_boxplot.jpg", width=1000, height=1000)
snps_boxplot <- boxplot(snps$count,main ="SNPs in VBSEQ", xlab="All individuals", ylab="SNPs count")
dev.off()

#to retrieve the outliers from boxplot
snps_outliers <- snps[which(snps$count %in%  snps_boxplot$out),]$sample

#male/female boxplot
jpeg("snps_boxplot_male_female.jpg", width=1000, height=1000)
snps_boxplot_mf <- boxplot(snps[which(snps$sample %in% male_inds$SEQ_ID),]$count,snps[which(snps$sample %in% female_inds$SEQ_ID),]$count,main ="SNPs in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="SNPs count")
dev.off()


###############SAMPLE TS/TV
#samples-tstv

s_tstv <- read.table(paste(vcfstats_filepath,paste(vcfstats_file_prefix,'samples-tstv',sep="."),sep="/"), header=F,skip=1,col.names=c('Transitions','Transversions','ts_tv_ratio','sample'))

s_tstv_summary <- summary(s_tstv)

#write.table(s_tstv_summary, file='s_tstv_summary.txt',append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NaN", dec =".", row.names = FALSE, col.names = TRUE)

jpeg("s_ts_density.jpg", width=1000, height=1000)
plot(density(s_tstv$Transitions),col='red',main ="Transitions density distribution in VBSEQ",xlab="Transition count")
dev.off()

jpeg("s_tv_density.jpg", width=1000, height=1000)
plot(density(s_tstv$Transversions),col='green',main ="Transversions density distribution in VBSEQ",xlab="Transversions count")
dev.off()

jpeg("s_tstv_density.jpg", width=1000, height=1000)
plot(density(s_tstv$ts_tv_ratio),col='blue',main ="Ts/Tv ratio density distribution in VBSEQ",xlab="Ts/Tv ratio")
dev.off()

#boxplots
jpeg("s_ts_boxplot.jpg", width=1000, height=1000)
s_ts_boxplot <- boxplot(s_tstv$Transitions,main ="Transition in VBSEQ",ylab="Transition count", xlab="All individuals", names=c("Transitions"))
dev.off()

jpeg("s_tv_boxplot.jpg", width=1000, height=1000)
s_tv_boxplot <- boxplot(s_tstv$Transversions,main ="Transversion in VBSEQ",ylab="Transversion count", xlab="All individuals", names=c("Transversions"))
dev.off()

jpeg("s_tstv_boxplot.jpg", width=1000, height=1000)
s_tstv_boxplot <- boxplot(s_tstv$ts_tv_ratio,main ="Ts/Tv ratio in VBSEQ",ylab="Ts/Tv ratio", xlab="All individuals", names=c("Ts/Tv"))
dev.off()

#male/female boxplot
jpeg("s_ts_boxplot_male_female.jpg", width=1000, height=1000)
s_ts_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$Transitions,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$Transitions,main ="Transitions in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Transitions count")
dev.off()

jpeg("s_tv_boxplot_male_female.jpg", width=1000, height=1000)
s_tv_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$Transversions,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$Transversions,main ="Transversions in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Transversions count")
dev.off()

jpeg("s_tstv_boxplot_male_female.jpg", width=1000, height=1000)
s_tstv_boxplot_mf <- boxplot(s_tstv[which(s_tstv$sample %in% male_inds$SEQ_ID),]$ts_tv_ratio,s_tstv[which(s_tstv$sample %in% female_inds$SEQ_ID),]$ts_tv_ratio,main ="Ts/Tv ratio in VBSEQ", names=c("MALE individuals","FEMALE individuals"), ylab="Ts/Tv ratio")
dev.off()

#to retrieve the outliers from boxplot
s_ts_outliers <- s_tstv[which(s_tstv$Transitions %in%  s_ts_boxplot$out),]$sample
s_tv_outliers <- s_tstv[which(s_tstv$Transversions %in%  s_tv_boxplot$out),]$sample
s_tstv_outliers <- s_tstv[which(s_tstv$ts_tv_ratio %in%  s_tstv_boxplot$out),]$sample


##############################################################################################
# plots stats for new seq and heterozygosys
cohort <- "CARL"
cohort <- "FVG"
rm(list=ls())
cohort <- "VBI"

het_all <- NULL
het_all_merged <- NULL
for (chr in 1:22){
print(chr)
current_chr <- paste(chr,".vcf.gz.snps.stats.psc.tab",sep="")
het_current <- read.table(current_chr,header=T)
if (chr==1){
het_all_merged <- het_current
}else{
het_all_merged <- merge(het_all_merged,het_current, by="sample")
}
}

het_all <- as.data.frame(het_all_merged$sample)
colnames(het_all) <- c("sample")
het_all$nRefHom <- apply(het_all_merged[,grep("nRefHom",colnames(het_all_merged))],1,sum) 
het_all$nNonRefHom <- apply(het_all_merged[,grep("nNonRefHom",colnames(het_all_merged))],1,sum) 
het_all$nHets <- apply(het_all_merged[,grep("nHets",colnames(het_all_merged))],1,sum)
het_all$nSingletons <- apply(het_all_merged[,grep("nSingletons",colnames(het_all_merged))],1,sum)
het_all$totVar <- apply(het_all[,c("nRefHom","nNonRefHom","nHets")],1,sum)
het_all$het_rate <- (het_all$nHets / (het_all$totVar))*100
het_all$singleton_rate <- (het_all$nSingletons / (het_all$totVar))
het_all$prop_alt_het<- (het_all$nHets/(het_all$nNonRefHom))
het_all$het_hom <- het_all$nHets/(het_all$nRefHom+het_all$nNonRefHom)

###############################PLOT ALTERNATIVE HET proportion

#calculate stdev to select outliers
sd_het <- sd(het_all$prop_alt_het)
mean_het <- mean(het_all$prop_alt_het)
sd3u_het <- mean_het + 3*sd_het
sd3l_het <- mean_het - 3*sd_het
het_all$sample <- as.character(het_all$sample)


summaries <- as.data.frame(cbind(val=c(mean_het,sd3u_het,sd3l_het),names=c("mean","+3sd","-3sd")))
summaries$val <- as.numeric(as.character(summaries$val))
# outliers <- (het_all[which(het_all$het_rate > (mean_het + 3*sd_het) | het_all$het_rate < (mean_het - 3*sd_het)),1])

# for (i in 1:length(het_all$outliers)){
# if (!het_all[i,]$outliers %in% outliers){
# het_all[i,]$outliers <- NA
#   }
# }

require(ggplot2)
pl <- ggplot(het_all)

pl <- pl + geom_point(size=3,colour="red")
pl <- pl + aes(x = factor(sample), y = prop_alt_het)
pl <- pl + geom_hline(data=summaries,aes(yintercept=val),linetype=8,colour=names)
pl <- pl + ylim(1,4.5)

pl <- pl + xlab("Samples")
pl <- pl + ylab("Het/nonRefHom sites")
pl <- pl + theme_bw()
pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
pl <- pl + geom_text(aes(label=ifelse(prop_alt_het > sd3u_het |prop_alt_het < sd3l_het,sample,'')),hjust=0.5,vjust=1.5)

base_folder <- getwd()
png(paste(base_folder,"/",cohort,"_het_nonREFH.png",sep=""),width=1400, height=500)
print(pl)
dev.off()

###############################PLOT STANDARD HET proportion
s_sd_het <- sd(het_all$het_rate)
s_mean_het <- mean(het_all$het_rate)
s_sd3u_het <- s_mean_het + 3*s_sd_het
s_sd3l_het <- s_mean_het - 3*s_sd_het

s_summaries <- as.data.frame(cbind(val=c(s_mean_het,s_sd3u_het,s_sd3l_het),names=c("mean","+3sd","-3sd")))
s_summaries$val <- as.numeric(as.character(s_summaries$val))

pl <- ggplot(het_all)

pl <- pl + geom_point(size=3,colour="red")
pl <- pl + aes(x = factor(sample), y = het_rate)
pl <- pl + geom_hline(data=s_summaries,aes(yintercept=val),linetype=8,colour=names)
pl <- pl + ylim(0,21)
pl <- pl + xlab("Samples")
pl <- pl + ylab("Proportion of heterozygous sites (%)")
pl <- pl + theme_bw()
pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
pl <- pl + geom_text(aes(label=ifelse(het_rate > s_sd3u_het |het_rate < s_sd3l_het,sample,'')),hjust=0.5,vjust=1.5)

base_folder <- getwd()
png(paste(base_folder,"/",cohort,"_het.png",sep=""),width=1400, height=500)
print(pl)
dev.off()


##############################################################################################
# Manhattan plot of something by chromosome: hard code the name of the column you want to plot against the position


# cohorts <- c("VBI","FVG","CARL")
# plot_columns <- c("lp_HWE","QUAL","VQSLOD")
cohort <- "VBI"
cohort <- "FVG"

require(ggplot2)
cohort <- "CARL"
# plot_columns <- c("O_HET")
plot_columns <- c("lp_HWE")

#we want to plot all by chromosome, but also to have some stats genome wide
#so we need to merge the data together, calculate the desired stats and then plot by chr
data_all <- NULL
for (chr in 1:22){
        # chr=22
        print(chr)
        current_chr <- paste(chr,".vcf.gz_qvh.tab",sep="")
        # current_chr <- paste(chr,".vcf.gz.snps.hardy.hwe.idxed",sep="")
        data_current <- read.table(current_chr,header=T) 
        data_current$HWE <- as.numeric(as.character(data_current$HWE))
        data_current$lp_HWE <- -log10(data_current$HWE)
        data_all <- rbind(data_all,data_current)

    if(chr%%6 == 0 | chr == 22){
        print("PLOT!!")
        for (plot_column in plot_columns){
            # pl <- ggplot(data_current)
            pl <- ggplot(data_all)

            pl <- pl + geom_point(size=2,aes(colour=factor(CHROM)))
            pl <- pl + aes_string(x = 'POS', y = plot_column)
            # pl <- pl + aes_string(x = 'IDX', y = plot_column)

            pl <- pl + xlab("Position")
            pl <- pl + ylab(plot_column)
            pl <- pl + theme_bw()
            pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
            # pl <- pl + scale_y_continuous()
            # pl <- pl + geom_text(aes(label=ifelse(het_rate > (mean_het + 3*sd_het) |het_rate < (mean_het - 3*sd_het),sample,'')),hjust=0.5,vjust=-0.5)
            pl <- pl + facet_wrap(~ CHROM,nrow=3,scales = "free_x")
            base_folder <- getwd()
            png(paste(base_folder,"/",cohort,"_manhattan_",plot_column,"_",chr,".png",sep=""),width=1400, height=500)
            print(pl)
            dev.off()
        }
      data_all <- NULL  
    }
}

################################################################################
## Plot 
# 1)MAF wgs vs MAF gwas 
# 2)MAF wgs vs NRD by site
# 3)MAF gwas vs NRD by site

# 1)MAF wgs vs MAF gwas 
# cohorts <- c("FVG","CARL")
rm(list=ls())
require(ggplot2)
cohorts <- c("VBI","FVG","CARL")
#cohorts <- c("FVG","CARL")
# cohort <- "VBI"
# cohort <- "FVG"
# cohort <- "CARL"
for(cohort in cohorts){
    print(cohort)
    if(cohort == "VBI"){
        gwas_set <- read.table("/lustre/scratch113/projects/esgi-vbseq/08092015/13102015_RELEASE/REL_CHECK/MAF/ALL_VB_20151013_gwas_def.s1.freq.frq", header=T)
        wgs_set <- read.table("/lustre/scratch113/projects/esgi-vbseq/08092015/13102015_RELEASE/REL_CHECK/MAF/ALL_VB_20151013_seq.s2.freq.frq", header=T)
        nrd_pop <- "/lustre/scratch113/projects/esgi-vbseq/08092015/13102015_RELEASE/REL_CHECK/MAF/site_concordance_discordance_table_chr"
    
    }else if(cohort == "FVG"){
        gwas_set <- read.table("/lustre/scratch113/projects/fvg_seq/16092015/13102015_RELEASE/REL_CHECK/MAF/ALL_FVG_20151013_gwas_last.s1.freq.frq", header=T)
        wgs_set <- read.table("/lustre/scratch113/projects/fvg_seq/16092015/13102015_RELEASE/REL_CHECK/MAF/ALL_FVG_20151013_seq_last_definitivo.s2.freq.frq", header=T)
        nrd_pop <- "/lustre/scratch113/projects/fvg_seq/16092015/13102015_RELEASE/REL_CHECK/MAF/site_concordance_discordance_table_chr"

    }else if(cohort == "CARL"){
        gwas_set <- read.table("/lustre/scratch113/projects/carl_seq/variant_refinement/13102015_RELEASE/REL_CHECK/MAF/ALL_CARL_20151013_gwas_last.s1.freq.frq",header=T)
        wgs_set <- read.table("/lustre/scratch113/projects/carl_seq/variant_refinement/13102015_RELEASE/REL_CHECK/MAF/ALL_CARL_20151013_seq_last_definitivo.s2.freq.frq",header=T)
        nrd_pop <- "/lustre/scratch113/projects/carl_seq/variant_refinement/13102015_RELEASE/REL_CHECK/MAF/site_concordance_discordance_table_chr"
    }

    gwas_wgs <- merge(gwas_set,wgs_set, by="SNP")
    
    #MAF.x > maf gwas, MAF.y > maf wgs
    gwas_wgs$CHR.y <- NULL
    gwas_wgs$A1.x <- as.character(gwas_wgs$A1.x)
    gwas_wgs$A1.y <- as.character(gwas_wgs$A1.y)
    gwas_wgs$A2.x <- as.character(gwas_wgs$A2.x)
    gwas_wgs$A2.y <- as.character(gwas_wgs$A2.y)
    gwas_wgs$NCHROBS.y <- NULL


    pl <- ggplot(gwas_wgs)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'MAF.x', y = 'MAF.y')
    pl <- pl + xlab("MAF gwas")
    pl <- pl + ylab("MAF wgs")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_mafs.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()

    # 2)MAF wgs vs NRD by site
    # get all NRD by site
    NRD_all <- NULL
    for(chr in 1:22){
        print(chr)
        NRD_current <- read.table(paste(nrd_pop,chr,".txt",sep=""),header=T)
        NRD_all <- rbind(NRD_all,NRD_current)
    }

    #output some stats
    #1) how many sites above the average NRD
    #1) how many sites above the median NRD
    # add lines for those on the plots
    gwas_wgs_nrd <- merge(gwas_wgs,NRD_all, by="SNP")
    gwas_wgs_nrd$maf_diff <- gwas_wgs_nrd$MAF.x - gwas_wgs_nrd$MAF.y

    gwas_wgs_nrd$CHR <- NULL
    
    #gwas_wgs_nrd_all_match <- gwas_wgs_nrd[which(gwas_wgs_nrd$A1.x == gwas_wgs_nrd$A1.y & gwas_wgs_nrd$A2.x == gwas_wgs_nrd$A2.y ),]
    #gwas_wgs_nrd_all_match[which(gwas_wgs_nrd_all_match$NRD > 0.5),]
    
    #mismatch <- gwas_wgs_nrd[!(gwas_wgs_nrd$SNP %in% gwas_wgs_nrd_all_match$SNP),]
    #mismatch[which(mismatch$NRD > 0.1),]

    NRD_median <- median(gwas_wgs_nrd$NRD,na.rm=T)
    NRD_mean <- mean(gwas_wgs_nrd$NRD,na.rm=T)
    NRD_sd <- sd(gwas_wgs_nrd$NRD,na.rm=T)
    NRD_3sdu <- NRD_mean + 3*NRD_sd
    NRD_3sdl <- NRD_mean - 3*NRD_sd

    NRD_summary <- summary(gwas_wgs_nrd$NRD)
    #Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    #0.000000 0.007634 0.012990 0.018620 0.020410 1.000000

    out_median <- dim(gwas_wgs_nrd[which(gwas_wgs_nrd$NRD > NRD_median),])[1]
    out_mean <- dim(gwas_wgs_nrd[which(gwas_wgs_nrd$NRD > NRD_mean),])[1]
    out_3sd <- dim(gwas_wgs_nrd[which(gwas_wgs_nrd$NRD > (NRD_mean + 3*NRD_sd) | gwas_wgs_nrd$NRD < (NRD_mean - 3*NRD_sd)),])[1]
    
    outliers_3sd <- (gwas_wgs_nrd[which(gwas_wgs_nrd$NRD > (NRD_mean + 3*NRD_sd) | gwas_wgs_nrd$NRD < (NRD_mean - 3*NRD_sd)),])

    ###write a file for outliers
    sink(paste(base_folder,"/",cohort,"_outliers_NRD.txt",sep=""))
    cat(paste('NRD mean:',NRD_mean,'\n',sep=''))
    cat(paste('NRD median:',NRD_median,'\n',sep=''))
    cat(paste('NRD sd:',NRD_sd,'\n',sep=''))
    cat('\nOutliers NRD > mean\n')
    cat(out_median)
    cat('\nOutliers NRD > median\n')
    cat(out_mean)
    cat('\nOutliers NRD > (<) 3 sd\n')
    cat(out_3sd)
    cat('\n')
    sink()

    summaries <- as.data.frame(cbind(val=c(NRD_mean,NRD_3sdu,NRD_3sdl),names=c("mean","+3sd","-3sd")))
    summaries$val <- as.numeric(as.character(summaries$val))

    #2)MAF gwas vs NRD by site

    pl <- ggplot(gwas_wgs_nrd)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'MAF.x', y = 'NRD')
    pl <- pl + geom_hline(data=summaries, aes(yintercept=val),linetype=8,colour=names)
    pl <- pl + xlab("MAF gwas")
    pl <- pl + ylab("NRD")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_gwas_maf_nrd.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()
    
    # 3)MAF wgs vs NRD by site
    # get all NRD by site

    pl <- ggplot(gwas_wgs_nrd)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'MAF.y', y = 'NRD')
    pl <- pl + geom_hline(data=summaries, aes(yintercept=val),linetype=8,colour=names)
    pl <- pl + xlab("MAF gwas")
    pl <- pl + ylab("NRD")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_seq_maf_nrd.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()

    # 4)diff MAF gwas/seq vs NRD by site
    # get all NRD by site

    pl <- ggplot(gwas_wgs_nrd)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'maf_diff', y = 'NRD')
    pl <- pl + geom_hline(data=summaries, aes(yintercept=val),linetype=8,colour=names)
    pl <- pl + xlab("Maff difference (Gwas - Seq")
    pl <- pl + ylab("NRD")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_diff_maf_nrd.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()

    #################### OUTLIERS PLOTS ###############
    #
    #1a)MAF gwas vs MAF Wgs by site OUTLIERS
    pl <- ggplot(outliers_3sd)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'MAF.x', y = 'MAF.y')
    pl <- pl + xlab("MAF gwas")
    pl <- pl + ylab("MAF wgs")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_mafs_out.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()

    #2a)MAF gwas vs NRD by site OUTLIERS
    pl <- ggplot(outliers_3sd)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'MAF.x', y = 'NRD')
    pl <- pl + xlab("MAF gwas")
    pl <- pl + ylab("NRD")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_gwas_maf_nrd_out.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()

    # 3a)MAF wgs vs NRD by site OUTLIERS
    # get all NRD by site
    pl <- ggplot(outliers_3sd)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'MAF.y', y = 'NRD')
    pl <- pl + xlab("MAF gwas")
    pl <- pl + ylab("NRD")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_seq_maf_nrd_out.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()

    # 4a)diff MAF gwas/seq vs NRD by site OUTLIERS
    # get all NRD by site
    pl <- ggplot(outliers_3sd)
    pl <- pl + geom_point(size=2,aes(colour=factor(CHR.x)))
    pl <- pl + aes_string(x = 'maf_diff', y = 'NRD')
    pl <- pl + xlab("Maff difference (Gwas - Seq")
    pl <- pl + ylab("NRD")
    pl <- pl + theme_bw()
    pl <- pl + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
    pl <- pl + theme(strip.text.x = element_text(size = 20))
    pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
    pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
    pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_blank())

    base_folder <- getwd()
    png(paste(base_folder,"/",cohort,"_diff_maf_nrd_out.png",sep=""),width=1400, height=500)
    print(pl)
    dev.off()

    ###write the table for OUTLIERS 
    outliers_3sd_table <- (within(outliers_3sd,POS <- as.data.frame(do.call('rbind',strsplit(as.character(POS),':',fixed=T)))$V2))
    colnames(outliers_3sd_table) <- c("SNP","CHR","A1.g","A2.g","MAF.g","NCHROBS.g","A1.w","A2.w","MAF.w","POS","OGC","NRD","maf_diff")

    write.table(outliers_3sd_table,file=paste(base_folder,"/",cohort,"_out_table.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)


}


###write a file for outliers
sink('Outliers.txt')
cat('Outliers for INDELS\n')
cat(indels_outliers)
cat('\n\nOutliers for PRIVATE SNPs\n')
cat(as.character(private_outliers))
cat('\n\nOutliers for SNPs\n')
cat(snps_outliers)
cat('\n\nOutliers for Transitions\n')
cat(as.character(s_ts_outliers))
cat('\n\nOutliers for Transversions\n')
cat(as.character(s_tv_outliers))
cat('\n\nOutliers for Ts/Tv ratio\n')
cat(as.character(s_tstv_outliers))
cat('\n')
sink()


