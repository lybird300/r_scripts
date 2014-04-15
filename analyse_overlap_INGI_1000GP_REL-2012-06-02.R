source("~/R_dir/R_scripts/load_nsta_functions.R")


##### Prepare data #####

## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/
## ftp-trace.wip.ncbi.nlm.nih.gov/1000genomes/ftp/technical/working/20120524_phase1_data/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz
## ftp-trace.wip.ncbi.nlm.nih.gov/1000genomes/ftp/technical/working/20120524_phase1_data/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz.tbi
## Same data on:
## http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/
  
## Split 1000GP Phase I sites into chroms
## echo '~/scripts/uk10k/split_phaseI_1000GP_into_chroms.sh ${LSB_JOBINDEX}' | bsub -J "split_chr[1-22]" -o split_chr%I.out

## Annotate UK10K sites with 1000GP genome mask
## echo '~/scripts/uk10k/annotate_vcf_strict_mask.sh ${LSB_JOBINDEX}' | bsub -J "chr[1-22]" -o chr%I.out

##### Strict genome mask from 1000GP #####
## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20120417_phase1_masks/20120515_phase1_masks.README
## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20120417_phase1_masks/StrictMask/
## N - the base is an N in the reference genome GRCh37
## L - depth of coverage is much lower than average
## H - depth of coverage is much higher than average
## Z - too many reads with zero mapping quality overlap this position
## Q - the average mapping quality at the position is too low
## P - the base passed all filters
## 0 - an overlapping base was never observed in aligned reads

## Compute overlap
## ~/R_dir/R_scripts_uk10k_2013/compute_overlap_uk10k_1000GP_REL-2012-06-02.R.R
## ~/R_dir/R_scripts_uk10k_2013/compute_overlap_uk10k_1000GP_P-sites_REL-2012-06-02.R.R


##### Description for 1000GP #####

## We computed the overlap of UK10K SNP sites of the release
## REL-2011-12-01 comprising 2432 samples with the Phase I release of the
## 1000 Genomes Project for different minor allele frequencies, 0.1%,
## 0.2%, 0.5%, 1%, 2%, 5%, 10%, each with a very narrow bin around it of
## + or - 10% (Table X, Fig X). We also repeated the same analysis (Table
## Y, Fig Y) but using only the sites where the base passed all filters
## (P sites).


#Analysis (plot mainly) of overlapping sites between INGI pop and TGP and UK10K
##### Overlap with 1000GP by given MAF 0.1%, 0.2%, 0.5%, 1%, 2%, 5%, 10%, each with a very narrow bin around it of + or - 10% #####
pop_table <- read.table("fvg_all_no_vqslod.txt",sep="\t",header=T)

#remove monomorphic sites
pop_table_no_mono <- pop_table[-which(pop_table$ALT_AF ==1 | pop_table$ALT_AF == 0),]
pop_table_mono <- pop_table[which(pop_table$ALT_AF ==1 | pop_table$ALT_AF == 0),]

source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
# maf_classes <- c(0,0.001,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
maf_classes <- c(0,0.005,0.01,0.02,0.05,0.10,0.20,0.30)
# maf_classes <- c(0,0.02,0.05,0.1,0.2,0.3,0.4,0.5)
population <- "FVG"
all_maf_classes <- split_bins(maf_classes,pop_table_no_mono,population)

maf_shared <- NULL
for(i in 2:(length(maf_classes) + 1)){
    class_maf_name <- paste(population,'_maf_lte_',maf_classes[i],sep='')
    filename <- paste(class_maf_name,'table.txt',sep='_')
    
    if (maf_classes[i - 1] == maf_classes[length(maf_classes)]){
      if (maf_classes[length(maf_classes)] != 0.5) {        
        class_maf_name <- paste(population,'_maf_gte_',maf_classes[i-1],sep='')
        filename <- paste(class_maf_name,'table.txt',sep='_')
        maf_table <- read.table(filename,sep="\t",header=T)
        i=i-1
      }else{ break}
    }else{
      maf_table <- read.table(filename,sep="\t",header=T)
    }
  #extract overlapping sites
  tgp_maf <- maf_table[which(maf_table$TGP==1),]
  uk10k_maf <- maf_table[which(maf_table$UK10K==1),]

  total <- length(maf_table$POS)
  over_tgp <- length(tgp_maf$POS)
  over_uk10k <- length(uk10k_maf$POS)

  p_tgp <- round(100*over_tgp/total,2)
  p_uk10k <- round(100*over_uk10k/total,2)

  maf_shared <- rbind(maf_shared,c(maf_classes[i],total,over_tgp,over_uk10k,p_tgp,p_uk10k))
  }
maf_shared <- as.data.frame(maf_shared)
colnames(maf_shared) <- c("MAF","total","overlap_tgp","overlap_uk10k","overlap_tgp(%)","overlap_uk10k(%)")

##### Plots ######

unuk10k <- maf_shared$total - maf_shared$overlap_uk10k
ung1k <- maf_shared$total - maf_shared$overlap_tgp
xmax <- 3000000
whitebar <- xmax - (maf_shared$overlap_uk10k + unuk10k)
fscale <- 100000
sxmax <- xmax/fscale

pg1k <- paste(round(100 * maf_shared$overlap_tgp / maf_shared$total,1),"%", sep="")
puk10k <- paste(round(100 * maf_shared$overlap_uk10k / maf_shared$total,1),"%", sep="")

#pdf("~/uk10k/manuscripts/main/plots/shared_variants_uk10k_1000GP_paper.pdf", width=8, height=4, pointsize=12)
par(mfrow=c(1,2), mai=c(1,1.2,0.1,0.1))
slab <- c(0,5,10,15,20)
# labnames=c("singletons    ","doubletons   ","AF<=5%     ","AF<=10%    ","AF>10%     ")
# maf_classes <- c(0,0.005,0.01,0.02,0.05,0.10,0.20,0.30)
labnames=c("AF<=0.005","1%  ","2% ","5% ","10%  ","20%  ","30%  ","AF>30% ")
m <- barplot(rbind(whitebar, unuk10k, maf_shared$overlap_uk10k)/fscale, horiz=T, col=c("white","red","black"), las=1, xlab="number of variants (x1000)",
             axes=F, xlim=c(0,sxmax+1), border=c("white","red","black"))
axis(1, labels=slab, at=sxmax-slab)
text(whitebar/200000+2, m, puk10k)
barplot(rbind(maf_shared$overlap_tgp, ung1k)/fscale, horiz=T, names.arg=labnames, col=c("black","red"), las=1,
        xlab="number of variants (x1000)             ",
        axes=F, xlim=c(0,sxmax+1), border=c("black","red"))
axis(1, labels=slab, at=slab)
text((maf_shared$overlap_tgp+ung1k)/80000+2, m, pg1k)
legend(14,2, legend=c("shared","unique"), fill=c("black","red"), cex=0.8)

dev.off()

##### Number of overlaps #####
c(sum(g1k$total),
  sum(g1k$overlap),
  sum(uk10k$total),
  sum(uk10k$overlap))
## exomes         wgs
##    total   shared    total   shared
## 36666891 15155075 40057860 15182676

## all shared                               
100*sum(g1k$overlap)/sum(g1k$total)        ## 41.3%
100*sum(uk10k$overlap)/sum(uk10k$total)    ## 37.9%


