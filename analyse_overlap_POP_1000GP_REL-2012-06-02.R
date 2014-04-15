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


##### Overlap with 1000GP by given MAF 0.1%, 0.2%, 0.5%, 1%, 2%, 5%, 10%, each with a very narrow bin around it of + or - 10% #####
### Combine all chroms just for SNPs at all sites
total <- over <- p <- NULL
for (chr in 1:22)  {
  fn <- paste("~/uk10k/analysis/overlaps_uk10k_1000GP_chroms_all_SNPs_REL-2012-06-02/chr", chr, ".txt", sep="")
  f <- read.delim(fn)
  total <- cbind(total, f[,3])
  over <- cbind(over, f[,4])
  p <- cbind(p, f[,5])
}
mtotal <- apply(total, 1, sum)
mover <- apply(over, 1, sum)
mp <- round(100 * mover / mtotal,2)
dat <- data.frame(c(0.1,0.2,0.5,1,2,5,10), mp, mover, mtotal)
names(dat) <- c("MAF","overlap(%)", "overlap","total")
#write.table(dat, "~/uk10k/analysis/overlaps_uk10k_1000GP_all_sites_REL-2012-06-02.txt", row.names=F, quote=F, sep="\t")

##    MAF overlap(%)   ## REL-2011-12-01
## 1  0.1      46.00
## 2  0.2      62.00
## 3  0.5      87.61
## 4  1.0      95.82
## 5  2.0      97.39
## 6  5.0      97.70
## 7 10.0      97.95

##    MAF overlap(%) overlap  total    ## REL-2012-06-02
## 1  0.1      42.94  347881 810139
## 2  0.2      55.24  247765 448505
## 3  0.5      79.08  265464 335702
## 4  1.0      89.79  289880 322845
## 5  2.0      92.74  289916 312624
## 6  5.0      94.86  307883 324567
## 7 10.0      96.57  365244 378209


### Combine all chroms just for SNPs at P-sites
total <- over <- p <- NULL
for (chr in 1:22)  {
  fn <- paste("~/uk10k/analysis/overlaps_uk10k_1000GP_chroms_P-SNPs_REL-2012-06-02/chr", chr, ".txt", sep="")
  f <- read.delim(fn)
  total <- cbind(total, f[,3])
  over <- cbind(over, f[,4])
  p <- cbind(p, f[,5])
}
mtotal <- apply(total, 1, sum)
mover <- apply(over, 1, sum)
mp <- round(100 * mover / mtotal,2)
dat <- data.frame(c(0.1,0.2,0.5,1,2,5,10), mp, mover, mtotal)
names(dat) <- c("MAF","overlap(%)", "overlap","total")
#write.table(dat, "~/uk10k/analysis/overlaps_uk10k_1000GP_P-sites_REL-2012-06-02.txt", row.names=F, quote=F, sep="\t")


### Plots for Nicole ###
f1 <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_all_sites.txt", stringsAsFactors=F)
f2 <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_P-sites.txt", stringsAsFactors=F)
f3 <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_all_sites_REL-2012-06-02.txt", stringsAsFactors=F)
f4 <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_P-sites_REL-2012-06-02.txt", stringsAsFactors=F)
#pdf("/lustre/scratch113/projects/uk10k/users/kw8/plots/overlap_uk10k_g1k_by_maf_REL-2012-06-02.pdf", width=5, height=5, pointsize=13)
plot(log10(f1$MAF), f1[,2], col="blue", ty="n", pch=15, axes=F, xlab="MAF in UK10K", ylab="Percentage found in 1000GP", ylim=c(0,105), cex=1.3)
points(log10(f1$MAF), f1[,2], col="darkblue", ty="b", pch=15, lty=2, cex=1.3)
points(log10(f2$MAF), f2[,2], col="darkred", ty="b", pch=16, lty=3, cex=1.3)
points(log10(f3$MAF), f3[,2], col="blue", ty="b", pch=15, lty=4, cex=1.3)
points(log10(f4$MAF), f4[,2], col="red", ty="b", pch=16, lty=5, cex=1.3)
legend(log10(0.16), 34, legend=c("all sites (REL-2011-12-01)", "easy access sites (REL-2011-12-01)", "all sites (REL-2012-06-02)", "easy access sites (REL-2012-06-02)"),
       col=c("darkblue","darkred","blue","red"), pch=c(15,16), lty=c(2:5), cex=0.85)
box(which = "plot", lty = "solid")
axis(1, labels=f1$MAF, at=log10(f1$MAF), cex.axis=1.1)
axis(2, las=1, cex.axis=1.1)
dev.off()

cbind(f1[,1:2], f3[,2])



##### Compare to Jie's numbers for chr20 #####

f <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_chroms_all_SNPs_bins_nums/chr20.txt")
lab <- c("singletons","(singletons,0.1%]","(0.1%,1%]","(1%,5%]","(5%,50%]")
#pdf("~/uk10k/plots/overlap_uk10k_g1k_nums.pdf", width=6.5, height=5, pointsize=10)
barplot(as.matrix(t(f[,3:4])), beside=T, col=c("darkblue","lightblue"), names.arg=lab)
legend(10,250000, legend=c("UK10K total calls","shared with 1000GP"), fill=c("darkblue","lightblue"))
dev.off()



##### Overlap 1000GP with UK10K #####

total <- over <- p <- NULL
for (chr in 1:22)  {
  fn <- paste("~/uk10k/analysis/overlaps_1000GP_uk10k_chroms_all_SNPs_REL-2012-06-02/chr", chr, ".txt", sep="")
  f <- read.delim(fn)
  total <- cbind(total, f[,3])
  over <- cbind(over, f[,4])
  p <- cbind(p, f[,5])
}
mtotal <- apply(total, 1, sum)
mover <- apply(over, 1, sum)
mp <- round(100 * mover / mtotal,2)
dat <- data.frame(c(0.1,0.2,0.5,1,2,5,10,20,30,40,50), mp, mover, mtotal)
names(dat) <- c("MAF","overlap(%)", "overlap","total")
#write.table(dat[-c(1:3),], "~/uk10k/analysis/overlaps_1000GP_uk10k_all_sites_REL-2012-06-02.txt", row.names=F, quote=F, sep="\t")


##### Plots #####

over1 <- read.delim("~/uk10k/analysis/overlaps_1000GP_uk10k_all_sites_REL-2012-06-02.txt")
over2 <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_all_sites_REL-2012-06-02.txt")
over1 <- over1[1:4,]
over2 <- over2[4:7,]
## question: use MAF from UK10K or MAF from 1000GP as reference?



##### For UK10K manuscript #####

### UK10K in 1000GP 
total <- over <- matrix(0, nrow=5, ncol=22)
for (chr in c(1:22))  {
  print(chr)
  fn <- paste("~/uk10k/analysis/overlaps_uk10k_1000GP_REL-2012-06-02_paper/chr", chr, ".txt", sep="")
  f <- read.delim(fn)
  print(f)
  total[,chr] <- f[,3]
  over[,chr] <- f[,4]
}
mtotal <- apply(total, 1, sum)
mover <- apply(over, 1, sum)
mp <- round(100 * mover / mtotal,2)
dat <- data.frame(f[,2], mp, mover, mtotal)
names(dat) <- c("AF","overlap(%)", "overlap","total")
#write.table(dat, "~/uk10k/analysis/overlaps_uk10k_1000GP_REL-2012-06-02_paper.txt", row.names=F, quote=F, sep="\t")

### 1000GP in UK10K
total <- over <- NULL
for (chr in 1:22)  {
  fn <- paste("~/uk10k/analysis/overlaps_1000GP_uk10k_REL-2012-06-02_paper/chr", chr, ".txt", sep="")
  f <- read.delim(fn)
  total <- cbind(total, f[,3])
  over <- cbind(over, f[,4])
}
mtotal <- apply(total, 1, sum)
mover <- apply(over, 1, sum)
mp <- round(100 * mover / mtotal,2)
dat <- data.frame(f[,2], mp, mover, mtotal)
names(dat) <- c("AF","overlap(%)","overlap","total")
#write.table(dat, "~/uk10k/analysis/overlaps_1000GP_uk10k_REL-2012-06-02_paper.txt", row.names=F, quote=F, sep="\t")

##### Table #####

uk10k <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_REL-2012-06-02_paper.txt", stringsAsFactors=F)
g1k <- read.delim("~/uk10k/analysis/overlaps_1000GP_uk10k_REL-2012-06-02_paper.txt", stringsAsFactors=F)
an <- 7562
uac <- round(c(1,2,0.05*an,0.1*an,an))
an <- 2184
gac <- round(c(1,2,0.05*an,0.1*an,an))
dat <- data.frame(uac,uk10k[c(1,4,3,2)], gac,g1k[,c(1,4,3,2)])
names(dat) <- c("AC","AF","total","shared with 1000GP","shared (%)","AC","AF","total","shared with UK10K","shared (%)")
fdat <- apply(dat, 2, sum)
fdat <- rbind(dat, fdat)
#write.table(fdat, "~/uk10k/manuscripts/main/tables/Table_Comparison_with_1000GP.txt", quote=F, row.names=F, sep="\t")



##### Plots ######

uk10k <- read.delim("~/uk10k/analysis/overlaps_uk10k_1000GP_REL-2012-06-02_paper.txt", stringsAsFactors=F)
g1k <- read.delim("~/uk10k/analysis/overlaps_1000GP_uk10k_REL-2012-06-02_paper.txt", stringsAsFactors=F)
unuk10k <- uk10k$total - uk10k$overlap
ung1k <- g1k$total - g1k$overlap
xmax <- 30000000
whitebar <- xmax - (uk10k$overlap + unuk10k)
fscale <- 1000000
sxmax <- xmax/fscale
pg1k <- paste(round(100 * g1k$overlap / g1k$total,1),"%", sep="")
puk10k <- paste(round(100 * uk10k$overlap / uk10k$total,1),"%", sep="")

#pdf("~/uk10k/manuscripts/main/plots/shared_variants_uk10k_1000GP_paper.pdf", width=8, height=4, pointsize=12)
par(mfrow=c(1,2), mai=c(1,1.2,0.1,0.1))
slab <- c(0,5,10,15,20)
labnames=c("singletons    ","doubletons   ","AF<=5%     ","AF<=10%    ","AF>10%     ")
m <- barplot(rbind(whitebar, unuk10k, uk10k$overlap)/fscale, horiz=T, col=c("white","red","black"), las=1, xlab="number of variants (x1000000)",
             axes=F, xlim=c(0,sxmax+1), border=c("white","red","black"))
text(whitebar/fscale-3.5, m, puk10k)
axis(1, labels=slab, at=sxmax-slab)
barplot(rbind(g1k$overlap, ung1k)/fscale, horiz=T, names.arg=labnames, col=c("black","red"), las=1,
        xlab="number of variants (x1000000)             ",
        axes=F, xlim=c(0,sxmax+1), border=c("black","red"))
text((g1k$overlap+ung1k)/1000000+3.5, m, pg1k)
axis(1, labels=slab, at=slab)
legend(14,5, legend=c("shared","unique"), fill=c("black","red"), cex=0.8)
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


