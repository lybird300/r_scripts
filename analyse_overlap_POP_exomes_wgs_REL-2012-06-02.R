source("~/R_dir/R_scripts/load_nsta_functions.R")


##  a = lower bin boundary
##  b = upper bin boundary
##  maf = MAF
##  x = id for first set
##  y = id for second set
##  s = bin size for MAF quantiles
get.overlap.by.maf.bin4 <- function(a, b, maf, x, y)  { 
  inds <- which(maf>a & maf<=b)
  total <- length(inds)
  if (total!=0)  {
    num <- length(intersect(x[inds], y))
  }  else  {
    num <- 0
  }
  c(total, num)
}

#write(x[inds], "~/wgs_chr1_pos.txt")   ##  6507
#write(y, "~/exome_chr1_pos.txt")       ## ##### Exomes found in WGS by MAF #####

total <- overlap <- rep(0,14)
for (chr in 1:22)  {
  fn <- paste("~/uk10k/analysis/overlaps_uk10k_exomes_wgs_all_SNPs_REL-2012-06-02/chr", chr,".txt", sep="")
  f <- read.delim(fn)
  total <- total + f$total
  overlap <- overlap + f$overlap
}
perc <- round(100 * overlap / total,2)
dat <- data.frame(f[,1:2], total, overlap, perc)
#write.table(dat, "~/uk10k/analysis/overlaps_uk10k_exomes_wgs_all_SNPs_REL-2012-06-02_unfiltered.txt", row.names=F, quote=F, sep="\t")
#write.table(dat, "~/uk10k/analysis/overlaps_uk10k_exomes_wgs_all_SNPs_REL-2012-06-02.txt", row.names=F, quote=F, sep="\t")

qu <- round(100*c(1/7032, 2/7032, 5/7032, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, seq(0.1,0.5,0.1)),2)
#pdf("/lustre/scratch113/projects/uk10k/users/kw8/plots3/overlap_uk10k_exomes_wgs_by_maf_unfiltered.pdf", width=6.5, height=5, pointsize=12)
#pdf("/lustre/scratch113/projects/uk10k/users/kw8/plots3/overlap_uk10k_exomes_wgs_by_maf.pdf", width=6.5, height=5, pointsize=12)
plot(log10(dat$maf2)[-c(12,13)], dat$perc[-c(12,13)], col="red", pch=15, ty="b", xlab="MAF in exomes (%)", axes=F, 
     ylab="Exonic SNPs in WGS (%)", ylim=c(0,100), cex=1.3)
box(which = "plot", lty = "solid")
axis(1, labels=qu[-c(3,12,13)], at=log10(dat$maf2)[-c(3,12,13)], cex.axis=0.8)
axis(2, las=1, cex.axis=0.8)
dev.off()

## comparison unfiltered and filtered
dat1 <- read.delim("~/uk10k/analysis/overlaps_uk10k_exomes_wgs_all_SNPs_REL-2012-06-02_unfiltered.txt")
dat2 <- read.delim("~/uk10k/analysis/overlaps_uk10k_exomes_wgs_all_SNPs_REL-2012-06-02.txt")
plot(dat1$overlap, dat2$overlap, col="blue")
plot(dat1$perc, dat2$perc, col="blue")
abline(0,1)



##### Variants found in WGS but not in exomes using the baits - chr20 unfiltered #####

f <- read.delim("/lustre/scratch113/projects/uk10k/cohorts/exomes/wgs_but_not_in_exomes/chr20_unfiltered.vcf.gz", stringsAsFactors=F)   ## 11507
## remove INDELs
indrm <- grep("INDEL",f$INFO)
f <- f[-indrm,]
## get VQSLOD and HWE
vqslod <- get.vcf.info(f$INFO, "VQSLOD")
hwe <- get.vcf.info(f$INFO, "HWE")
100 * mean(vqslod < -0.6804)                         ## 8.9% 
100 * mean(hwe < 1e-6, na.rm=T)                      ## 0.3% 
100 * mean(vqslod < -0.6804 | hwe < 1e-6, na.rm=T)   ## 9.1% are filtered out with the new exclusion list 
## filtered variants
f <- f[hwe >= 1e-6 & vqslod >= -0.6804 ,]             ## 10468
ac <- get.vcf.info(f$INFO, "AC")
af <- 100 * ac/7562
summary(af)
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
##  0.01322   0.01322   0.01322   2.30300   0.02645 100.00000       130 
sum(af>0.1, na.rm=T)      ## 1258
sum(af>1, na.rm=T)        ##  896
sum(af>5, na.rm=T)        ##  673

## How many variants in WGS but not exomes
bins <- cbind(c(0,0.1,0.5,1,5), c(0.1,0.5,1,5,100))
indb <- apply(bins, 1, function(v) )


###### Extract all baits regions from WGS and calculate overlap both ways ######

## echo 'bash ~/scripts/uk10k/extract_wgs_in_exome_baits_REL-2012-06-02.sh ${LSB_JOBINDEX}' | bsub -J "chr[1-22]" -o ~/lsfout/extract_wgs_in_exome_baits_REL-2012-06-02_chr%I.out -R "select[mem>=1000] rusage[mem=1000]" -M1000

## UK10K exomes data
load("/lustre/scratch113/projects/uk10k/cohorts/exomes/exomes_final_SNVs_freqs_in_final_samples.Rdata")    ## counts
## select columns
exomes <- counts[,c(1,2,13,14)]
## remove non-polymorphic
indnp <- which(exomes$all_alt_allele_count==0)   ## 329,214 = 28%
exomes <- exomes[-indnp,]
## loop trough WGS
d <- dir("/lustre/scratch113/projects/uk10k/cohorts/exomes/wgs_in_baits", full.names=T)
qu <- cbind(c(0, 1/7032,2/7032,0.001,0.005,0.01,0.05), c(1/7032,2/7032,0.001,0.005,0.01,0.05,1))
n <- nrow(qu)
allexo <- allwgs <- matrix(0, nrow=n, ncol=2)
for (chr in 1:length(d))  {
  print(paste("chr =", chr))
  ## UK10K WGS baits
  fn <- paste("/lustre/scratch113/projects/uk10k/cohorts/exomes/wgs_in_baits/chr",chr,".vcf.gz", sep="")
  wgs <- read.delim(fn, header=F, stringsAsFactors=F)    ## 70217 for chr1
  ## get chrom for exomes
  indchr <- which(exomes$chrom==chr)
  exs <- exomes[indchr,]                                 ## 86620 for chr1
  ## combine AF for multi-allelic by replacing AC of multi-allelic in bi-allelic
  indd <- which(duplicated(exs$pos))                     ##  698 ~ 1% for chr1
  multiallelic <- exs[indd,]
  biallelic <- exs[-indd,]
  m <- merge(biallelic, multiallelic, by.x=2, by.y=2, all.x=T)    ## 85922
  m[,3] <- apply(m[,c(3,6)], 1, sum, na.rm=T)
  m <- m[,c(2,1,3,4)]
  names(m) <- names(exs)
  exs <- m                                               ## 85922
  rm(m)
  ## filter for SNPs in UK10K
  indrm <- grep("INDEL", wgs[,8])                        ##  2789
  wgs <- wgs[-indrm,]                                    ## 67428
  ## AF and MAF for exomes
  afe <- exs$all_alt_allele_count / exs$all_total_allele_count
  ## AF and MAF for WGS
  ac <- get.vcf.info(wgs[,8], "AC", num=F)
  indm <- grep(",", ac)
  multiallelic <- ac[indm]
  multiallelic <- strsplit(multiallelic, ",")
  multiallelic <- sapply(multiallelic, function(v) sum(as.numeric(v)))
  ac[indm] <- multiallelic
  ac <- as.numeric(ac)
  afw <- ac / 7562
  ## get overlap
  overexo <- overwgs <- matrix(0, nrow=n, ncol=2)
  for (i in 1:n)  {
    print(i)
    a <- qu[i,1]
    b <- qu[i,2]
    overexo[i,] <- get.overlap.by.maf.bin4(a=a, b=b, maf=afe, x=exs[,2], y=wgs[,2])
    overwgs[i,] <- get.overlap.by.maf.bin4(a=a, b=b, maf=afw, x=wgs[,2], y=exs[,2])
  }
  allexo <- allexo + overexo
  allwgs <- allwgs + overwgs
}
dat <- data.frame(100*qu, allexo, allwgs)
dat[,1] <- round(dat[,1],3)
dat[,2] <- round(dat[,2],3)
names(dat) <- c("AF1","AF2","total.ex","overlap.ex","total.wg","overlap.wg")
fn <- "~/uk10k/analysis/overlaps_uk10k_exomes_wgs_REL-2012-06-02.txt"
#write.table(dat, fn, sep="\t", row.names=F, quote=F)


##### Plots ######

dat <- read.delim("~/uk10k/analysis/overlaps_uk10k_exomes_wgs_REL-2012-06-02.txt", stringsAsFactors=F)
unwg <- dat$total.wg - dat$overlap.wg
unex <- dat$total.ex - dat$overlap.ex
xmax <- 570000
whitebar <- xmax - (dat$overlap.wg + unwg)
sxmax <- 5.7
perce <- paste(round(100 * dat$overlap.ex / dat$total.ex,1),"%", sep="")
percw <- paste(round(100 * dat$overlap.wg / dat$total.wg,1),"%", sep="")

#pdf("~/uk10k/manuscripts/main/plots/shared_variants_exomes_wgs.pdf", width=8, height=4, pointsize=12)
par(mfrow=c(1,2), mai=c(1,1.2,0.3,0.1))
labnames=c("singletons  ","doubletons ","AF<=0.1%  ","AF<=0.5%  ","AF<=1%    ","AF<=5%    ","AF>5%     ")
m <- barplot(rbind(whitebar, unwg, dat$overlap.wg)/100000, horiz=T, col=c("white","red","black"), las=1, xlab="number of variants (x100000)",
             axes=F, xlim=c(0,sxmax), border=c("white","red","black"))
text(whitebar/100000-0.6, m, percw)
axis(1, labels=0:5, at=sxmax-0:5)
barplot(rbind(dat$overlap.ex, unex)/100000, horiz=T, names.arg=labnames, col=c("black","red"), las=1, xlab="number of variants (x100000)",
        axes=F, xlim=c(0,sxmax), border=c("black","red"), legend.text=c("shared","unique"))
text((dat$overlap.ex+unex)/100000+0.6, m, perce)
axis(1)
dev.off()

##### Number of overlaps #####
sdat <- dat[,c(1:4,7,8)]

c(sum(dat$total.ex),
  sum(dat$overlap.ex),
  sum(dat$total.wg),
  sum(dat$overlap.wg))
## total  shared  total shared
## 847637 322445 662891 322445

## all shared
100*sum(dat$overlap.ex)/sum(dat$total.ex)    ## 38.0%
100*sum(dat$overlap.wg)/sum(dat$total.wg)    ## 46.6%
## shared for AF <= 0.1%
100*sum(dat$overlap.ex[1:3])/sum(dat$total.ex[1:3])    ## 180366 / 699891 = 25.8%,  699891 - 180366 = 519525
100*sum(dat$overlap.wg[1:3])/sum(dat$total.wg[1:3])    ## 179391 / 500764 = 35.8%,  500764 - 179391 = 321373
## shared for AF > 0.1%
100*sum(dat$overlap.ex[4:7])/sum(dat$total.ex[4:7])    ## 142079 / 147746 = 96.2%,  147746 - 142079 =   5667
100*sum(dat$overlap.wg[4:7])/sum(dat$total.wg[4:7])    ## 143054 / 162127 = 88.2%,  162127 - 143054 =  19073

