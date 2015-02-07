source("~/R_dir/R_scripts/load_nsta_functions.R")


### Recode genotypes
## v = gt vector
recode.gwas.gts <- function(ref, alt, v)  {
  homRef <- paste(ref, ref, sep="")
  het <- paste(ref, alt, sep="")
  homAlt <- paste(alt, alt, sep="")
  w <- NULL
  w[v==homRef] <- 0
  w[v==het] <- 1
  w[v==homAlt] <- 2
  w[v=="00"] <- NA
  as.numeric(w)
}



##### Compare low-coverage genotypes with genotypes from exome sequencing #####
##Extract only GENOTYPES from VCF files..
##
## cp /nfs/team151/ym3/exomes/stat/pfizer-pain/PfizerExomes_Discovery_UK10K_overlap.recode.vcf /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome
## cp /nfs/team151/ym3/exomes/stat/pfizer-pain/PfizerExomes_Replication_UK10K_overlap.recode.vcf /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome
## grep #CHROM /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap.recode.vcf | cut -f10- > ~/uk10k/data/pfizer_exomes_discovery_samples.txt
## grep #CHROM /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap.recode.vcf  | cut -f10- > ~/uk10k/data/pfizer_exomes_replication_samples.txt
## zcat /lustre/scratch106/projects/uk10k/RELEASE/UK10K_COHORT/REL-2011-12-01/v3/20.beagle.impute2.anno.csq.vcf.gz | head -n100 | grep #CHROM | cut -f10- > ~/uk10k/data/samples_REL-2011-12-01_final.txt

## In total 190 samples overlap
sdat <- read.delim("~/uk10k/data/uk10k_pfizer_overlaps.txt", stringsAsFactors=F)
#write.table(sdat[,2], "~/uk10k/data/uk10k_pfizer_overlaps_uk10k_id.txt", row.names=F, col.names=F, sep="\t", quote=F)

######  USELESS for VBI   ######################################################################################################################################
##### Problems with BGI exome samples #####
## 32355, 10142, 32355, 51656, 28652, 35438, 35623
s <- c("32355","10142","32355","51656","28652","35438","35623")
s1 <- scan("~/uk10k/data/pfizer_exomes_discovery_samples.txt", what="char")
s2 <- scan("~/uk10k/data/pfizer_exomes_replication_samples.txt", what="char")

#intersect() like %in%
intersect(s, s1)    ## "32355" "10142"
intersect(s, s2)    ## "28652"
samples <- read.delim("~/uk10k/data/uk10k_pfizer_overlaps.txt")       ## 190
ms <- merge(s, samples, by.x=1, by.y=3)
## 32355     59022 QTL218676              ## to be removed, was already identified as a problem
indrm <- which(samples[,2]=="QTL218676")
#write.table(samples[-indrm,], "~/uk10k/data/uk10k_pfizer_overlaps_filtered.txt", row.names=F, sep="\t", quote=F)        ## 189
#write.table(samples[-indrm,2], "~/uk10k/data/uk10k_pfizer_overlaps_id_filtered.txt", row.names=F, col.names=F sep="\t", quote=F)    ## 189


## Check sample overlap 
sdat <- read.delim("~/uk10k/data/uk10k_pfizer_overlaps_filtered.txt", stringsAsFactors=F)   ##  189
s1 <- scan("~/uk10k/data/pfizer_exomes_discovery_samples.txt", what="char")        ##   81
s2 <- scan("~/uk10k/data/pfizer_exomes_replication_samples.txt", what="char")      ##  108
sall <- scan("~/uk10k/data/samples_REL-2011-12-01_final.txt", what="char")         ## 2432
length(intersect(sdat[,1], s1))   ##   0
length(intersect(sdat[,1], s2))   ##   0
length(intersect(sdat[,3], s1))   ##   80
length(intersect(sdat[,3], s2))   ##  108
setdiff(sdat[,3], union(s1, s2))
setdiff(union(s1, s2), sdat[,3])  ## "32355"
length(unique(sdat[,3]))          ##  188
sdat[duplicated(sdat[,3]),]
sdat[which(sdat[,1]==17012),]
## 17012   QTL212453  30419
## 17012 UK10K_30419  30419
length(intersect(sdat[,2], sall))           ## 142
any(sall=="QTL212453")      ## FALSE
any(sall=="UK10K_30419")    ## TRUE
## samples common in REL-2012-12-01 and Exome
inds <- sapply(intersect(sdat[,2], sall), function(v) which(sdat[,2]==v))
fdat <- sdat[inds,]
#write.table(fdat, "~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01.txt", row.names=F, quote=F, sep="\t")                                                    ## 142
#write.table(intersect(fdat[,3],s1), "~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01_bgi_id_discovery.txt", row.names=F, col.names=F, quote=F, sep="\t")    ##  61
#write.table(intersect(fdat[,3],s2), "~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01_bgi_id_replication.txt", row.names=F, col.names=F, quote=F, sep="\t")  ##  81
#write.table(fdat[,2], "~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01_uk10k_id.txt", row.names=F, col.names=F, quote=F, sep="\t")                          ## 142
#################################################################################################################################################################


##### Select samples which are in common #####

## bsub -o subset_discovery.out "vcf-subset -c ~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01_bgi_id_discovery.txt /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap.recode.vcf.gz | bgzip -c > /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples.recode.vcf.gz"
## bsub -o subset_replication.out "vcf-subset -c ~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01_bgi_id_replication.txt /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap.recode.vcf.gz | bgzip -c > /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples.recode.vcf.gz"
## bsub -o subset_uk10k.out "vcf-subset -c ~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01_uk10k_id.txt /lustre/scratch106/projects/uk10k/RELEASE/UK10K_COHORT/REL-2011-12-01/v3/20.beagle.impute2.anno.csq.vcf.gz | bgzip -c > /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples.vcf.gz"


##### Select chr20 from Exomes #####

## tabix -p vcf /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples.recode.vcf.gz
## tabix -p vcf /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples.recode.vcf.gz
## tabix /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples.recode.vcf.gz 20 > temp_disc_genotypes.vcf
## tabix /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples.recode.vcf.gz 20 > temp_repl_genotypes.vcf
## zcat /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples.recode.vcf.gz | head -n100 | grep "#CHROM" > temp_disc_head.txt
## zcat /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples.recode.vcf.gz | head -n100 | grep "#CHROM" > temp_repl_head.txt
## cat temp_disc_head.txt temp_disc_genotypes.vcf | bgzip -c > /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples_chr20.recode.vcf.gz
## cat temp_repl_head.txt temp_repl_genotypes.vcf | bgzip -c > /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples_chr20.recode.vcf.gz


##### Select sites which are in common #####

f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples_chr20.recode.vcf.gz", stringsAsFactor=F)   ## 4920 x 70
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples_chr20.recode.vcf.gz", stringsAsFactor=F) ## 6141 x 90
sites <- union(f1[,2], f2[,2])   ## 8436
dat <- data.frame(rep(20, length(sites)), sites)
names(dat) <- c("chr","pos")
#write.table(dat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr_sites_union_discovery_replication.txt", row.names=F, quote=F, sep="\t")     ## 8436
#write.table(dat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr_sites_union_discovery_replication_no-header.txt", row.names=F, col.names=F, quote=F, sep="\t")     ## 8436

## extract Pfizer exome sites from UK10K low-coverage using vcftools -  did not work for keeping INFO
## bsub -o subset_sites.out "vcftools --gzvcf /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples.vcf.gz --positions /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr_sites_union_discovery_replication.txt --recode-INFO-all --out /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_test"

## test file for vcftools
# head /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr_sites_union_discovery_replication.txt > test_positions_chr20.txt
# zcat /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples.vcf.gz | head -n100000 | awk '$2<=138165' | bgzip -c > test.vcf.gz
# vcftools --gzvcf test.vcf.gz --positions test_positions_chr20.txt --recode-INFO-all --recode --out test.out


### Remove INDELs
f <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples.vcf.gz", skip=55)   ## 802480, unique 
ind <- which(duplicated(f[,2]))
indd <- unlist(sapply(ind, function(v) which(f[v,2]==f[,2])))
indrm <- grep("INDEL", f[ind,][indd,8])
fdat <- rbind(f[-ind,], f[ind,][-indrm,])
ord <- order(fdat[,1], fdat[,2])
fdat <- fdat[ord,]

## this one was used
## bsub -o select_common_samples_uk10k.out "zcat /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples.vcf.gz | perl scripts/uk10k/select_common_sites_uk10k.pl - /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr_sites_union_discovery_replication_no-header.txt > /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples_sites_test.vcf"
## bgzip /lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples_sites.vcf



##### Compare genotypes in discovery file #####

samples <- read.delim("~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01.txt")      ## 143 => 142
#lc <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples_sites.vcf.recode.vcf.gz", skip=56, stringsAsFactors=F)            ## 5106 x 152
lc <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples_sites.vcf.gz", skip=56, stringsAsFactors=F)                        ## 5106 x 152
ex <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples_chr20.recode.vcf.gz", stringsAsFactor=F)    ## 4920 x  71
## sanity checks
mdat <- merge(lc[,1:8], ex[,1:8], by.x=2, by.y=2)    ## 3445
all(mdat$ID.x==mdat$ID.y)     ## FALSE
all(mdat$REF.x==mdat$REF.y)   ## FALSE    4 not equal, INDELs vs SNPs
all(mdat$ALT.x==mdat$ALT.y)   ## FALSE  462 not equal
#mdat[mdat$REF.x!=mdat$REF.y,]
#mdat[mdat$ALT.x!=mdat$ALT.y,][1:10,]
sum((mdat$ALT.x=="." | mdat$ALT.y==".") & mdat$ALT.x!=mdat$ALT.y)                                                                  ## 452 sites where one of ALT == .
#mdat[setdiff(which(mdat$ALT.x!=mdat$ALT.y), which((mdat$ALT.x=="." | mdat$ALT.y==".") & mdat$ALT.x!=mdat$ALT.y)),c(3:5,10:12)]     ## 7 INDELs, 2 multi-allelic, 1 where one LC-ALT=T and EX-ALT=A
## mark multi-allelic sites
mlcref <- grep(",", lc$REF)    ## 0
mlcalt <- grep(",", lc$ALT)    ## 0
mexref <- grep(",", ex$REF)    ## 0
mexalt <- grep(",", ex$ALT)    ## 2 
## exclude INDELs
ilcref <- which(nchar(lc$REF)>1)    ## 11
ilcalt <- which(nchar(lc$ALT)>1)    ## 18   => one in common, so 28 INDELs
iexref <- which(nchar(ex$REF)>1)    ##  0
iexalt <- which(nchar(ex$ALT)>1)    ##  2  (= multi-allelic)
lc <- lc[nchar(lc$REF)==1 & nchar(lc$ALT)==1,]
ex <- ex[nchar(ex$REF)==1 & nchar(ex$ALT)==1,]
## match sites for genotypes
mgeno <- merge(lc[,-c(1,3:9)], ex[,-c(1,3:9)], by.x=1, by.y=1)    ## 3433 x 205, check which ones belong to LC and which ones to EXOME
lcgeno <- mgeno[,c(1:144)]
exgeno <- mgeno[,c(1,145:205)]
## match samples and remove QTL218676 (= X32355)
exgeno <- exgeno[,!names(exgeno)=="X32355"]
tlcgeno <- data.frame(names(lcgeno), t(lcgeno))
texgeno <- data.frame(names(exgeno), t(exgeno))
texgeno[,1] <- gsub("X", "", texgeno[,1])
texgeno[-1,1] <- as.character(samples[unlist(sapply(as.character(texgeno[,1]), function(v) which(v==as.character(samples[,3])))),2])
tmgeno <- merge(tlcgeno, texgeno, by.x=1, by.y=1)
fgeno <- t(tmgeno)
k <- (((nrow(fgeno)-1)/2)+1)
all(fgeno[2:k,1]==fgeno[(k+1):nrow(fgeno),1])     ## TRUE
flcgeno <- fgeno[2:k,-1]                          ## 3433 x 61
fexgeno <- fgeno[(k+1):nrow(fgeno),-1]            ## 3433 x 61
## extract genotype quality and depth
exdp <- fexgeno
exgq <- fexgeno
for (i in 1:ncol(fexgeno))  {
  v <- vecstr2liststr(fexgeno[,i], ":")
  exgq[,i] <- sapply(v, function(w) w[2])
  exgq[exgq[,i]==".",i] <- NA
  exdp[,i] <- sapply(v, function(w) w[3])
  exdp[exdp[,i]==".",i] <- NA
}
fexgq <- data.frame(lcgeno[,1], exgq)
fexdp <- data.frame(lcgeno[,1], exdp)
names(fexgq) <- names(fexdp) <- fgeno[1,]
#write.table(fexgq, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_gq_discovery.txt", quote=F, row.names=F, sep="\t")        ## 3433 sites x 61 samples
#write.table(fexdp, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_dp_discovery.txt", quote=F, row.names=F, sep="\t")        ## 3433 sites x 61 samples

## extract genotype and edit to 00, 01, 11
for (i in 1:ncol(flcgeno))  {
  v <- vecstr2liststr(flcgeno[,i], ":")
  flcgeno[,i] <- sapply(v, function(w) w[1])
}
for (i in 1:ncol(fexgeno))  {
  v <- vecstr2liststr(fexgeno[,i], ":")
  fexgeno[,i] <- sapply(v, function(w) w[1])
}
flcgeno <- apply(flcgeno, 2, function(w) gsub("\\|","", w))
flcgeno <- apply(flcgeno, 2, function(w) gsub("/","", w))
flcgeno <- apply(flcgeno, 2, function(w) gsub("10","01", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("\\|","", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("/","", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("10","01", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("\\.\\.",NA, w))    
table(flcgeno)
##     00     01     11 
## 174359  24521  13966
## 171537  24138  13738    ## after removing sample QTL218676 (= X32355)
table(fexgeno)
##     ..     00     01     11 
##   5507 169820  24056  13463    ## change .. to NA
##        167087  23638  13259    ## after removing sample QTL218676 (= X32355)
## type of discordance
gmat <- matrix(NA, nrow=nrow(flcgeno), ncol=ncol(flcgeno))
for (i in 1:ncol(flcgeno))  {
  gmat[,i] <- paste(flcgeno[,i], fexgeno[,i], sep="_")
}
sort(table(gmat))
##  00_11  11_01  01_00  01_11  11_NA  00_01  01_NA  00_NA  11_11  01_01  00_00 
##     26     68    151    196    633    660    881   3915  13037  22910 166936       ##    24 + 196 + 660 = 880   missed in LC,     68 + 151 = 219   overcall in LC
sort(round(100 * table(gmat)/length(gmat),2)) 
##  11_01  00_11  01_00  01_11  00_01  11_NA  01_NA  11_11  00_NA  01_01  00_00 
##   0.01   0.03   0.07   0.09   0.30   0.32   0.42   1.87   6.23  10.94  79.72       ##                     0.4% missed in LC,                 0.01% overcall in LC

## save genotypes with position
flcgeno <- apply(flcgeno, 2, function(v) gsub("00", 0, v))
flcgeno <- apply(flcgeno, 2, function(v) gsub("01", 1, v))
flcgeno <- apply(flcgeno, 2, function(v) gsub("11", 2, v))
fexgeno <- apply(fexgeno, 2, function(v) gsub("00", 0, v))
fexgeno <- apply(fexgeno, 2, function(v) gsub("01", 1, v))
fexgeno <- apply(fexgeno, 2, function(v) gsub("11", 2, v))
flcdat <- data.frame(lcgeno[,1], flcgeno)
fexdat <- data.frame(lcgeno[,1], fexgeno)
fdat <- data.frame(lcgeno[,1], gmat)
names(fdat) <- names(flcdat) <- names(fexdat) <- fgeno[1,]
#write.table(flcdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_gts_discovery.txt", row.names=F, quote=F, sep="\t")        ## 3433 sites x 62 samples
#write.table(fexdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_exome_gts_discovery.txt", row.names=F, quote=F, sep="\t")     ## 3433 sites x 62 samples
#write.table(fdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_gts_discovery.txt", row.names=F, quote=F, sep="\t")    ## 3433 sites x 62 samples

## compare genotypes
gtc <- matrix(NA, nrow=nrow(flcgeno), ncol=ncol(flcgeno))
for (i in 1:ncol(flcgeno))  {
  gtc[,i] <- flcgeno[,i]==fexgeno[,i]
}
## save genotypes comparison with position
fdat <- data.frame(lcgeno[,1], gtc)
names(fdat) <- fgeno[1,]
#write.table(fdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", row.names=F, quote=F, sep="\t")    ## 3433 sites x 62 samples

table(gtc)
##  FALSE   TRUE 
##   1715 205624     ## 99.2% concordance

## number and percentage of concordant sites and samples
csites <- apply(gtc, 1, sum, na.rm=T)
csamples <- apply(gtc, 2, sum, na.rm=T)
msites <- 100*apply(gtc, 1, mean, na.rm=T)
msamples <- 100*apply(gtc, 2, mean, na.rm=T)
round(summary(csites))
summary(csamples)
summary(msites)
summary(msamples)
round(summary(100-msamples),1)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       1      61      62      60      62      62
##    2741    3281    3342    3317    3367    3397  
##       1.60   98.40  100.00   96.61  100.00  100.00     ## sites
##      79.80   95.58   97.30   96.60   98.10   99.00     ## samples
## after removing NAs
##       1      61      62      60      62      62 
##    2741    3281    3342    3317    3367    3397 
##      17.70   98.40  100.00   99.13  100.00  100.00     ## sites
##      81.70   99.40   99.50   99.18   99.60   99.80     ## samples
##       0.2     0.4     0.5     0.8     0.6    18.3      ## discordance by sample
sum(msites<90)      ## 232 (= 6.8%)     ##  => 27 (=0.8%) after removing NAs
sum(msamples<90)    ##   1 (= 1.6%)
round(summary(100-msamples),2)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    0.24    0.42    0.51    0.54    0.63    0.94 

## number and percentage of N/A genotypes
nacsites <- apply(gtc, 1, function(v) sum(is.na(v)))
nacsamples <- apply(gtc, 2, function(v) sum(is.na(v)))
namsites <- round(100*apply(gtc, 1, function(v) mean(is.na(v))),1)
namsamples <- round(100*apply(gtc, 2, function(v) mean(is.na(v))),1)
round(summary(nacsites))
round(summary(nacsamples))
round(summary(namsites),1)
round(summary(namsamples),1)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0       0       0       2       0      61 
##      15      52      73      89     125     216 
##     0.0     0.0     0.0     2.6     0.0    98.4     ## sites
##     0.4     1.5     2.1     2.6     3.6     6.3     ## samples

#pdf("~/uk10k/plots/discordance_lc_exome_by_site_chr20_discovery.pdf", width=4, height=4, pointsize=10)
plot(sort(100-msites), col="red", pch=15, ylim=c(0,100), main="Discordance by Site", xlab="sites", ylab="discordance (%)", cex.lab=1.3, cex.main=1.3)
dev.off()
#pdf("~/uk10k/plots/discordance_lc_exome_by_sample_chr20_discovery.pdf", width=4, height=4, pointsize=10)
plot(sort(100-msamples), col="red", pch=15, main="Discordance by Sample", xlab="samples", ylab="discordance (%)", cex.lab=1.3, cex.main=1.3)
abline(h=0.9, lty=2, col="darkgrey")
text(15, 2, "discordance = 0.9%", col="darkgrey")
text(45, 16.5, "discordance = 18.3%", col="darkgrey")
text(44, 18, "sample QTL218676", col="darkgrey")
dev.off()
#pdf("~/uk10k/plots/discordance_lc_exome_by_sample_chr20_discovery_v2.pdf", width=4, height=4, pointsize=10)
plot(sort(100-msamples), col="red", pch=15, main="Discordance by Sample", xlab="samples", ylab="discordance (%)", cex.lab=1.3, cex.main=1.3)
dev.off()
#pdf("~/uk10k/plots/na_lc_exome_by_site_chr20_discovery.pdf", width=4, height=4, pointsize=10)
plot(sort(namsites), col="blue", pch=15, main="N/A by Site", xlab="sites", ylab="N/A (%)", cex.lab=1.3, cex.main=1.3)
dev.off()
#pdf("~/uk10k/plots/na_lc_exome_by_sample_chr20_discovery.pdf", width=4, height=4, pointsize=10)
plot(sort(namsamples), col="blue", pch=15, main="N/A by Sample", xlab="samples", ylab="N/A (%)", ylim=c(0,16), cex.lab=1.3, cex.main=1.3)
dev.off()

## check high discordance sample
fgeno[1,-1][which(msamples<90)]    ## QTL218676
samples[samples[,2]=="QTL218676",]
##  X.twin_id    uk10k_id   bgi_id
##      59022   QTL218676    32355

## check samples with high NAs
fgeno[1,-1][which(namsamples>3)]
##  [1] "QTL_RP333422" "QTL191237"    "QTL191786"    "QTL192008"    "QTL210326"   
##  [6] "QTL210500"    "QTL210566"    "QTL210576"    "QTL210586"    "QTL210779"   
## [11] "QTL210781"    "QTL211629"    "QTL211981"    "QTL212089"    "QTL212364"   
## [16] "QTL212935"    "QTL218573"    "QTL218997"    "UK10K_119845"


##### Compare genotypes in replication file ######

samples <- read.delim("~/uk10k/data/uk10k_pfizer_overlaps_REL-2012-12-01.txt")       ## 143
lc <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples_sites.vcf.recode.vcf.gz", skip=56, stringsAsFactors=F)                      ## 5106 x 152
ex <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples_chr20.recode.vcf.gz", stringsAsFactor=F)           ## 6141 x 90
## sanity checks
mdat <- merge(lc[,1:8], ex[,1:8], by.x=2, by.y=2)    ## 3926
all(mdat$ID.x==mdat$ID.y)     ## FALSE
all(mdat$REF.x==mdat$REF.y)   ## FALSE    9 not equal, i.e. INDEL vs SNP
all(mdat$ALT.x==mdat$ALT.y)   ## FALSE  382 not equal, either ALT is unknown or INDEL vs SNP
#mdat[mdat$REF.x!=mdat$REF.y,]
#mdat[mdat$ALT.x!=mdat$ALT.y,][1:10,]
sum((mdat$ALT.x=="." | mdat$ALT.y==".") & mdat$ALT.x!=mdat$ALT.y)                                                                  ## 368 sites where one of ALT == .
#mdat[setdiff(which(mdat$ALT.x!=mdat$ALT.y), which((mdat$ALT.x=="." | mdat$ALT.y==".") & mdat$ALT.x!=mdat$ALT.y)),c(3:5,10:12)]     ## 11 INDELs, 2 multi-allelic, 1 where one LC-ALT=A and EX-ALT=C,
## mark multi-allelic sites                                                                                                        ## but is A in dbSNP !!!
mlcref <- grep(",", lc$REF)    ## 0
mlcalt <- grep(",", lc$ALT)    ## 0
mexref <- grep(",", ex$REF)    ## 0
mexalt <- grep(",", ex$ALT)    ## 4 
## exclude INDELs
ilcref <- which(nchar(lc$REF)>1)    ## 11
ilcalt <- which(nchar(lc$ALT)>1)    ## 18
iexref <- which(nchar(ex$REF)>1)    ##  0
iexalt <- which(nchar(ex$ALT)>1)    ##  4  (= multi-allelic)
lc <- lc[nchar(lc$REF)==1 & nchar(lc$ALT)==1,]
ex <- ex[nchar(ex$REF)==1 & nchar(ex$ALT)==1,]
## match sites for genotypes
mgeno <- merge(lc[,-c(1,3:9)], ex[,-c(1,3:9)], by.x=1, by.y=1)    ## 3903 x 225, check which ones belong to LC and which ones to EXOME
lcgeno <- mgeno[,c(1:144)]
exgeno <- mgeno[,c(1,145:225)]
## match samples
tlcgeno <- data.frame(names(lcgeno), t(lcgeno))
texgeno <- data.frame(names(exgeno), t(exgeno))
texgeno[,1] <- gsub("X", "", texgeno[,1])
texgeno[-1,1] <- as.character(samples[unlist(sapply(as.character(texgeno[,1]), function(v) which(v==as.character(samples[,3])))),2])
tmgeno <- merge(tlcgeno, texgeno, by.x=1, by.y=1)
fgeno <- t(tmgeno)
k <- (((nrow(fgeno)-1)/2)+1)
all(fgeno[2:k,1]==fgeno[(k+1):nrow(fgeno),1])     ## TRUE
flcgeno <- fgeno[2:k,-1]                          ## 3903 x 81
fexgeno <- fgeno[(k+1):nrow(fgeno),-1]            ## 3903 x 81
## extract genotype quality and depth
exdp <- fexgeno
exgq <- fexgeno
for (i in 1:ncol(fexgeno))  {
  v <- vecstr2liststr(fexgeno[,i], ":")
  exgq[,i] <- sapply(v, function(w) w[2])
  exgq[exgq[,i]==".",i] <- NA
  exdp[,i] <- sapply(v, function(w) w[3])
  exdp[exdp[,i]==".",i] <- NA
}
fexgq <- data.frame(lcgeno[,1], exgq)
fexdp <- data.frame(lcgeno[,1], exdp)
names(fexgq) <- names(fexdp) <- fgeno[1,]
#write.table(fexgq, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_gq_replication.txt", quote=F, row.names=F, sep="\t")        ## 3903 sites x 82 samples
#write.table(fexdp, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_dp_replication.txt", quote=F, row.names=F, sep="\t")        ## 3903 sites x 82 samples

## extract genotype and edit to 00, 01, 11
for (i in 1:ncol(flcgeno))  {
  v <- vecstr2liststr(flcgeno[,i], ":")
  flcgeno[,i] <- sapply(v, function(w) w[1])
}
for (i in 1:ncol(fexgeno))  {
  v <- vecstr2liststr(fexgeno[,i], ":")
  fexgeno[,i] <- sapply(v, function(w) w[1])
}
flcgeno <- apply(flcgeno, 2, function(w) gsub("\\|","", w))
flcgeno <- apply(flcgeno, 2, function(w) gsub("/","", w))
flcgeno <- apply(flcgeno, 2, function(w) gsub("10","01", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("\\|","", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("/","", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("10","01", w))
fexgeno <- apply(fexgeno, 2, function(w) gsub("\\.\\.",NA, w))    
table(flcgeno)
##     00     01     11
## 174359  24521  13966           ## LC het / homAlt = 1.756 (discovery)  and 1.720 (replication)
## 256543  37686  21914           ## change from discovery to replication set by factor 1.471, 1.537, 1.569
table(fexgeno)
##     ..     00     01     11 
##   5507 169820  24056  13463    ## change .. to NA,  EX het / homAlt = 1.787 (discovery)  and 1.780 (replication)
##        236403  33812  18993    ## change from discovery to replication set by factor 1.392, 1.406, 1.411

## type of discordance
gmat <- matrix(NA, nrow=nrow(flcgeno), ncol=ncol(flcgeno))
for (i in 1:ncol(flcgeno))  {
  gmat[,i] <- paste(flcgeno[,i], fexgeno[,i], sep="_")
}
sort(table(gmat))
##  11_01  00_11  01_00  01_11  00_01  11_NA  01_NA  11_11  00_NA  01_01  00_00 
##     60     91    200    213   1147   3165   4668  18689  19102  32605 236203         =>    91 + 213 + 1147 = 1451  missed in LC,      60 + 200 = 260   overcall in LC
sort(round(100 * table(gmat)/length(gmat),2)) 
##  11_01  00_11  01_00  01_11  00_01  11_NA  01_NA  11_11  00_NA  01_01  00_00 
##   0.02   0.03   0.06   0.07   0.36   1.00   1.48   5.91   6.04  10.31  74.71         => 0.03 + 0.07 + 0.36 = 0.46% missed in LC,   0.02 + 0.06 = 0.08% overcall in LC

## save genotypes with position
flcgeno <- apply(flcgeno, 2, function(v) gsub("00", 0, v))
flcgeno <- apply(flcgeno, 2, function(v) gsub("01", 1, v))
flcgeno <- apply(flcgeno, 2, function(v) gsub("11", 2, v))
fexgeno <- apply(fexgeno, 2, function(v) gsub("00", 0, v))
fexgeno <- apply(fexgeno, 2, function(v) gsub("01", 1, v))
fexgeno <- apply(fexgeno, 2, function(v) gsub("11", 2, v))
flcdat <- data.frame(lcgeno[,1], flcgeno)
fexdat <- data.frame(lcgeno[,1], fexgeno)
fdat <- data.frame(lcgeno[,1], gmat)
names(fdat) <- names(flcdat) <- names(fexdat) <- fgeno[1,]
#write.table(flcdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_gts_replication.txt", quote=F, row.names=F, sep="\t")        ## 3903 sites x 82 samples
#write.table(fexdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_exome_gts_replication.txt", quote=F, row.names=F, sep="\t")     ## 3903 sites x 82 samples
#write.table(fdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_gts_replication.txt", row.names=F, quote=F, sep="\t")    ## 3903 sites x 82 samples

## compare genotypes
gtc <- matrix(NA, nrow=nrow(flcgeno), ncol=ncol(flcgeno))
for (i in 1:ncol(flcgeno))  {
  gtc[,i] <- flcgeno[,i]==fexgeno[,i]
}
## save genotypes with position
fdat <- data.frame(lcgeno[,1], gtc)
names(fdat) <- fgeno[1,]
#write.table(fdat, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", row.names=F, quote=F, sep="\t")    ##  3903 sites x 81 samples

table(gtc)
##  FALSE   TRUE 
##   1715 205624     ## 99.2% concordance (discovery)     =>  0.8% discordance
##   1711 287497     ## 99.4% concordance (replication)   =>  0.6% discordance

## number and percentage of concordant sites and samples
csites <- apply(gtc, 1, sum, na.rm=T)
csamples <- apply(gtc, 2, sum, na.rm=T)
msites <- 100*apply(gtc, 1, mean, na.rm=T)
msamples <- 100*apply(gtc, 2, mean, na.rm=T)
round(summary(csites))
summary(csamples)
summary(msites)
summary(msamples)
round(summary(100-msamples),2)
## Discovery
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##       1      61      62      60      62      62 
##    2741    3281    3342    3317    3367    3397 
##      17.70   98.40  100.00   99.13  100.00  100.00     ## sites
##      81.70   99.40   99.50   99.18   99.60   99.80     ## samples
## Replication
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##       0      79      81      74      81      81 
##    3262    3504    3559    3549    3613    3756 
##       0.00  100.00  100.00   99.21  100.00  100.00     ## sites 
##      98.90   99.30   99.40   99.41   99.50   99.80     ## samples 
##       0.20    0.50    0.60    0.59    0.70    1.10     ## discordance by samples
sum(msites<90)      ##  48 (= 1.2%) after removing NAs
sum(msamples<90)    ##   0

## number and percentage of N/A genotypes
nacsites <- apply(gtc, 1, function(v) sum(is.na(v)))
nacsamples <- apply(gtc, 2, function(v) sum(is.na(v)))
namsites <- round(100*apply(gtc, 1, function(v) mean(is.na(v))),1)
namsamples <- round(100*apply(gtc, 2, function(v) mean(is.na(v))),1)
round(summary(nacsites))
round(summary(nacsamples))
round(summary(namsites),1)
round(summary(namsamples),1)
round(summary(namsamples),1)
## Dicoverys
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0.0     0.0     0.0     2.6     0.0    98.4     ## sites
##       0.4     1.5     2.1     2.6     3.6     6.3     ## samples
## Replication
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0.0     0.0     0.0     8.5     2.5    98.8     ## sites
##       3.3     6.9     8.2     8.5     9.6    15.7     ## samples
round(summary(100-msamples),2)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    0.25    0.49    0.58    0.59    0.70    1.05 

#pdf("~/uk10k/plots/discordance_lc_exome_by_site_chr20_replication.pdf", width=4, height=4, pointsize=10)
plot(sort(100-msites), col="red", pch=15, ylim=c(0,100), main="Discordance by Site", xlab="sites", ylab="discordance (%)", cex.lab=1.3, cex.main=1.3)
dev.off()
#pdf("~/uk10k/plots/discordance_lc_exome_by_sample_chr20_replication.pdf", width=4, height=4, pointsize=10)
plot(sort(100-msamples), col="red", pch=15, ylim=c(0,2), main="Discordance by Sample", xlab="samples", ylab="discordance (%)", cex.lab=1.3, cex.main=1.3)
abline(h=1.1, lty=2, col="darkgrey")
text(20, 1.2, "discordance = 1.1%", col="darkgrey")
dev.off()
#pdf("~/uk10k/plots/discordance_lc_exome_by_sample_chr20_replication_v2.pdf", width=4, height=4, pointsize=10)
plot(sort(100-msamples), col="red", pch=15, main="Discordance by Sample", xlab="samples", ylab="discordance (%)", cex.lab=1.3, cex.main=1.3)
dev.off()
#pdf("~/uk10k/plots/na_lc_exome_by_site_chr20_replication.pdf", width=4, height=4, pointsize=10)
plot(sort(namsites), col="blue", pch=15, main="N/A by Site", xlab="sites", ylab="N/A (%)", cex.lab=1.3, cex.main=1.3)
dev.off()
#pdf("~/uk10k/plots/na_lc_exome_by_sample_chr20_replication.pdf", width=4, height=4, pointsize=10)
plot(sort(namsamples), col="blue", pch=15, main="N/A by Sample", xlab="samples", ylab="N/A (%)", ylim=c(0,16), cex.lab=1.3, cex.main=1.3)
dev.off()

## samples with high NAs
fgeno[1,-1][which(namsamples>3)]




##### Prepare data for investigating attributes by concordance #####

## bsub -o info.out "perl ~/scripts/uk10k/get_vcf_chr_pos_info.pl /lustre/scratch106/projects/uk10k/RELEASE/UK10K_COHORT/REL-2011-12-01/v1/20.vqsr.anno.csq.sites.vcf.gz AF1,DP,HWE,VQSLOD > /lustre/scratch107/projects/uk10k/cohorts/REL-2011-12-01/vqsr/chr_pos_af_dp_hwe_vqslod_chr20.txt"

## get VQSR scores
#lc.sites <- read.delim("/lustre/scratch106/projects/uk10k/RELEASE/UK10K_COHORT/REL-2011-12-01/v1/20.vqsr.anno.csq.sites.vcf.gz", skip=48, stringsAsFactors=F)   ## 1031589  
lc.sites <- read.delim("/lustre/scratch107/projects/uk10k/cohorts/REL-2011-12-01/vqsr/chr_pos_af_dp_hwe_vqslod_chr20.txt", stringsAsFactors=F)
lc <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples_sites.vcf.gz", skip=56, stringsAsFactors=F)                        ## 5106 x 152
## remove INDELs from LC
indk1 <- nchar(lc$REF)==1 & nchar(lc$ALT)==1
## filter for VQSR [SNPs:-3.1951, truth sensitivity 99.85; indels:1.6334, truth sensitivity 95.0]
indk2 <- lc.sites$VQSLOD>=-3.1951 & nchar(lc.sites$REF)==1 & nchar(lc.sites$ALT)==1
mlc <- merge(lc[indk1,1:8], lc.sites[indk2,], by.x=2, by.y=2)   ##  5058
#write.table(mlc, "/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_info_selected.txt", row.names=F, quote=F, sep="\t")    ## 5058



##### Investigate attributes by concordance #####

info <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_info_selected.txt", stringsAsFactors=F)                        ## 5058
lc <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_uk10k_exome_samples_sites.vcf.gz", skip=56, stringsAsFactors=F)  ## 5106 x 152
f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", stringsAsFactors=F)                ## 3433 x 63
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", stringsAsFactors=F)              ## 3903 x 82
gts1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_gts_discovery.txt", stringsAsFactors=F)          ## 3433 x 62
gts2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_gts_replication.txt", stringsAsFactors=F)        ## 3903 x 82
lc1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_gts_discovery.txt", stringsAsFactors=F)                 ## 3433 x 62
lc2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_gts_replication.txt", stringsAsFactors=F)               ## 3903 x 82
ex1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_exome_gts_discovery.txt", stringsAsFactors=F)              ## 3433 x 62
ex2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_exome_gts_replication.txt", stringsAsFactors=F)            ## 3903 x 82
## remove sample QTL218676
indrm <- names(f1)=="QTL218676"
f1 <- f1[,!indrm]
### remove AC==0 ###
sac1 <- apply(mlc1[,-c(1:18)], 1, sum)
sac2 <- apply(mlc2[,-c(1:18)], 1, sum)
sum(sac1==0)     ## 1205 = 35%
sum(sac2==0)     ## 1216 = 31%
## f1 <- f1[sac1!=0,]
## f2 <- f2[sac2!=0,]
## gts1 <- gts1[sac1!=0,]
## gts2 <- gts2[sac2!=0,]
## lc1 <- lc1[sac1!=0,]
## lc2 <- lc2[sac2!=0,]
## ex1 <- ex1[sac1!=0,]
## ex2 <- ex2[sac2!=0,]
## remove INDELs from LC
indk <- nchar(lc$REF)==1 & nchar(lc$ALT)==1
mf1 <- merge(lc[indk,c(2,8)], f1, by.x=1, by.y=1)   ## 3433
mf2 <- merge(lc[indk,c(2,8)], f2, by.x=1, by.y=1)   ## 3903
## merge with INFO
mf1 <- merge(info, f1, by.x=1, by.y=1)              ## 2216
mf2 <- merge(info, f2, by.x=1, by.y=1)              ## 2670
mgts1 <- merge(info, gts1, by.x=1, by.y=1)          ## 2216
mgts2 <- merge(info, gts2, by.x=1, by.y=1)          ## 2670
mlc1 <- merge(info, lc1, by.x=1, by.y=1)            ## 2216
mlc2 <- merge(info, lc2, by.x=1, by.y=1)            ## 2670
mex1 <- merge(info, ex1, by.x=1, by.y=1)            ## 2216
mex2 <- merge(info, ex2, by.x=1, by.y=1)            ## 2670
## get concordance by sites
msites1 <- 100*apply(mf1[,19:79], 1, mean, na.rm=T)
msites2 <- 100*apply(mf2[,19:99], 1, mean, na.rm=T)
summary(msites1)
summary(msites2)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   18.03  100.00  100.00   99.42  100.00  100.00 
##    0.00  100.00  100.00   99.22  100.00  100.00
mean(msites1==100)  ##    84% have full concordance
mean(msites2==100)  ##    82% have full concordance
mean(msites1<99)    ## 348 = 16%
mean(msites2<99)    ## 475 = 18%
mean(msites1<98)    ## 123 =  5.6%%
mean(msites2<98)    ## 167 =  6.3%
mean(msites1<97)    ## 108 =  4.9%
mean(msites2<97)    ## 106 =  4.0%
mean(msites1<95)    ##  42 =  1.9%
mean(msites2<95)    ##  59 =  2.2%
mean(msites1<90)    ##  17 =  0.8%
mean(msites2<90)    ##  31 =  1.2%
## concordance for AC==0 in LC
masites1 <- 100*apply(mf1[ac1==0,19:79], 1, mean, na.rm=T)
masites2 <- 100*apply(mf2[ac2==0,19:99], 1, mean, na.rm=T)
summary(masites1)
summary(masites2)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   36.07  100.00  100.00   99.47  100.00  100.00 
##   22.73  100.00  100.00   99.21  100.00  100.00 

s <- seq(90,100,0.1)
disc1 <- 100*sapply(s, function(v) mean(msites1<v))
disc2 <- 100*sapply(s, function(v) mean(msites2<v))

#pdf("~/uk10k/plots/sites_min_concordance_lc_exome_chr20.pdf", width=8, height=4, pointsize=10)
par(mfrow=c(1,2))
plot(s, disc1, col="red", pch=16, xlab="minimum concordance (%)", ylab="sites (%)", xlim=c(90,100.3), main="Discovery", cex.main=1.3, cex.lab=1.2)
abline(v=s[84], lty=2, col="darkgrey")
abline(h=disc1[84], lty=2, col="darkgrey")
text(s[84]-5, disc1[84]+1, paste(round(disc1[84],1),"%"), col="darkgrey")
text(s[84]+1, disc1[84]-3.5, paste(s[84],"%"), col="darkgrey")
plot(s, disc2, col="red", pch=16, xlab="minimum concordance (%)", ylab="sites (%)", xlim=c(90,100.3), main="Replication", cex.main=1.3, cex.lab=1.2)
abline(v=s[88], lty=2, col="darkgrey")
abline(h=disc2[88], lty=2, col="darkgrey")
text(s[88]-5, disc2[88]+1, paste(round(disc2[88],1),"%"), col="darkgrey")
text(s[88]+1, disc2[88]-3.5, paste(s[88],"%"), col="darkgrey")
par(mfrow=c(1,1))
dev.off()

## genotype composition for low concordance
indc1 <- msites1<s[84]     ## 6.4% of sites
indc2 <- msites2<s[88]     ## 8.3% of sites

table(as.matrix(mgts1[,-c(1:18)]))
table(as.matrix(mgts1[indc1,-c(1:18)]))
##  00_00  00_01  00_11  00_NA  01_00  01_01  01_11  01_NA  11_01  11_11  11_NA 
## 166908    649     25   3894    151  22910    196    881     68  12407    592 
##   7705    410     25   1082    104   2063    172    353     52   1016    316
table(as.matrix(mgts2[,-c(1:18)]))
table(as.matrix(mgts2[indc2,-c(1:18)]))
##  00_00  00_01  00_11  00_NA  01_00  01_01  01_11  01_NA  11_01  11_11  11_NA 
## 236108   1144     91  19038    200  32605    213   4668     58  17930   2792 
##  13169    879     91   4415    148   3475    180   1146     44   1895    640
dt1 <- round(100 * table(as.matrix(mgts1[,-c(1:18)]))/length(as.matrix(mgts1[,-c(1:18)])),2)
dt2 <- round(100 * table(as.matrix(mgts1[indc1,-c(1:18)]))/length(as.matrix(mgts1[indc1,-c(1:18)])),2)
##  00_00 00_01 00_11 00_NA 01_00 01_01 01_11 01_NA 11_01 11_11 11_NA 
##  79.98  0.31  0.01  1.87  0.07 10.98  0.09  0.42  0.03  5.95  0.28 
##  57.94  3.08  0.19  8.14  0.78 15.51  1.29  2.65  0.39  7.64  2.38    ## fewer RR calls, more RA and AA calls; 10x as many 00_01, 00_11, 01_00, 01_11, 5x as many NA calls
rt1 <- round(100 * table(as.matrix(mgts2[,-c(1:18)]))/length(as.matrix(mgts2[,-c(1:18)])),2)
rt2 <- round(100 * table(as.matrix(mgts2[indc2,-c(1:18)]))/length(as.matrix(mgts2[indc2,-c(1:18)])),2)
##  00_00 00_01 00_11 00_NA 01_00 01_01 01_11 01_NA 11_01 11_11 11_NA 
##  74.99  0.36  0.03  6.05  0.06 10.36  0.07  1.48  0.02  5.69  0.89 
##  50.49  3.37  0.35 16.93  0.57 13.32  0.69  4.39  0.17  7.27  2.45    ## fewer RR calls, more RA and AA calls; 10x as many 00_01, 00_11, 01_00, 01_11, 5x as many NA calls 

#pdf("~/uk10k/plots/gts_lc_exome_distribution_chr20_discovery.pdf", width=11, height=4, pointsize=14)
par(mfrow=c(1,3))
barplot(rbind(dt1[c(1,6,10)], dt2[c(1,6,10)]), beside=T, col=c("darkred","red"), las=1,
        ylab="percentage of genotype combination", xlab="genotype combination", main="Concordance",
        cex.lab=1.3, cex.axis=1.3, cex.main=1.5, cex=1.1)
legend(3.5,80, legend=c("all", "discordant sites (<98.3%)"), fill=c("darkred","red"))
barplot(rbind(dt1[c(2,3,5,7,9)], dt2[c(2,3,5,7,9)]), beside=T, col=c("blue","lightblue"), las=1,
        ylab="percentage of genotype combination", xlab="genotype combination", main="Discordance", cex.lab=1.3, cex.axis=1.3, cex.main=1.5, cex=1.1)
legend(5,3, legend=c("all", "discordant sites (<98.3%)"), fill=c("blue","lightblue"))
barplot(rbind(dt1[c(4,8,11)], dt2[c(4,8,11)]), beside=T, col=c("darkgreen","lightgreen"), las=1,
        ylab="percentage of genotype combination", xlab="genotype combination", main="Missing genotypes in exomes", cex.lab=1.3, cex.axis=1.5, cex.main=1.5, cex=1.1)
legend(3.5,8, legend=c("all", "discordant sites (<98.3%)"), fill=c("darkgreen","lightgreen"))
par(mfrow=c(1,1))
dev.off()

#pdf("~/uk10k/plots/gts_lc_exome_distribution_chr20_replication.pdf", width=11, height=4, pointsize=14)
par(mfrow=c(1,3))
barplot(rbind(rt1[c(1,6,10)], rt2[c(1,6,10)]), beside=T, col=c("darkred","red"), las=1,
        ylab="percentage of genotype combination", xlab="genotype combination", main="Concordance",
        cex.lab=1.3, cex.axis=1.3, cex.main=1.5, cex=1.1)
legend(3.5,70, legend=c("all", "discordant sites (<98.7%)"), fill=c("darkred","red"))
barplot(rbind(rt1[c(2,3,5,7,9)], rt2[c(2,3,5,7,9)]), beside=T, col=c("blue","lightblue"), las=1,
        ylab="percentage of genotype combination", xlab="genotype combination", main="Discordance", cex.lab=1.3, cex.axis=1.3, cex.main=1.5, cex=1.1)
legend(5,3, legend=c("all", "discordant sites (<98.7%)"), fill=c("blue","lightblue"))
barplot(rbind(rt1[c(4,8,11)], rt2[c(4,8,11)]), beside=T, col=c("darkgreen","lightgreen"), las=1,
        ylab="percentage of genotype combination", xlab="genotype combination", main="Missing genotypes in exomes", cex.lab=1.3, cex.axis=1.5, cex.main=1.5, cex=1.1)
legend(3.5,15, legend=c("all", "discordant sites (<98.7%)"), fill=c("darkgreen","lightgreen"))
par(mfrow=c(1,1))
dev.off()


### INFO by concordance ###  ##           removing AC==0 in LC
ind1.99 <- msites1<s[84]     ## 218       46 
ind2.99 <- msites2<s[88]     ## 322       61
ind1.95 <- msites1<95        ##  59       12
ind2.95 <- msites2<95        ##  84       20
ind1.90 <- msites1<90        ##  25        4
ind2.90 <- msites2<90        ##  46       15

## for the plots
mycols <- c("darkgrey","blue","darkred","red","orange")
mytype <- c(1,3,2,4,5)

### VQSR scores ###
summary(mf1$VQSLOD)
summary(mf2$VQSLOD)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -3.126   4.010   5.920   5.673   7.724  11.600 
##  -3.194   4.342   6.157   5.900   7.850  11.600 

#pdf("~/uk10k/plots/vqsr_by_concordance_lc_exome_chr20_discovery.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(mf1$VQSLOD, from=min(mf1$VQSLOD), to=max(mf1$VQSLOD)), col=mycols[1], lty=mytype[1], lwd=1.6, main="VQSLOD", xlim=c(-5,13), ylim=c(0,0.21), cex.main=1.3)
lines(density(mf1$VQSLOD[!ind1.99], from=min(mf1$VQSLOD), to=max(mf1$VQSLOD)), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(mf1$VQSLOD[ind1.99], from=min(mf1$VQSLOD), to=max(mf1$VQSLOD)), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(mf1$VQSLOD[ind1.95], from=min(mf1$VQSLOD), to=max(mf1$VQSLOD)), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(mf1$VQSLOD[ind1.90], from=min(mf1$VQSLOD), to=max(mf1$VQSLOD)), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(-5, 0.21, legend=c("all (n=3421)", ">=98 (n=3244)", "<98 (n=544)", "<95 (n=59)", "<90 (n=25)"), col=mycols, lty=mytype, lwd=1.6, cex=0.8)
dev.off()
#pdf("~/uk10k/plots/vqsr_by_concordance_lc_exome_chr20_replication.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(mf2$VQSLOD, from=min(mf2$VQSLOD), to=max(mf2$VQSLOD)), col=mycols[1], lty=mytype[1], lwd=1.6, main="VQSLOD", xlim=c(-5,13), ylim=c(0,0.21), cex.main=1.3)
lines(density(mf2$VQSLOD[!ind2.99], from=min(mf2$VQSLOD), to=max(mf2$VQSLOD)), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(mf2$VQSLOD[ind2.99], from=min(mf2$VQSLOD), to=max(mf2$VQSLOD)), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(mf2$VQSLOD[ind2.95], from=min(mf2$VQSLOD), to=max(mf2$VQSLOD)), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(mf2$VQSLOD[ind2.90], from=min(mf2$VQSLOD), to=max(mf2$VQSLOD)), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(-5, 0.21, legend=c("all (n=3887)", ">=98 (n=3565)", "<98 (n=686)", "<95 (n=84)", "<90 (n=46)"), col=mycols, lty=mytype, lwd=1.6, cex=0.8)
dev.off()
#pdf("~/uk10k/plots/vqsr_by_concordance_lc_exome_chr20_discovery_v2.pdf", width=3.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/vqsr_by_concordance_lc_exome_chr20_discovery_v2.png", width=720, height=720, pointsize=25)
plot(mf1$VQSLOD, msites1, xlab="VQSLOD score", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
dev.off()
#pdf("~/uk10k/plots/vqsr_by_concordance_lc_exome_chr20_replication_v2.pdf", width=3.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/vqsr_by_concordance_lc_exome_chr20_replication_v2.png", width=720, height=720, pointsize=25)
plot(mf2$VQSLOD, msites2, xlab="VQSLOD score", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
dev.off()

### AC and MAF recalculated from AC for the samples used - highly correlated with AC from INFO and with MAF EUR ###
ac1 <- get.vcf.info(mf1$INFO, "AC")
ac2 <- get.vcf.info(mf2$INFO, "AC")
sac1 <- apply(mlc1[,-c(1:18)], 1, sum)
sac2 <- apply(mlc2[,-c(1:18)], 1, sum)
sum(sac1==0)     ## 1205 = 35%
sum(sac2==0)     ## 1216 = 31%
n1 <- ncol(mlc1[,-c(1:18)])
n2 <- ncol(mlc2[,-c(1:18)])
maf1 <- sac1 / (2*n1)
maf2 <- sac2 / (2*n2)
maf1[maf1>0.5] <- 1-maf1[maf1>0.5]
maf2[maf2>0.5] <- 1-maf2[maf2>0.5]
round(summary(maf1),2)
round(summary(maf2),2)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    0.00    0.00    0.01    0.08    0.11    0.50 
##    0.00    0.00    0.01    0.08    0.12    0.50 

#pdf("~/uk10k/plots/ac_by_concordance_lc_exome_chr20_discovery.pdf", width=6.5, height=4.5, pointsize=10)
boxplot(ac1, ac1[!ind1.99], ac1[ind1.99], ac1[ind1.95], ac1[ind1.90], col=mycols, main="AC",
        names=c("all (n=3421)", ">=98 (n=3244)", "<98 (n=544)", "<95 (n=59)", "<90 (n=25)"), xlab="concordance", ylab="AC distribution", cex.lab=1.3, cex.main=1.5)
dev.off()
#pdf("~/uk10k/plots/ac_by_concordance_lc_exome_chr20_replication.pdf", width=6.5, height=4.5, pointsize=10)
mycols <- c("darkgrey","blue","darkred","red","orange")
boxplot(ac2, ac2[!ind2.99], ac2[ind2.99], ac2[ind2.95], ac2[ind2.90], col=mycols, main="AC",
        names=c("all (n=3887)", ">=98 (n=3565)", "<98 (n=686)", "<95 (n=84)", "<90 (n=46)"), xlab="concordance", ylab="AC distribution", cex.lab=1.3, cex.main=1.5)
dev.off()

ind1 <- cut(maf1, breaks=c(0, 0.01, 0.05, 0.1, 0.5), include.lowest=T, labels=F)
ind2 <- cut(maf2, breaks=c(0, 0.01, 0.05, 0.1, 0.5), include.lowest=T, labels=F)
table(ind1)
table(ind2)
##  [0,0.01] (0.01,0.05]  (0.05,0.1]   (0.1,0.5] 
##         1           2           3           4 
##      1758         521         264         878 
##      1888         680         276        1043 

#pdf("~/uk10k/plots/concordance_lc_exome_chr20_discovery.pdf", width=7.5, height=3.5, pointsize=9)
par(mfrow=c(1,2))
hist(msites1, breaks=50, xlab="concordance (%)", cex.lab=1.3, cex.axis=1.1, main="")
plot(density(msites1, to=100), col=mycols[1], cex.axis=1.1, las=1, main="", xlab="concordance (%)", cex.lab=1.3)
lines(density(msites1[ind1==1], to=100, bw=1), col=mycols[2])
lines(density(msites1[ind1==2], to=100, bw=1), col=mycols[3])
lines(density(msites1[ind1==3], to=100, bw=1), col=mycols[4])
lines(density(msites1[ind1==4], to=100, bw=1), col=mycols[5])
legend(20,0.55, legend=c("all","MAF=[0,0.01]","MAF=(0.01,0.05]","MAF=(0.05,0.1]","MAF=(0.1,0.5]"), fill=mycols)
par(mfrow=c(1,1))
dev.off()
#pdf("~/uk10k/plots/concordance_lc_exome_chr20_replication.pdf", width=7.5, height=3.5, pointsize=9)
par(mfrow=c(1,2))
hist(msites2, breaks=50, xlab="concordance (%)", cex.lab=1.3, cex.axis=1.1, main="")
plot(density(msites2, to=100), col=mycols[1], cex.axis=1.1, las=1, main="", xlab="concordance (%)", cex.lab=1.3)
lines(density(msites2[ind2==1], to=100, bw=1), col=mycols[2])
lines(density(msites2[ind2==2], to=100, bw=1), col=mycols[3])
lines(density(msites2[ind2==3], to=100, bw=1), col=mycols[4])
lines(density(msites2[ind2==4], to=100, bw=1), col=mycols[5])
legend(0,0.4, legend=c("all","MAF=[0,0.01]","MAF=(0.01,0.05]","MAF=(0.05,0.1]","MAF=(0.1,0.5]"), fill=mycols)
par(mfrow=c(1,1))
dev.off()
#pdf("~/uk10k/plots/maf_by_concordance_lc_exome_chr20_discovery.pdf", width=7.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/maf_by_concordance_lc_exome_chr20_discovery.png", width=1500, height=700, pointsize=23)
#par(mfrow=c(1,2))
par(mfrow=c(1,2), mai=c(1.8,1.8,0.5,0.1))
plot(maf1, msites1, xlab="MAF", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
boxplot(msites1, msites1[ind1==1], msites1[ind1==2], msites1[ind1==3], msites1[ind1==4], 
        names=c("all","[0,0.01]","(0.01,0.05]","(0.05,0.1]","(0.1,0.5]"), cex.lab=1.3, las=1, xlab="MAF", ylab="concordance (%)", cex.axis=0.8)
par(mfrow=c(1,1))
dev.off()
#pdf("~/uk10k/plots/maf_by_concordance_lc_exome_chr20_replication.pdf", width=7.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/maf_by_concordance_lc_exome_chr20_replication.png", width=1500, height=700, pointsize=23)
#par(mfrow=c(1,2))
par(mfrow=c(1,2), mai=c(1.8,1.8,0.5,0.1))
plot(maf2, msites2, xlab="MAF", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
boxplot(msites2, msites2[ind2==1], msites2[ind2==2], msites2[ind2==3], msites2[ind2==4], 
        names=c("all","[0,0.01]","(0.01,0.05]","(0.05,0.1]","(0.1,0.5]"), cex.lab=1.3, las=1, xlab="MAF", ylab="concordance (%)", cex.axis=0.8)
par(mfrow=c(1,1))
dev.off()

### AF_MAX ###
af.max1 <- get.vcf.info(mf1$INFO, "AF_MAX")
af.max2 <- get.vcf.info(mf2$INFO, "AF_MAX")
#pdf("~/uk10k/plots/af_max_by_concordance_lc_exome_chr20_discovery.pdf", width=6.5, height=4.5, pointsize=10)
boxplot(af.max1, af.max1[!ind1.99], af.max1[ind1.99], af.max1[ind1.95], af.max1[ind1.90], col=mycols, main="AF MAX",
        names=c("all (n=3421)", ">=98 (n=3244)", "<98 (n=544)", "<95 (n=59)", "<90 (n=25)"), xlab="concordance", ylab="AF MAX distribution", cex.lab=1.3, cex.main=1.5)
dev.off()
#pdf("~/uk10k/plots/af_max_by_concordance_lc_exome_chr20_replication.pdf", width=6.5, height=4.5, pointsize=10)
boxplot(af.max2, af.max2[!ind2.99], af.max2[ind2.99], af.max2[ind2.95], af.max2[ind2.90], col=mycols, main="AF MAX",
        names=c("all (n=3887)", ">=98 (n=3565)", "<98 (n=686)", "<95 (n=84)", "<90 (n=46)"), xlab="concordance", ylab="AF MAX distribution", cex.lab=1.3, cex.main=1.5)
dev.off()

### DP ###
dp1 <- get.vcf.info(mf1$INFO, "DP")
dp2 <- get.vcf.info(mf2$INFO, "DP")
#pdf("~/uk10k/plots/depth_by_concordance_lc_exome_chr20_discovery.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(dp1, bw=900, from=0), col=mycols[1], lty=mytype[1], lwd=1.6, main="Depth", ylim=c(0,0.000115), xlim=c(0,37000), cex.main=1.5)
lines(density(dp1[!ind1.99], bw=900, from=0), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(dp1[ind1.99], bw=3000, from=0), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(dp1[ind1.95], bw=3000, from=0), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(dp1[ind1.90], bw=3000, from=0), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(22000, 0.000115, legend=c("all (n=3421)", ">=98 (n=3244)", "<98 (n=544)", "<95 (n=59)", "<90 (n=25)"), col=mycols, lty=mytype, lwd=1.6, cex=0.9)
dev.off()
#pdf("~/uk10k/plots/depth_by_concordance_lc_exome_chr20_replication.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(dp2, bw=900, from=0), col=mycols[1], lty=mytype[1], lwd=1.6, main="Depth", ylim=c(0,0.000115), xlim=c(0,37000), cex.main=1.5)
lines(density(dp2[!ind2.99], bw=900, from=0), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(dp2[ind2.99], bw=3000, from=0), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(dp2[ind2.95], bw=3000, from=0), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(dp2[ind2.90], bw=3000, from=0), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(22000, 0.000115, legend=c("all (n=3887)", ">=98 (n=3565)", "<98 (n=686)", "<95 (n=84)", "<90 (n=46)"), col=mycols, lty=mytype, lwd=1.6, cex=0.9)
dev.off()

### HWE ###
hwe1 <- get.vcf.info(mf1$INFO, "HWE")
hwe2 <- get.vcf.info(mf2$INFO, "HWE")
lhwe1 <- -log10(hwe1+0.0000001)
lhwe2 <- -log10(hwe2+0.0000001)
ind1 <- cut(hwe1, breaks=c(0, 0.01, 0.05, 0.1, 1), include.lowest=T, labels=F)
ind2 <- cut(hwe2, breaks=c(0, 0.01, 0.05, 0.1, 1), include.lowest=T, labels=F)
indl1 <- cut(lhwe1, breaks=c(min(lhwe1), 0.1, 0.5, 1, 2, 4, max(lhwe1)), include.lowest=T, labels=F)
indl2 <- cut(lhwe2, breaks=c(min(lhwe2), 0.1, 0.5, 1, 2, 4, max(lhwe2)), include.lowest=T, labels=F)
table(ind1)
table(ind2)
## [0,0.01] (0.01,0.05] (0.05,0.1] (0.1,1]
##        1          2           3       4 
##       46         60          71    3244 
##       60         76          92    3659 
table(indl1)
table(indl2)
##  [-4.34e-08,0.1] (0.1,0.5] (0.5,1] (1,2] (2,4] (4,7]
##                1        2        3     4     5     6 
##             2151      765      328   131    31    15 
##             2464      845      350   168    34    26 

#pdf("~/uk10k/plots/hwe_by_concordance_lc_exome_chr20_discovery.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(hwe1, from=0, to=1), col=mycols[1], lty=mytype[1], lwd=1.6, main="HWE p-value", cex.main=1.5, ylim=c(0,4.5))
lines(density(hwe1[!ind1.99], from=0, to=1), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(hwe1[ind1.99], from=0, to=1), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(hwe1[ind1.95], from=0, to=1), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(hwe1[ind1.90], from=0, to=1), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(0, 4.5, legend=c("all (n=3421)", ">=98 (n=3244)", "<98 (n=544)", "<95 (n=59)", "<90 (n=25)"), col=mycols, lty=mytype, lwd=1.6)
dev.off()
#pdf("~/uk10k/plots/hwe_by_concordance_lc_exome_chr20_replication.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(hwe2, from=0, to=1), col=mycols[1], lty=mytype[1], lwd=1.6, main="HWE p-value", cex.main=1.5, ylim=c(0,4.5))
lines(density(hwe2[!ind2.99], from=0, to=1), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(hwe2[ind2.99], from=0, to=1), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(hwe2[ind2.95], from=0, to=1), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(hwe2[ind2.90], from=0, to=1), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(0, 4.5, legend=c("all (n=3887)", ">=98 (n=3565)", "<98 (n=686)", "<95 (n=84)", "<90 (n=46)"), col=mycols, lty=mytype, lwd=1.6)
dev.off()
#pdf("~/uk10k/plots/hwe_by_concordance_lc_exome_chr20_discovery_v2.pdf", width=8.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/hwe_by_concordance_lc_exome_chr20_discovery_v2.png", width=1700, height=700, pointsize=25)
par(mfrow=c(1,2))
plot(lhwe1, msites1, xlab="-log10(HWE)", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
boxplot(msites1, msites1[indl1==1], msites1[indl1==2], msites1[indl1==3], msites1[indl1==4], msites1[indl1==5], msites1[indl1==6], 
        names=c("all","(0,0.1]","(0.1,0.5]","(0.5,1]","(1,2]","(2,4]","(4,7]"), cex.lab=1.3, las=1, xlab="-log10(HWE)", ylab="concordance (%)", cex.axis=0.8)
par(mfrow=c(1,1))
dev.off()
#pdf("~/uk10k/plots/hwe_by_concordance_lc_exome_chr20_replication_v2.pdf", width=8.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/hwe_by_concordance_lc_exome_chr20_replication_v2.png", width=1700, height=700, pointsize=25)
par(mfrow=c(1,2))
plot(lhwe2, msites2, xlab="-log10(HWE)", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
boxplot(msites2, msites2[indl2==1], msites2[indl2==2], msites2[indl2==3], msites2[indl2==4], msites2[indl2==5], msites2[indl2==6],  
        names=c("all","(0,0.1]","(0.1,0.5]","(0.5,1]","(1,2]","(2,4]","(4,7]"), cex.lab=1.3, las=1, xlab="-log10(HWE)", ylab="concordance (%)", cex.axis=0.8)
par(mfrow=c(1,1))
dev.off()

### ICF ###
icf1 <- get.vcf.info(mf1$INFO, "ICF")
icf2 <- get.vcf.info(mf2$INFO, "ICF")
#pdf("~/uk10k/plots/icf_by_concordance_lc_exome_chr20_discovery.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(icf1, from=-0.15, to=0.15, bw=0.005), col=mycols[1], lty=mytype[1], lwd=1.6, main="Inbreeding coefficient", cex.main=1.5)
lines(density(icf1[!ind1.99], from=-0.15, to=0.15, bw=0.005), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(icf1[ind1.99], from=-0.15, to=0.15, bw=0.01), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(icf1[ind1.95], from=-0.15, to=0.15, bw=0.01), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(icf1[ind1.90], from=-0.15, to=0.15, bw=0.01), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(-0.15, 44, legend=c("all (n=3421)", ">=98 (n=3244)", "<98 (n=544)", "<95 (n=59)", "<90 (n=25)"), col=mycols, lty=mytype, lwd=1.6, cex=0.9)
dev.off()
#pdf("~/uk10k/plots/icf_by_concordance_lc_exome_chr20_replication.pdf", width=4.5, height=4.5, pointsize=10)
plot(density(icf2, from=-0.15, to=0.15, bw=0.005), col=mycols[1], lty=mytype[1], lwd=1.6, main="Inbreeding coefficient", cex.main=1.5)
lines(density(icf2[!ind2.99], from=-0.15, to=0.15, bw=0.01), col=mycols[2], lty=mytype[2], lwd=1.6)
lines(density(icf2[ind2.99], from=-0.15, to=0.15, bw=0.01), col=mycols[3], lty=mytype[3], lwd=1.6)
lines(density(icf2[ind2.95], from=-0.15, to=0.15, bw=0.01), col=mycols[4], lty=mytype[4], lwd=1.6)
lines(density(icf2[ind2.90], from=-0.15, to=0.15, bw=0.01), col=mycols[5], lty=mytype[5], lwd=1.6)
legend(-0.15, 43, legend=c("all (n=3887)", ">=98 (n=3565)", "<98 (n=686)", "<95 (n=84)", "<90 (n=46)"), col=mycols, lty=mytype, lwd=1.6, cex=0.9)
dev.off()

### IMP2 ###
imp1 <- get.vcf.info(mf1$INFO, "IMP2", num=F)
imp2 <- get.vcf.info(mf2$INFO, "IMP2", num=F)
imp1 <- as.numeric(vecstr2matstr(imp1,",")[,2])
imp2 <- as.numeric(vecstr2matstr(imp2,",")[,2])
ind1 <- cut(imp1, breaks=c(0, 0.4, 1), include.lowest=T, labels=F)
ind2 <- cut(imp2, breaks=c(0, 0.4, 1), include.lowest=T, labels=F)
table(ind1)
table(ind2)
##    1    2 
##    1 3420 
##    9 3878 

#pdf("~/uk10k/plots/imp_by_concordance_lc_exome_chr20_discovery.pdf", width=6.5, height=4.5, pointsize=10)
boxplot(imp1, imp1[!ind1.99], imp1[ind1.99], imp1[ind1.95], imp1[ind1.90], col=mycols, main="IMPUTE2 info score",
        names=c("all (n=3421)", ">=98 (n=3244)", "<98 (n=544)", "<95 (n=59)", "<90 (n=25)"), xlab="concordance", ylab="info score distribution", cex.lab=1.3, cex.main=1.5)
dev.off()
#pdf("~/uk10k/plots/imp_by_concordance_lc_exome_chr20_replication.pdf", width=6.5, height=4.5, pointsize=10)
boxplot(imp2, imp2[!ind2.99], imp2[ind2.99], imp2[ind2.95], imp2[ind2.90], col=mycols, main="IMPUTE2 info score",
        names=c("all (n=3887)", ">=98 (n=3565)", "<98 (n=686)", "<95 (n=84)", "<90 (n=46)"), xlab="concordance", ylab="info score distribution", cex.lab=1.3, cex.main=1.5)
dev.off()
#pdf("~/uk10k/plots/imp_by_concordance_lc_exome_chr20_discovery_v2.pdf", width=3.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/imp_by_concordance_lc_exome_chr20_discovery_v2.png", width=700, height=700, pointsize=25)
plot(imp1, msites1, xlab="IMPUTE2 info score", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
dev.off()
#pdf("~/uk10k/plots/imp_by_concordance_lc_exome_chr20_replication_v2.pdf", width=3.5, height=3.5, pointsize=9)
#png("~/uk10k/plots/imp_by_concordance_lc_exome_chr20_replication_v2.png", width=700, height=700, pointsize=25)
plot(imp2, msites2, xlab="IMPUTE2 info score", ylab="concordance (%)", cex.lab=1.3, cex.axis=1.1, las=1)
dev.off()


### Plot concordance by chromosomal position ###
chrsize <- read.delim("~/data/genome_info/human_chr_size.txt", header=F, stringsAsFactors=F)
chr20size <- chrsize[20,2] / 1000000
dsites1 <- 100 - msites1
dsites2 <- 100 - msites2

#pdf("~/uk10k/plots/concordance_by_position_lc_exome_chr20.pdf", width=10, height=6, pointsize=10)
#png("~/uk10k/plots/concordance_by_position_lc_exome_chr20.png", width=1000, height=600, pointsize=18)
#par(mfrow=c(2,1), mai=c(0.7,0.7,0.6,0.1))
par(mfrow=c(2,1), mai=c(1.3,1.3,0.6,0.1))
plot(c(0,chr20size), c(0,100), ty="n", ylab="discordance (%)", main="Discovery", xlab="", cex.main=1.3, cex.lab=1.3)
points(mf1[,1]/1000000, dsites1, col="blue")
plot(c(0,chr20size), c(0,100), ty="n", xlab="position (chr20)", ylab="discordance (%)", main="Replication", cex.main=1.3, cex.lab=1.3)
points(mf2[,1]/1000000, dsites2, col="blue")
par(mfrow=c(1,1))
dev.off()

### Plot NAs by chromosomal position ###
chrsize <- read.delim("~/data/genome_info/human_chr_size.txt", header=F, stringsAsFactors=F)
chr20size <- chrsize[20,2] / 1000000
nasites1 <- round(100*apply(gts1, 1, function(v) mean(is.na(v))),1)
nasites2 <- round(100*apply(gts2, 1, function(v) mean(is.na(v))),1)
nasites1 <- round(100*apply(mf1[,19:79], 1, function(v) mean(is.na(v))),1)
nasites2 <- round(100*apply(mf2[,19:99], 1, function(v) mean(is.na(v))),1)

#pdf("~/uk10k/plots/na_by_position_lc_exome_chr20.pdf", width=10, height=6, pointsize=10)
#png("~/uk10k/plots/na_by_position_lc_exome_chr20.png", width=1000, height=600, pointsize=18)
#par(mfrow=c(2,1), mai=c(0.7,0.7,0.6,0.1))
par(mfrow=c(2,1), mai=c(1.3,1.3,0.6,0.1))
plot(c(0,chr20size), c(0,100), ty="n", ylab="discordance (%)", main="Discovery", xlab="", cex.main=1.3, cex.lab=1.3)
points(mf1[,1]/1000000, nasites1, col="blue")
plot(c(0,chr20size), c(0,100), ty="n", xlab="position (chr20)", ylab="discordance (%)", main="Replication", cex.main=1.3, cex.lab=1.3)
points(mf2[,1]/1000000, nasites2, col="blue")
par(mfrow=c(1,1))
dev.off()

### Overlap between discovery and replication
length(intersect(f1$POS, f2$POS))     ## 2259 = 66%,  2259/nrow(f2) = 58%
length(setdiff(f1$POS, f2$POS))       ## 1174 = 34%
length(setdiff(f2$POS, f1$POS))       ## 1644 = 42%

### Discordance by overlap of positions
indb <- intersect(f1$POS, f2$POS)
indo1 <- sapply(indb, function(v) which(v==f1$POS))
indo2 <- sapply(indb, function(v) which(v==f2$POS))
mdat1 <- f1[indo1,]
mdat2 <- f2[indo2,]
ddat1 <- f1[-indo1,]
ddat2 <- f2[-indo2,]
msites1 <- round(100*apply(mdat1[,-1], 1, mean, na.rm=T),1)
msites2 <- round(100*apply(mdat2[,-1], 1, mean, na.rm=T),1)
dsites1 <- round(100*apply(ddat1[,-1], 1, mean, na.rm=T),1)
dsites2 <- round(100*apply(ddat2[,-1], 1, mean, na.rm=T),1)
msamples1 <- round(100*apply(mdat1[,-1], 2, mean, na.rm=T),1)
msamples2 <- round(100*apply(mdat2[,-1], 2, mean, na.rm=T),1)
dsamples1 <- round(100*apply(ddat1[,-1], 2, mean, na.rm=T),1)
dsamples2 <- round(100*apply(ddat2[,-1], 2, mean, na.rm=T),1)

boxplot(msites1, msites2, dsites1, dsites2)
plot(msites1, msites2)
boxplot(msamples1, msamples2, dsamples1, dsamples2)


##### Compare with GWAS genotypes #####
info <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_info_selected.txt", stringsAsFactors=F)                        ## 5058
gwas <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/data_2453-samples/chr20/gwas_common_sites_samples_uk10k_ord.txt", stringsAsFactors=F)    ## 20927 x 1596
lc1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_gts_discovery.txt", stringsAsFactors=F)                 ## 3433 x 62
lc2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_gts_replication.txt", stringsAsFactors=F)               ## 3903 x 82
ex1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_exome_gts_discovery.txt", stringsAsFactors=F)              ## 3433 x 62
ex2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_exome_gts_replication.txt", stringsAsFactors=F)            ## 3903 x 82
f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", stringsAsFactors=F)                ## 3433 x 63
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", stringsAsFactors=F)              ## 3903 x 82
## remove sample QTL218676
indrm <- names(f1)=="QTL218676"
f1 <- f1[,!indrm]
## select samples
indc1 <- sapply(names(f1)[-1], function(v) which(v==names(gwas)))
indc2 <- sapply(names(f2)[-1], function(v) which(v==names(gwas)))
gwas1 <- gwas[,c(2,4,5,indc1)]
gwas2 <- gwas[,c(2,4,5,indc2)]
## select sites
mgwas1 <- merge(f1[,1], gwas1, by.x=1, by.y=1)     ## 603 x 64
mgwas2 <- merge(f2[,1], gwas2, by.x=1, by.y=1)     ## 684 x 84
mlc1 <- merge(mgwas1[,1], lc1, by.x=1, by.y=1)     ## 603 x 62
mlc2 <- merge(mgwas2[,1], lc2, by.x=1, by.y=1)     ## 684 x 82 
mex1 <- merge(mgwas1[,1], ex1, by.x=1, by.y=1)     ## 603 x 62
mex2 <- merge(mgwas2[,1], ex2, by.x=1, by.y=1)     ## 684 x 82 
## edit GWAS data
fgwas1 <- apply(mgwas1[,-c(1:3)], 2, function(v) recode.gwas.gts(ref=mgwas1[,2], alt=mgwas1[,3], v))   ## 603 x 61
fgwas2 <- apply(mgwas2[,-c(1:3)], 2, function(v) recode.gwas.gts(ref=mgwas2[,2], alt=mgwas2[,3], v))   ## 684 x 81
## compare GWAS with LC and exome
lgts1 <- matrix(NA, nrow=nrow(fgwas1), ncol=ncol(fgwas1))
lgts2 <- matrix(NA, nrow=nrow(fgwas2), ncol=ncol(fgwas2))
egts1 <- matrix(NA, nrow=nrow(fgwas1), ncol=ncol(fgwas1))
egts2 <- matrix(NA, nrow=nrow(fgwas2), ncol=ncol(fgwas2))
for (j in 1:ncol(fgwas1)) {
  lgts1[,j] <- fgwas1[,j]==mlc1[,-1][,j]
}
for (j in 1:ncol(fgwas2)) {
  lgts2[,j] <- fgwas2[,j]==mlc2[,-1][,j]
}
for (j in 1:ncol(fgwas1)) {
  egts1[,j] <- fgwas1[,j]==mex1[,-1][,j]
}
for (j in 1:ncol(fgwas2)) {
  egts2[,j] <- fgwas2[,j]==mex2[,-1][,j]
}
## LC with GWAS
round(100 * mean(lgts1, na.rm=T),2)    ##  99.76
round(100 * mean(lgts2, na.rm=T),2)    ##  99.75  
round(100 * mean(c(as.vector(lgts1), as.vector(lgts2)), na.rm=T),2)    ## 99.75
## Exome with GWAS
round(100 * mean(egts1, na.rm=T),2)    ##  99.97
round(100 * mean(egts2, na.rm=T),2)    ##  99.97     ## concordance of GWAS with exomes is slightly better than with LC  
round(100 * mean(c(as.vector(egts1), as.vector(egts2)), na.rm=T),2)    ## 99.97
## LC with exomes over the same ~600 sites
tgts1 <- matrix(NA, nrow=nrow(fgwas1), ncol=ncol(fgwas1))
tgts2 <- matrix(NA, nrow=nrow(fgwas2), ncol=ncol(fgwas2))
for (j in 1:ncol(fgwas1)) {
  tgts1[,j] <- mex1[,j]==mlc1[,j]
}
for (j in 1:ncol(fgwas2)) {
  tgts2[,j] <- mex2[,j]==mlc2[,j]
}
round(100 * mean(tgts1, na.rm=T),2)     ## 99.67
round(100 * mean(tgts2, na.rm=T),2)     ## 99.64
round(100 * mean(c(as.vector(tgts1), as.vector(tgts2)), na.rm=T),2)    ## 99.65
## MAF for ~600common sites
lc1.info <- merge(mlc1, info, by.x=1, by.y=1)   ## 603
lc2.info <- merge(mlc2, info, by.x=1, by.y=1)   ## 684
ex1.info <- merge(mex1, info, by.x=1, by.y=1)   ## 603
ex2.info <- merge(mex2, info, by.x=1, by.y=1)   ## 684
lc1.maf <- get.vcf.info(lc1.info$INFO, "AF_EUR")
lc2.maf <- get.vcf.info(lc2.info$INFO, "AF_EUR")
ex1.maf <- get.vcf.info(ex1.info$INFO, "AF_EUR")
ex2.maf <- get.vcf.info(ex2.info$INFO, "AF_EUR")
lc1.maf[lc1.maf>0.5 & !is.na(lc1.maf)] <- 1 - lc1.maf[lc1.maf>0.5 & !is.na(lc1.maf)]
lc2.maf[lc2.maf>0.5 & !is.na(lc2.maf)] <- 1 - lc2.maf[lc2.maf>0.5 & !is.na(lc2.maf)]
ex1.maf[ex1.maf>0.5 & !is.na(ex1.maf)] <- 1 - ex1.maf[ex1.maf>0.5 & !is.na(ex1.maf)]
ex2.maf[ex2.maf>0.5 & !is.na(ex2.maf)] <- 1 - ex2.maf[ex2.maf>0.5 & !is.na(ex2.maf)]

round(summary(lc1.maf),4)
round(summary(lc2.maf),4)
round(summary(ex1.maf),4)
round(summary(ex2.maf),4)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##  0.0026  0.0879  0.1864  0.2148  0.3281  0.4987     2
##  0.0026  0.0919  0.2021  0.2211  0.3386  0.4987     1
##  0.0026  0.0879  0.1864  0.2148  0.3281  0.4987     2
##  0.0026  0.0919  0.2021  0.2211  0.3386  0.4987     1
sum(lc1.maf<0.1, na.rm=T)     ## 174 = 29%
sum(lc1.maf<0.05, na.rm=T)    ##  63 = 10%
sum(lc1.maf<0.01, na.rm=T)    ##   4 =  0.7%
sum(lc2.maf<0.1, na.rm=T)     ## 189 = 28%
sum(lc2.maf<0.05, na.rm=T)    ##  70 = 10%
sum(lc2.maf<0.01, na.rm=T)    ##   4 =  0.6%



##### SNPs on target #####
dat1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Discovery_UK10K_overlap_143-samples_chr20.recode.vcf.gz", stringsAsFactor=F)   ## 4920 x 70
dat2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/PfizerExomes_Replication_UK10K_overlap_143-samples_chr20.recode.vcf.gz", stringsAsFactor=F) ## 6141 x 90
ann1 <- read.delim("/nfs/team151/ym3/exomes/stat/pfizer-pain/functions/Samtools_GATK_11-Nov-10.severe_effects.chr20.txt", skip=14, stringsAsFactors=F)           ## 4922
ann2 <- read.delim("/nfs/team151/ym3/exomes/stat/pfizer-pain/functions/PainReplExomes-Sept11_Rel.severe_effects.chr20.txt", skip=14, stringsAsFactors=F)         ## 6161 
f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", stringsAsFactors=F)                ## 3433 x 63
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", stringsAsFactors=F)              ## 3903 x 82
## remove sample QTL218676
indrm <- names(f1)=="QTL218676"
f1 <- f1[,!indrm]
## edit position
ann1$Location <- gsub("20:", "", ann1$Location) 
ann2$Location <- gsub("20:", "", ann2$Location) 
## overlap between annotations
length(intersect(ann1[,1], ann2[,1]))    ## 2633
## overlap between annotation and exomes
length(intersect(ann1[,2], dat1[,2]))    ## 4920
length(intersect(ann2[,2], dat2[,2]))    ## 6141
## merge annotations with gts concordance
mf1 <- merge(ann1, f1, by.x=2, by.y=1)   ## 3433
mf2 <- merge(ann2, f2, by.x=2, by.y=1)   ## 3903
## gts concordance by site
gts1 <- apply(mf1[,-c(1:14)], 1, mean, na.rm=T)
gts2 <- apply(mf2[,-c(1:14)], 1, mean, na.rm=T)
## consequences
as.matrix(table(mf1$Consequence))
as.matrix(table(mf2$Consequence))
## consequences for discordant sites
as.matrix(table(mf1$Consequence[gts1<0.90]))
## 5PRIME_UTR               1
## INTRONIC                14
## NON_SYNONYMOUS_CODING    3
## SPLICE_SITE              2
## SYNONYMOUS_CODING        6
as.matrix(table(mf2$Consequence[gts2<0.90]))
## INTERGENIC               1
## INTRONIC                20
## NON_SYNONYMOUS_CODING   12
## SPLICE_SITE              3
## SYNONYMOUS_CODING       12

## NAs by consequence
na1 <- apply(mf1[,-c(1:14)], 1, function(v) mean(is.na(v)))
na2 <- apply(mf2[,-c(1:14)], 1, function(v) mean(is.na(v)))
as.matrix(table(mf1$Consequence[na1>0.1]))
## 3PRIME_UTR                2
## 5PRIME_UTR               18
## INTRONIC                 39
## NON_SYNONYMOUS_CODING    63
## SPLICE_SITE               7
## SYNONYMOUS_CODING        63
## UPSTREAM                  3
## WITHIN_NON_CODING_GENE    2
as.matrix(table(mf2$Consequence[na2>0.1]))
## 3PRIME_UTR               17
## 5PRIME_UTR               10
## ESSENTIAL_SPLICE_SITE     2
## INTERGENIC                1
## INTRONIC                158
## NON_SYNONYMOUS_CODING   217
## SPLICE_SITE              28
## STOP_GAINED               2
## SYNONYMOUS_CODING       214
## UPSTREAM                  1
## WITHIN_NON_CODING_GENE    2
as.matrix(table(mf1$Consequence[na1>0.5]))
## 3PRIME_UTR               1
## 5PRIME_UTR               9
## INTRONIC                11
## NON_SYNONYMOUS_CODING   18
## SPLICE_SITE              3
## SYNONYMOUS_CODING       13
## UPSTREAM                 3
as.matrix(table(mf2$Consequence[na2>0.5]))
## 3PRIME_UTR               5
## 5PRIME_UTR               4
## ESSENTIAL_SPLICE_SITE    1
## INTRONIC                58
## NON_SYNONYMOUS_CODING   92
## SPLICE_SITE             14
## SYNONYMOUS_CODING      101


##### Genotype quality of discordant exomes #####
gq1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_gq_discovery.txt", stringsAsFactors=F)
gq2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_gq_replication.txt", stringsAsFactors=F)
f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", stringsAsFactors=F)                ## 3433 x 63
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", stringsAsFactors=F)              ## 3903 x 82
## remove sample QTL218676
indrm <- names(f1)=="QTL218676"
f1 <- f1[,!indrm]
## GQ of exomes for con- and discordant sites
dgq1 <- cgq1 <- NULL
for (j in 1:(ncol(gq1)-1))  {
  dgq1 <- c(dgq1, gq1[,-1][,j][!f1[,-1][,j]])
  cgq1 <- c(cgq1, gq1[,-1][,j][f1[,-1][,j]])
}
dgq2 <- cgq2 <- NULL
for (j in 1:(ncol(gq2)-1))  {
  dgq2 <- c(dgq2, gq2[,-1][,j][!f2[,-1][,j]])
  cgq2 <- c(cgq2, gq2[,-1][,j][f2[,-1][,j]])
}
round(summary(cgq1))
round(summary(dgq1))
round(summary(cgq2))
round(summary(dgq2))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##       0      99      99      93      99      99    5429 
##       0      99      99      90      99      99    5429 
##       0      75      99      85      99      99   26935 
##       0      84      99      86      99      99   26935 
boxplot(cgq1, dgq1, cgq2, dgq2)
mean(cgq1<99, na.rm=T)    ## 18%    ## discovery
mean(dgq1<99, na.rm=T)    ## 21%
mean(cgq1<90, na.rm=T)    ## 16%
mean(dgq1<90, na.rm=T)    ## 18%
mean(cgq2<99, na.rm=T)    ## 33%    ## replication
mean(dgq2<99, na.rm=T)    ## 29%
mean(cgq2<90, na.rm=T)    ## 30%
mean(dgq2<90, na.rm=T)    ## 27%
## concordance by site
sconc1 <- 100*apply(f1[,-1], 1, mean, na.rm=T)
sconc2 <- 100*apply(f2[,-1], 1, mean, na.rm=T)
mgq1 <- apply(gq1[,-1], 1, mean, na.rm=T)
mgq2 <- apply(gq2[,-1], 1, mean, na.rm=T)
#pdf("~/uk10k/plots/gq_concordance_lc_exome_chr20.pdf", width=8, height=4, pointsize=10)
par(mfrow=c(1,2))
boxplot(mgq1[sconc1>=90], mgq1[sconc1<90], ylab="mean GQ", xlab="concordance", names=c(">=90%","<90%"),
        main="Discovery", col=c("darkblue","darkred"), las=1, cex.lab=1.3, cex.main=1.5, cex.axis=1.1)
boxplot(mgq2[sconc2>=90], mgq2[sconc2<90], xlab="concordance", names=c(">=90%","<90%"), col=c("darkblue","darkred"),
        main="Replication", las=1, cex.lab=1.3, cex.main=1.5, cex.axis=1.1)
par(mfrow=c(1,1))
dev.off()


##### Depth of discordant exomes #####
dp1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_dp_discovery.txt", stringsAsFactors=F)
dp2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_dp_replication.txt", stringsAsFactors=F)
f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", stringsAsFactors=F)                ## 3433 x 63
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", stringsAsFactors=F)              ## 3903 x 82
## remove sample QTL218676
indrm <- names(f1)=="QTL218676"
f1 <- f1[,!indrm]
## DP of exomes for con- and discordant sites
ddp1 <- cdp1 <- NULL
for (j in 1:(ncol(dp1)-1))  {
  ddp1 <- c(ddp1, dp1[,-1][,j][!f1[,-1][,j]])
  cdp1 <- c(cdp1, dp1[,-1][,j][f1[,-1][,j]])
}
ddp2 <- cdp2 <- NULL
for (j in 1:(ncol(dp2)-1))  {
  ddp2 <- c(ddp2, dp2[,-1][,j][!f2[,-1][,j]])
  cdp2 <- c(cdp2, dp2[,-1][,j][f2[,-1][,j]])
}
round(summary(c(as.vector(as.matrix(dp1[,-1])), as.vector(as.matrix(dp2[,-1])))))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##       8      29      54      65      88     300   32364
boxplot(c(as.vector(as.matrix(dp1[,-1])), as.vector(as.matrix(dp2[,-1]))))
round(summary(cdp1))
round(summary(ddp1))
round(summary(cdp2))
round(summary(ddp2))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##       8      38      61      69      91     300    5429   ##   discovery  concordant
##       8      22      40      49      65     280    5429                   discordant
##       8      24      48      63      86     300   26935   ## replication  concordant 
##       8      18      38      55      67     295   26935                   discordant
boxplot(cdp1, ddp1, cdp2, ddp2)
mean(cdp1<20, na.rm=T)    ## 8%
mean(ddp1<20, na.rm=T)    ## 21%
mean(cdp2<20, na.rm=T)    ## 19%
mean(ddp2<20, na.rm=T)    ## 28%
## concordance by site
sconc1 <- 100*apply(f1[,-1], 1, mean, na.rm=T)
sconc2 <- 100*apply(f2[,-1], 1, mean, na.rm=T)
mdp1 <- apply(dp1[,-1], 1, mean, na.rm=T)
mdp2 <- apply(dp2[,-1], 1, mean, na.rm=T)
#pdf("~/uk10k/plots/dp_concordance_lc_exome_chr20.pdf", width=8, height=4, pointsize=10)
par(mfrow=c(1,2))
boxplot(mdp1[sconc1>=90], mdp1[sconc1<90], ylab="mean depth", xlab="concordance", names=c(">=90%","<90%"),
        main="Discovery", col=c("darkblue","darkred"), las=1, cex.lab=1.3, cex.main=1.5, cex.axis=1.1, ylim=c(0,275))
boxplot(mdp2[sconc2>=90], mdp2[sconc2<90], xlab="concordance", names=c(">=90%","<90%"), col=c("darkblue","darkred"),
        main="Replication", las=1, cex.lab=1.3, cex.main=1.5, cex.axis=1.1, ylim=c(0,275))
par(mfrow=c(1,1))
dev.off()


##### Do sites with more NAs have lower GQ and lower DP #####
dp1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_dp_discovery.txt", stringsAsFactors=F)              ## 3433 x 62
dp2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_dp_replication.txt", stringsAsFactors=F)            ## 3903 x 82
gq1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_gq_discovery.txt", stringsAsFactors=F)              ## 3433 x 62
gq2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_gq_replication.txt", stringsAsFactors=F)            ## 3903 x 82
f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", stringsAsFactors=F)         ## 3433 x 62
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", stringsAsFactors=F)       ## 3903 x 82
## remove sample QTL218676
indrm <- names(f1)=="QTL218676"
f1 <- f1[,!indrm]
## mean GQ of exomes per site
mgq1 <- apply(gq1[,-1], 1, mean, na.rm=T)
mgq2 <- apply(gq2[,-1], 1, mean, na.rm=T)
## mean DP of exomes per site
mdp1 <- apply(dp1[,-1], 1, mean, na.rm=T)
mdp2 <- apply(dp2[,-1], 1, mean, na.rm=T)
## summaries
round(summary(mgq1))
round(summary(mgq2))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      21      93      99      92      99      99 
##      21      65      97      82      99      99 
round(summary(mdp1))
round(summary(mdp2))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       8      40      64      67      90     212 
##       8      20      44      59      83     274 
boxplot(mgq1, mgq2)
boxplot(mdp1, mdp2)
## number of NAs per site
mna1 <- apply(gq1[,-1], 1, function(v) 100*mean(is.na(v)))
mna2 <- apply(gq2[,-1], 1, function(v) 100*mean(is.na(v)))
## concordance by site
sconc1 <- 100*apply(f1[,-1], 1, mean, na.rm=T)
sconc2 <- 100*apply(f2[,-1], 1, mean, na.rm=T)
## plots
#png("~/uk10k/plots/exome_gq_dp_na_correlations_chr20.png", width=1000, height=600, pointsize=18)
par(mfrow=c(2,3), mai=c(0.9,0.9,0.3,0.1))
plot(mdp1, mna1, xlab="mean depth per site", ylab="percentage of N/A per site", col="blue", cex.lab=1.3, ty="n")
points(mdp1[sconc1>=90], mna1[sconc1>=90], col="blue")
points(mdp1[sconc1<90], mna1[sconc1<90], col="red", pch=15)

plot(mgq1, mna1, xlab="mean GQ per site", ylab="percentage of N/A per site", col="blue", cex.lab=1.3, ty="n")
points(mgq1[sconc1>=90], mna1[sconc1>=90], col="blue")
points(mgq1[sconc1<90], mna1[sconc1<90], col="red", pch=15)

plot(mdp1, mgq1, xlab="mean depth per site", ylab="mean GQ per site", col="blue", cex.lab=1.3, ty="n")
points(mdp1[sconc1>=90], mgq1[sconc1>=90], col="blue")
points(mdp1[sconc1<90], mgq1[sconc1<90], col="red", pch=15)

plot(mdp2, mna2, xlab="mean depth per site", ylab="percentage of N/A per site", col="blue", cex.lab=1.3, ty="n")
points(mdp2[sconc2>=90], mna2[sconc2>=90], col="blue")
points(mdp2[sconc2<90], mna2[sconc2<90], col="red", pch=15)

plot(mgq2, mna2, xlab="mean GQ per site", ylab="percentage of N/A per site", col="blue", cex.lab=1.3, ty="n")
points(mgq2[sconc2>=90], mna2[sconc2>=90], col="blue")
points(mgq2[sconc2<90], mna2[sconc2<90], col="red", pch=15)

plot(mdp2, mgq2, xlab="mean depth per site", ylab="mean GQ per site", col="blue", cex.lab=1.3, ty="n")
points(mdp2[sconc2>=90], mgq2[sconc2>=90], col="blue")
points(mdp2[sconc2<90], mgq2[sconc2<90], col="red", pch=15)
legend(120, 40, legend=c("concordance>=90%","concordance<90%"), col=c("blue","red"), pch=c(1,15))
par(mfrow=c(1,1))
dev.off()


##### Cumulative plot of overall discordance by site discordance #####

## From slide 6 it looks like a small number of sites are very bad.
## So maybe if we dealt with these extremely bad sites we would deal
## with a significant fraction of all the discordances.  I guess it
## would be nice to know what fraction of all the discordances are in
## sites with overall concordance less than 90% or 95%?  Maybe best
## would be a cumulative plot of the fraction of all discordances in
## sites with site discordance below X %.

f1 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_discovery.txt", stringsAsFactors=F)                ## 3433 x 63
f2 <- read.delim("/lustre/scratch101/sanger/kw8/uk10k/pfizer_exome/chr20_comp_lc_exome_replication.txt", stringsAsFactors=F)              ## 3903 x 82
## remove sample QTL218676
indrm <- names(f1)=="QTL218676"
f1 <- f1[,!indrm]
bf1 <- as.matrix(f1[,-1])
bf2 <- as.matrix(f2[,-1])
## overall concordance and discordance
mean(as.vector(as.matrix(bf1)), na.rm=T)     ## 99.5%
mean(!as.vector(as.matrix(bf1)), na.rm=T)    ##  0.5%
## mean concordance per site
m1 <- 100*apply(bf1, 1, mean, na.rm=T)
m2 <- 100*apply(bf2, 1, mean, na.rm=T)
round(summary(m1))
round(summary(m2))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      18     100     100      99     100     100 
##       0     100     100      99     100     100 
## mean discordance per site
d1 <- 100*apply(!bf1, 1, mean, na.rm=T)
d2 <- 100*apply(!bf2, 1, mean, na.rm=T)
round(summary(d1))
round(summary(d2))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##       0       0       0       1       0      82 
##       0       0       0       1       0     100

## overall discordance versus per site concordance
s1 <- s2 <- seq(0,100,1)
disc1 <- sapply(s1, function(v) 100*mean(!as.matrix(bf1[d1>=v,]), na.rm=T))
disc2 <- sapply(s2, function(v) 100*mean(!as.matrix(bf2[d2>=v,]), na.rm=T))
p1 <- sapply(s1, function(v) 100*mean(d1>=v))
p2 <- sapply(s2, function(v) 100*mean(d2>=v))
n1 <- sapply(s1, function(v) sum(d1>=v))
n2 <- sapply(s2, function(v) sum(d2>=v))

#pdf("~/uk10k/plots/discordance_lc_exome_discovery_chr20.pdf", width=4.5, height=4.5, pointsize=11)
plot(s1, disc1, ty="b", col="red", xlab="discordance by site", ylab="percentage", pch=17, cex.lab=1.3, las=1, ylim=c(0,100))
points(s1, p1, ty="b", col="blue", pch=18)
legend(45,40, legend=c("overall discordance","sites"), col=c("red","blue"), pch=c(17,18))
dev.off()
#pdf("~/uk10k/plots/discordance_lc_exome_replication_chr20.pdf", width=4.5, height=4.5, pointsize=11)
plot(s2, disc2, ty="b", col="red", xlab="discordance by site", ylab="percentage", pch=17, cex.lab=1.3, las=1, ylim=c(0,100))
points(s2, p2, ty="b", col="blue", pch=18)
legend(45,40, legend=c("overall discordance","sites"), col=c("red","blue"), pch=c(17,18))
dev.off()

cbind(s1, round(disc1,1), n1, round(p1,1))
cbind(s2, round(disc2,1), n2, round(p2,1))
