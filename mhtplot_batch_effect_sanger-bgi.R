#source("/nfs/users/nfs_k/kw8/R_dir/R_scripts/load_nsta_functions.R")


## bsub -o mhtplot_batch_effect_sanger-bgi.out -q hugemem -R "select[mem>=28000] rusage[mem=28000]" -M28000000 R CMD BATCH ~/R_dir/R_scripts_uk10k_2013/mhtplot_batch_effect_sanger-bgi.R


mht.plot <- function(chr, pos, p, chrlengths)  {
  ## remove negative p-values
  indrm <- p<0
  chr <- chr[!indrm]
  pos <- pos[!indrm]
  p <- p[!indrm]
  ## allocate max p-value to p==0
  ind <- p==0
  p[ind] <- max(p[p>0])
  ## get re-scaled x-coordinates and colours
  xchr <- c(0, cumsum(chrlengths/1000000))
  xpos <- pos/1000000
  mycols <- rep("darkblue", length(pos))
  for (i in 1:22)  {
    xpos[chr==i] <- xchr[i] + pos[chr==i]/1000000
    if (is.element(i, seq(2,22,2)))  {
      mycols[chr==i] <- "lightblue"
    }
  }
  ## get plot limits
  xmax <- max(xpos)
  ymax <- round(max(-log10(p[p>0]), na.rm=T))
  ## plot -log10 p-values
  plot(xpos, -log10(p), axes=F, xlab="chromosome", cex.lab=1.5, cex.main=1.9,
       xlim=c(0,xmax), ylim=c(0,ymax), col=mycols, pch=18, cex=1.2,)
  box()
  axis(1, at=xchr[-23], labels=c(1:22))
  yticks <- round(seq(min(-log10(p[p>0]), na.rm=T), max(-log10(p[p>0]), na.rm=T), length.out=5))
  axis(2, at=yticks, labels=yticks, las=1)
}


##### Manhattan plots for Plink case-control Sanger-BGI, separate by chromosome #####

chrsize <- read.delim("/nfs/users/nfs_k/kw8/data/chr_size/human_chr_length_ncbi37_hg19.txt", stringsAsFactors=F, header=F)
#d <- dir("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/output", full.names=T)
d <- dir("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/gzip", full.names=T)
ind <- grep("assoc", d)
d <- d[ind]
## loop through chroms
for (i in 2:length(d))  {
  print(i)
  f <- read.table(d[i], stringsAsFactors=F, header=T)
  ## Manhattan plots for p-values
  #fn <- paste("/lustre/scratch113/projects/uk10k/users/kw8/plots2/mhtplots_batch/mhtplot_batch_sanger-bgi_REL-2012-06-02_chr",i,".png", sep="")
  fn <- paste("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PLOTS/mhtplot_batch_sanger-bgi_REL-2012-06-02_chr",i,".png", sep="")
  png(fn, width=1200, height=600, pointsize=16)
  m <- paste("Chr",i)
  xl <- chrsize[i,2]/1000000
  plot(f$BP/1000000, -log10(f$P), col="blue", main=m, xlab="position (in Mb)", ylab="-log10(p-value)", xlim=c(0,xl))
  dev.off()
rm(f)
gc()
}


##### Manhattan plots for Plink case-control Sanger-BGI only for p< 1e-5 #####

#d <- dir("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/output", full.names=T)
d <- dir("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/gzip", full.names=T)
ind <- grep("assoc", d)
d <- d[ind]
## loop through chroms
fdat <- NULL
for (i in 1:length(d))  {
  print(i)
  f <- read.table(d[i], stringsAsFactors=F, header=T)
  ## collate p-values
  ind <- which(-log10(f$P)>5)
  dat <- cbind(f$CHR[ind], f$BP[ind], f$P[ind])
  fdat <- rbind(fdat, dat)          ## 115019
rm(dat)
gc() 
}

## plot
chrsize <- read.delim("/nfs/users/nfs_k/kw8/data/chr_size/human_chr_length_ncbi37_hg19.txt", stringsAsFactors=F, header=F)
#png("/lustre/scratch113/projects/uk10k/users/kw8/plots2/mhtplots_batch/mhtplot_batch_sanger-bgi_REL-2012-06-02.png", width=1500, height=600, pointsize=16)
mht.plot(chr=fdat[,1], pos=fdat[,2], p=fdat[,3], chrlengths=chrsize[1:22,2]) 
dev.off()


## loop through chroms
for (i in 1:22)  {
  print(i)
  ind <- which(fdat[,1]==i)
 # fn <- paste("/lustre/scratch113/projects/uk10k/users/kw8/plots2/mhtplots_batch/mhtplot_batch_sanger-bgi_REL-2012-06-02_chr",i,".png", sep="")
  fn <- paste("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PLOTS/mhtplot_batch_sanger-bgi_REL-2012-06-02_chr",i,".png", sep="")
  png(fn, width=1200, height=600, pointsize=16)
  m <- paste("Chr",i)
  xl <- chrsize[i,2]/1000000
  plot(fdat[ind,2]/1000000, -log10(fdat[ind,3]), col="blue", main=m, xlab="position (in Mb)", ylab="-log10(p-value)", xlim=c(0,xl))
  dev.off()
}

