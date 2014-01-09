source("~/R_dir/R_scripts/load_nsta_functions.R")


## bsub -o qqplot_batch_effect_sanger-bgi.out -q hugemem -R "select[mem>=28000] rusage[mem=28000]" -M28000000 R CMD BATCH ~/R_dir/R_scripts_uk10k_2013/qqplot_batch_effect_sanger-bgi.R


## Take qqplot from Yasin's code:                 
## /software/rarevar/plotting/QQPlot.R
QQplot <- function (pval, xl, yl) {
  ## prepare data
  ord <- order(pval, na.last=NA)
  P <- pval[ord]
  ## allocate min to P==0
  P[P==0] <- min(P[P>0])
  obspval <- P
  logobspval <- -(log10(obspval))
  exppval <- c(1:length(obspval))
  logexppval <- -(log10( (exppval-0.5)/length(exppval)))
  obsmax <- trunc(max(logobspval))+1
  expmax <- trunc(max(logexppval))+1
  ## plot
  plot(c(0,expmax), c(0,expmax), col="gray", lwd=1, type="l", xlab="Expected -log10 P-value", ylab="Observed -log10 P-value",
       #xlim=c(0,expmax), ylim=c(0,obsmax), las=1, xaxs="i", yaxs="i", bty="l", cex.lab=1.2)
       xlim=c(0,xl), ylim=c(0,yl), las=1, xaxs="i", yaxs="i", bty="l", cex.lab=1.2)
  points(logexppval, logobspval, pch=23, cex=0.4, bg="black")
}


##### Plink case-control Sanger-BGI #####

d <- dir("/lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/output", full.names=T)
ind <- grep("assoc", d)
d <- d[ind]

### all p-values
## p <- NULL
## for (i in 1:length(d))  {
##   print(i)
##   f <- read.table(d[i], stringsAsFactors=F, header=T)
##   ## QQplots for p-values
##   fn <- paste("/lustre/scratch113/projects/uk10k/users/kw8/plots2/qqplots_batch/qqplot_batch_sanger-bgi_REL-2012-06-02_chr",i,".png", sep="")
##   png(fn, width=600, height=600, pointsize=16)
##   QQplot(f$P)
##   dev.off()
## }

### all p-values, zoomed into p<1e-5
## p <- NULL
## for (i in 1:length(d))  {
##   print(i)
##   f <- read.table(d[i], stringsAsFactors=F, header=T)
##   ## QQplots for p-values
##   fn <- paste("/lustre/scratch113/projects/uk10k/users/kw8/plots2/qqplots_batch/qqplot_batch_sanger-bgi_p1-9_chr",i,".png", sep="")
##   png(fn, width=600, height=600, pointsize=16)
##   QQplot(f$P, xl=9, yl=9)
##   dev.off()
## }

p <- NULL
for (i in 1:length(d))  {
  print(i)
  f <- read.table(d[i], stringsAsFactors=F, header=T)
  ## QQplots for p-values
  fn <- paste("/lustre/scratch113/projects/uk10k/users/kw8/plots2/qqplots_batch/qqplot_batch_sanger-bgi_p1-5_chr",i,".png", sep="")
  png(fn, width=600, height=600, pointsize=16)
  QQplot(f$P, xl=5, yl=5)
  dev.off()
}

