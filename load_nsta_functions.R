# set options so that 10000 is not diplayed as 1e+05 !!!
options(scipen=3)

## sample ids to be excluded for BEAGLE analysis
diffname <- c("NA10851","NA12004","NA12414","NA12717")


binary <- function(x, left=TRUE, dp=10) {
  b <- NULL
  for (i in 2^(0:15)) {
    x.new <- floor(x/2)
    b <- c(b,x - 2*x.new)
    x <- x.new
  }
  if (left)  {
    r <- b[1:dp]   ## for BAM
  }  else  {
    r <- b[dp:1]
  }
  r
}
#flag40 <- binary(4*16)
#flag800 <- binary(8*16^2)
#if (any(binary(x) & flag40)) {
#  # flag hex 40 is set
#}

hexadecimal <- function(x) {
  b <- NULL
  for (i in 16^(0:4)) {
    x.new <- floor(x/16)
    b <- c(b,x - 16*x.new)
    x <- x.new
  }
  #b[4:1]
  b[1:3]   ## for BAM
}

## bin <- t(sapply(1:400, binary))
## bin <- as.vector((apply(bin, 1, function(v) paste(v, sep="", collapse=""))))
## hex <- t(sapply(1:400, hexadecimal))
## hex <- t(apply(hex, 1, function(v) gsub(10, "A", v)))
## hex <- t(apply(hex, 1, function(v) gsub(11, "B", v)))
## hex <- t(apply(hex, 1, function(v) gsub(12, "C", v)))
## hex <- t(apply(hex, 1, function(v) gsub(13, "D", v)))
## hex <- t(apply(hex, 1, function(v) gsub(14, "E", v)))
## hex <- t(apply(hex, 1, function(v) gsub(15, "F", v)))
## hex <- as.vector((apply(hex, 1, function(v) paste(v, sep="", collapse=""))))
## dat <- data.frame(1:400, hex, bin)
## names(dat) <- c("decimal","hexadecimal","binary")
## ## write.table(dat, "~/project_1000_genomes/analysis/dec2hex2bin.txt", quote=F, row.names=F, sep="\t")
## write.table(dat, "~/project_1000_genomes/analysis/dec2hex2bin_bam.txt", quote=F, row.names=F, sep="\t")

#### Edit ASCII table
## http://en.wikipedia.org/wiki/ASCII#ASCII_printable_characters

## infile <- "~/data/ascii/ascii_table.txt"
## com <- paste("cut -f1-5", infile, "> ~/data/ascii/bin_oct_dec_hex.txt")
## system(com)
## com <- paste("cut -f6", infile, "> ~/data/ascii/ascii.txt")
## system(com)

bin.oct.dec.hex <- read.delim("~/data/ascii/bin_oct_dec_hex.txt")
ascii <- scan("~/data/ascii/ascii.txt", what="character", na.strings="", sep="\n")
ascii <- c(" ", ascii) 
assign("bin.oct.dec.hex", bin.oct.dec.hex, env = .GlobalEnv)
assign("ascii", ascii, env = .GlobalEnv)

ascii2dec <- function(s)  {
  v <- unlist(strsplit(s, ""))
  ind <- sapply(v, function(v) which(v==as.character(ascii)))
  dec <- as.numeric(as.character(bin.oct.dec.hex[ind,3]))
  dec
}

##### convert string to a vector of real numbers

bp2real.n <- function(str)
{
   charvec <- unlist(strsplit(str, ""))
   vec <- gsub("A", "0", charvec)
   vec <- gsub("C", "1", vec)
   vec <- gsub("G", "2", vec)
   vec <- gsub("T", "3", vec)
   vec <- gsub("N", "4", vec)

   vec <- gsub("a", "0", vec)
   vec <- gsub("c", "1", vec)
   vec <- gsub("g", "2", vec)
   vec <- gsub("t", "3", vec)
   vec <- gsub("n", "4", vec)

   vec <- as.real(vec)
   vec
}

##### gets reverse complement
# returns integer vector
# A=0, C=1, G=2, T=3, N=4

get.reverse.compl <- function(str)  
{
  v <- bp2real.n(str)
  l <- length(v)
  ord <- seq(l, 1, -1)
  rv <- v[ord]
  
  ind0 <- rv==0
  ind1 <- rv==1  
  ind2 <- rv==2  
  ind3 <- rv==3  
  ind4 <- rv==4  
  
  rv[ind0] <- 3
  rv[ind1] <- 2
  rv[ind2] <- 1
  rv[ind3] <- 0
  rv[ind4] <- 4
  rv
}

##### changes real into bp

real2bp <- function(vec)
{
   vec <- gsub("0", "A", vec)
   vec <- gsub("1", "C", vec)
   vec <- gsub("2", "G", vec)
   vec <- gsub("3", "T", vec)
   vec <- gsub("4", "N", vec)

   fvec <- paste(vec, sep="", collapse="")
   fvec
}

##### changes a vector of sequence strings into reverse complements
##### gets reverse complement as bp

seq2revcompl <- function(sequ)
{
  rc <- get.reverse.compl(sequ)
  rcs <- real2bp(rc)
  rcs
}

### vecstr2matstr combined with str2vec ########################
### converts a vector of strings into a matrix of strings by ###
### splitting the string by a given pattern ####################

# string = string that needs splitting
# pat = pattern by which to split

str2vec <- function(string, pat) 
{
  charvec <- unlist(strsplit(string, pat))
  charvec
}

# vecstring = vector of string 
# pat = pattern by which to split
# returns matrix/list that contains split up strings

vecstr2matstr <- function(vecstring, pat)
{
  matstring <- t(apply(as.matrix(vecstring), 1, str2vec, pat))
  matstring
}

# vecstring = vector of string 
# pat = pattern by which to split
# returns matrix that contains split up strings

vecstr2liststr <- function(vecstring, pat)
{
  liststring <- lapply(as.matrix(vecstring), str2vec, pat)
  liststring
}

#### extract header/identifier
# fn = file name
# table = logical to indicate whether file is read by lines or as a table

get.header <- function(fn,table)  {
  if(table==F)  {  
    dna <- as.character(fn)
    h.ind <- grep(">",dna)
    dna.h <- dna[h.ind]
  }  else  {
    dna.h <- as.character(fn[,1])
  }
  dna.h
}

#### extract sequence
# fn = file name
# table = logical to indicate whether file is read by lines or as a table

get.sequence <- function(fn,table)  {
  if(table==F)  {  
    dna <- toupper(as.character(fn))
    h.ind <- grep(">",dna)
    s.ind <- h.ind+1
    dna.s <- dna[s.ind]
  }  else  {
    dna.s <- toupper(as.character(fn[,2]))
  }
  dna.s
}

### Change from fasta format to table format
# fasta = data in fasta format
fasta2tab2 <- function(fasta)  {
  indh <- grep(">",fasta)  
  header <- as.character(fasta[indh])
  indmat <- cbind(indh, c(indh[2:length(indh)], (length(fasta)+1)))
  inds <- apply(indmat, 1, function(v) (v[1]+1):(v[2]-1))
  sequence <- sapply(inds, function(ind) paste(fasta[ind], collapse=""))
  tab.df <- cbind(header,sequence)
  tab.df
}

#### Remove duplicate entries of column cn, in dataframe, needs to be sorted
# x.df = data frame
# cn = column number that contains duplicates
get.unique.df <- function(x.df, cn)  {
  ## sort data rame by cn
  ord <- order(x.df[,cn])
  x.df <- x.df[ord,]
  x <- as.character(x.df[,cn])
  y <- c(-1,x[1:(length(x)-1)])
  ind <- x == y
  r.df <- x.df[!ind,]
  r.df
}
## return row indices
get.unique.df.id <- function(x.df, cn)  { 
  x <- as.character(x.df[,cn])
  y <- c(-1,x[1:(length(x)-1)])
  ind <- x == y
  !ind
}

#### Remove duplicate entries of columns (!) cn, in dataframe, needs to be sorted
get.unique.df2 <- function(x.df, cn)  { 
  x <- apply(x.df[,cn], 1, paste, collapse="_")
  x <- gsub(" ", "", x)
  y <- c(-1,x[1:(length(x)-1)])
  ind <- x == y
  r.df <- x.df[!ind,]
  r.df
}

### transforms list with equal number of elements into matrix
# wlist = list
# cnum = number of columns

list2mat <- function(wlist, cnum)  {
  rnum <- length(wlist)
  v <- unlist(wlist)
  mat <- matrix(v, ncol=cnum, nrow=rnum, byrow=T)
  mat
}

### Create ID by pasting chr, start, end together

create.id <- function(dat)  {
  chr.pos <- t(apply(dat[,1:3], 1, as.character))
  pos <- t(apply(chr.pos[,2:3], 1, as.numeric))
  id <- apply(pos, 1, paste, sep="", collapse="-")
  id <- apply(data.frame(chr.pos[,1], id), 1, paste, sep="", collapse="-")
  id
}


### Use empirical distribution of inserts e.g. from chr20 to estimate
### a cutoff (e.g. 0.001) for anomalous RPs 
# insert = vector of inserts
# thresh = threshold for cutoff as proportion of all reads, e.g. 0.001=0.1%

get.insert.cutoff <- function(insert, thresh)  {
  n <- length(insert)
  tab <- table(insert)
  ftab <- cbind(as.numeric(names(tab)), tab)
  ## cumulative sum starting from end  (= largest distance)
  s <- cumsum(as.numeric(ftab[nrow(ftab):1,2]))
  p <- s/n
  ## find insert length where proportion of anomalous RPs > thresh
  ind <- length(which(p > thresh))
  ## in case there is not an entry for each insert size 
  ftab[ind,1]+1
}

### Plot in chunks
# x1, x2 = interval on chrX to be plotted
# p1, p2 = start and end position of reads (leftmost positions)
# dist = distance (=offset) between read pairs (leftmost positions)
# pscore = single or paired read alignment score
# maxlim = upper limit for distance that is plotted

plot.interval <- function(x1,x2, p1, p2, dist, pscore, maxlim=NULL)  {
  intvl <- c(x1:x2)
  ## select positions, offset and scores so that position of first or second read is within interval
  ind <- (p1 >= x1 & p1 <= x2) | (p2 >= x1 & p2 <= x2)
  sp1 <- p1[ind]
  sp2 <- p2[ind] + 35     # add read length
  sdist <- dist[ind]
  spscore <- pscore[ind]
  ## index reads with offset > 300
  if (is.null(maxlim))  {
    indoff <- abs(sdist) > 300
  }  else  {
    indoff <- abs(sdist) > 300 & abs(sdist) < maxlim
  }
  nbp1 <- sp1[!indoff]
  nbp2 <- sp2[!indoff]
  abp1 <- sp1[indoff]
  abp2 <- sp2[indoff]
  ## create random number for y-coord of abnormal read pairs
  n <- length(sp1)
  m <- length(abp1)
  r <- trunc(runif(n, 1, 200))
  ## plot
  sp1 <- sp1/1000000
  nbp1 <- nbp1/1000000
  nbp2 <- nbp2/1000000
  abp1 <- abp1/1000000
  abp2 <- abp2/1000000
  x1 <- x1/1000000
  x2 <- x2/1000000
  plot(sp1, spscore, ylim=c(min(pscore[ind]),max(pscore[ind])), xlim=c(x1,x2),
       xlab="position", ylab="read score", ty="n", cex.lab=1.3)
  for (i in 1:n) {
    lines(c(nbp1[i],nbp2[i]),c(spscore[i],spscore[i]), col="black", lwd=1)
  }
  for (i in 1:m) {
    #lines(c(abp1[i],abp2[i]),c(r[i],r[i]), col="red", lwd=1)
    lines(c(abp1[i],abp2[i]),c(spscore[i],spscore[i]), col="red", lwd=1.2)
  }
}


## Physical coverage
# frag.mean = mean insert size
# seq.cov = sequence coverage
# total.reads = number of total reads
# bp = length of reads
# hg.len = human genome size
get.pcov <- function(insertmean, seq.cov=NULL, total.reads=NULL, bp=36, hg.len=3000000000)  {
  ## Total number of reads depends on sequence coverage and vice-versa
  if (is.null(seq.cov))  {
    seq.cov <-  (total.reads * bp) / hg.len          ## for haploid
    #print(seq.cov)
  }  else if (is.null(total.reads))  {
    total.reads <- seq.cov * hg.len / bp
    #print(total.reads)
  }  else  {
    print("Error: seq.cov or total.reads must be specified")
  }
  ## Physical coverage of region where breakpoint can be discovered,
  ## depends on total number of reads and fragment length
  (pcov <- (total.reads * (insertmean - 2*bp)) / (hg.len * 2))    ## haploid, fragments=reads/2
  list(pcov=pcov, scov=seq.cov, reads=total.reads)
}


### Calculate depth within an interval for a single interval
# x = vector with start and end
# p = all positions available
# insertsize = length of actually sequenced DNA

get.depth <- function(x, p, insertsize)  {
  p1 <- min(x)
  p2 <- max(x)
  d <- p2 - p1 + 1
  ind <- (p >= p1) & (p <= p2)
  np <- length(p[ind]) * insertsize
  dep <- np/d
  dep
}


### Calculate depth within an interval for multiple intervals
# x = matrix with starts and ends
# p = all positions available
# insertsize = length of actually sequenced DNA

get.all.depth <- function(x, p, insertsize)  {
  dep <- apply(x, 1, get.depth, p, insertsize)
  dep
}


### Same as get.depth, just normalised

get.density <- function(x, p)  {
  p1 <- min(x)
  p2 <- max(x)
  d <- p2 - p1 + 1
  ind <- (p >= p1) & (p <= p2)
  np <- length(p[ind])
  dep <- np/d
  dep
}


### Same as get.all.depth, just normalised

get.all.density <- function(x, p)  {
  dep <- apply(x, 1, get.density, p)
  dep
}

### Same as get.depth, just counts for fixed interval

get.counts <- function(x, p)  {
  p1 <- min(x)
  p2 <- max(x)
  d <- p2 - p1 + 1
  ind <- (p >= p1) & (p <= p2)
  np <- length(p[ind])
  np
}


### Same as get.all.depth, just counts for fixed interval

get.all.counts <- function(x, p)  {
  dep <- apply(x, 1, get.counts, p)
  dep
}



get.abnormal.offset2 <- function(p1, p2, thresh)  {
  ind <- abs(p2-p1) > thresh
  abp1 <- p1[ind]
  abp2 <- p2[ind]

  dep1 <- NULL
  for (i in 1:length(abp1))  {
    dep1[i] <- get.density(abp1[i], abp2[i], p1)
  }
  dep2 <- NULL
  for (i in 1:length(abp1))  {
    dep2[i] <- get.density(abp1[i], abp2[i], p2)
  }
  dep <- (dep1+dep2)/2
  dep
}

get.abnormal.offset <- function(p1, p2, thresh)  {
  ind <- abs(p2-p1) > thresh
  abmat <- cbind(p1[ind], p2[ind])
  n <- nrow(abmat)
  
  dep <- apply(abmat, 1, get.density, p1)
  dep
}

get.depth.abnormal.offset <- function(p1, p2, thresh, insertsize)  {
  ind <- abs(p2-p1) > thresh
  abmat <- cbind(p1[ind], p2[ind])
  n <- nrow(abmat)
  
  dep <- apply(abmat, 1, get.depth, p1, insertsize)
  dep
}


get.overlap <- function(p1, p2, thresh)  {
  pover <- NULL
  for (i in 1:length(p1))  {
    ind1 <- which(abs(p1-p1[i]) <= thresh)
    ind2 <- which(abs(p2-p2[i]) <= thresh)
    ind <- intersect(ind1,ind2)
    pover[i] <- list(ind)
  }
  pover
}

## remove positions with identical start and end positions
rm.identical.pos <- function(p1, p2)  {
  pover <- NULL
  for (i in 1:length(p1))  {
    ind1 <- p1-p1[i]==0
    ind2 <- p2-p2[i]==0
    ind <- which(ind1 & ind2)
    pover[i] <- list(ind)
  }
  pover
}

## thresh.low=0 necessary for merge.clusters.trans !!!!!
get.overlap.0 <- function(p1, p2, thresh.low=0, thresh.up)  {
  pover <- NULL
  for (i in 1:length(p1))  {
    ind1 <- abs(p1-p1[i])<=thresh.up & abs(p1-p1[i])>=thresh.low
    ind2 <- abs(p2-p2[i])<=thresh.up & abs(p2-p2[i])>=thresh.low
    ind <- which(ind1 & ind2)
    pover[i] <- list(ind)
  }
  pover
}

get.overlap.0.mat <- function(p1, p2, thresh.low=0, thresh.up)  {
  n <- length(p1)
  over <- matrix(0, nrow=n, ncol=n)
  diag(over) <- 1
  for (i in 1:length(p1))  {
    ind1 <- abs(p1-p1[i])<=thresh.up & abs(p1-p1[i])>=thresh.low
    ind2 <- abs(p2-p2[i])<=thresh.up & abs(p2-p2[i])>=thresh.low
    ind <- which(ind1 & ind2)
    if (length(ind)!=0)  {
      over[i,ind] <- 1
    }
  }
  over
}

### within one matrix
### store result in matrix instead of list as in get.overlaps4
get.overlap.starts <- function(p1,p2,thresh)  {
  n <- length(p1)
  over <- matrix(0, nrow=n, ncol=n)
  diag(over) <- 1
  for (i in 1:n)  {
    ind1 <- abs(p1-p1[i])<=thresh
    ind2 <- abs(p2-p2[i])<=thresh
    ind <- which(ind1 & ind2)
    if (length(ind)!=0)  {
      over[i,ind] <- 1
    }
  }
  over
}

### mat1 is the reference
get.overlaps2 <- function(mat1, mat2)  {
  bg1 <- mat1[,1]
  ed1 <- mat1[,2]
  bg2 <- mat2[,1]
  ed2 <- mat2[,2]  
  pover <- NULL
  for (i in 1:length(bg1))  {
    ind1 <- which((bg2 >= bg1[i]) & (bg2 <= ed1[i]))
    ind2 <- which((ed2 >= bg1[i]) & (ed2 <= ed1[i]))
    ind3 <- which((bg2 <= bg1[i]) & (ed2 >= bg1[i]))
    #ind4 <- which((bg2 <= ed1[i]) & (ed2 >= ed1[i]))
    ind12 <- union(ind1, ind2)
    ind <- union(ind12, ind3)
    if (length(ind)==0)
      ind <- NA
    pover[i] <- list(ind)
  }
  pover
}

### get e.g. thresh=75% overlap
get.overlaps.thresh2 <- function(mat1, mat2, thresh)  {
  bg1 <- mat1[,1]
  ed1 <- mat1[,2]
  bg2 <- mat2[,1]
  ed2 <- mat2[,2]
  pover <- NULL
  for (i in 1:length(bg1))  {
    ## [bg1,ed1] overlaps [bg2,ed2] but [bg1,ed1] is shifted to the left
    ind1 <- which(bg2 >= bg1[i] & bg2 <= ed1[i] & ed2 >= ed1[i] & (abs(ed1[i]-bg2) / abs(ed1[i]-bg1[i]) >= thresh))
    ## [bg1,ed1] contains [bg2,ed2]
    ind2 <- which(bg2 >= bg1[i] & bg2 <= ed1[i] & ed2 < ed1[i] & (abs(ed2-bg2) / abs(ed1[i]-bg1[i]) >= thresh))
    ## [bg2,ed2] contains [bg1,ed1]
    ind3 <- which(bg2 < bg1[i] & ed2 >= ed1[i])
    ## [bg1,ed1] overlaps [bg2,ed2] but [bg1,ed1] is shifted to the right
    ind4 <- which(bg2 < bg1[i] & ed2 >= bg1[i]  & ed2 < ed1[i] & (abs(ed2-bg1[i]) / abs(ed1[i]-bg1[i]) >= thresh))
    ## combine
    ind12 <- union(ind1, ind2)
    ind13 <- union(ind12, ind3)
    ind <- union(ind13, ind4)
    if (length(ind)==0)
      ind <- NA
    pover[i] <- list(ind)
  }
  pover
}

#### Gets reciprocal overlap using "thresh" parameter (default=0.5), returns a vector of indices
## v = vector of start and end
## mat = matrix of starts and ends
get.reciprocal.overlaps <- function(v, mat, thresh=0.5)  {
  bg1 <- as.numeric(v[1])
  ed1 <- as.numeric(v[2])
  bg2 <- as.numeric(mat[,1])
  ed2 <- as.numeric(mat[,2])
  ind0 <- which(bg1==bg2 & ed1==ed2)
  ind1 <- which((bg1 <= bg2) & (ed1 <= ed2) & ((ed1-bg2)/(ed1-bg1) > thresh) & ((ed1-bg2)/(ed2-bg2) > thresh))
  ## [bg1,ed1] shifted left
  ind2 <- which((bg1 >= bg2) & (ed1 >= ed2) & ((ed2-bg1)/(ed1-bg1) > thresh) & ((ed2-bg1)/(ed2-bg2) > thresh))
  ## [bg1,ed1] shifted right
  ind3 <- which((bg1 <= bg2) & (ed1 >= ed2) & (ed2-bg2)/(ed1-bg1) > thresh)
  ## [bg1,ed1] includes
  ind4 <- which((bg1 >= bg2) & (ed1 <= ed2) & (ed1-bg1)/(ed2-bg2) > thresh)
  ## [bg1,ed1] is included
  ind01 <- union(ind0, ind1)
  ind02 <- union(ind01, ind2)
  ind03 <- union(ind02, ind3)
  ind <- union(ind03, ind4)
  ind
}


### Returns reciprocal overlaps as percentages
get.reciprocal.overlaps.p <- function(v, mat, thresh=0.5)  {
  overlap <- NULL
  bg1 <- as.numeric(v[1])
  ed1 <- as.numeric(v[2])
  bg2 <- as.numeric(mat[,1])
  ed2 <- as.numeric(mat[,2])
  ## bgs and eds are identical
  ind0 <- which(bg1==bg2 & ed1==ed2)
  overlap[ind0] <- 100
  ## [bg1,ed1] shifted left
  ind1 <- which((bg1 <= bg2) & (ed1 <= ed2) & ((ed1-bg2)/(ed1-bg1) > thresh) & ((ed1-bg2)/(ed2-bg2) > thresh))
  overlap[ind1] <-  100 * min(c((ed1[ind1]-bg2[ind1])/(ed1[ind1]-bg1[ind1]), ((ed1[ind1]-bg2[ind1])/(ed2[ind1]-bg2[ind1]))))
  ## [bg1,ed1] shifted right
  ind2 <- which((bg1 >= bg2) & (ed1 >= ed2) & ((ed2-bg1)/(ed1-bg1) > thresh) & ((ed2-bg1)/(ed2-bg2) > thresh))
  overlap[ind2] <-  100 * min(c((ed2[ind2]-bg1[ind2])/(ed1[ind2]-bg1[ind2]), (ed2[ind2]-bg1[ind2])/(ed2[ind2]-bg2[ind2])))
  ## [bg1,ed1] includes
  ind3 <- which((bg1 <= bg2) & (ed1 >= ed2) & (ed2-bg2)/(ed1-bg1) > thresh)
  overlap[ind3] <- 100 * (ed2[ind3]-bg2[ind3])/(ed1[ind3]-bg1[ind3])
  ## [bg1,ed1] is included
  ind4 <- which((bg1 >= bg2) & (ed1 <= ed2) & (ed1-bg1)/(ed2-bg2) > thresh)
  overlap[ind4] <- 100 * (ed1[ind4]-bg1[ind4])/(ed2[ind4]-bg2[ind4])
  ind01 <- union(ind0, ind1)
  ind02 <- union(ind01, ind2)
  ind03 <- union(ind02, ind3)
  ind <- union(ind03, ind4)
  list(ind=ind, overlap=overlap)
}

get.reciprocal.overlaps.old <- function(v, mat, thresh=0.5, d=500)  {
  bg1 <- as.numeric(v[1])
  ed1 <- as.numeric(v[2])
  bg2 <- as.numeric(mat[,1])
  ed2 <- as.numeric(mat[,2])
  ind0 <- bg1==bg2 & ed1==ed2
  ## [bg1,ed1] shifted left
  ind1 <- (bg1 <= bg2) & (ed1 <= ed2) & ((ed1-bg2)/(ed1-bg1) > thresh) & ((ed1-bg2)/(ed2-bg2) > thresh)
  ## [bg1,ed1] shifted right
  ind2 <- (bg1 >= bg2) & (ed1 >= ed2) & ((ed2-bg1)/(ed1-bg1) > thresh) & ((ed2-bg1)/(ed2-bg2) > thresh)
  ## [bg1,ed1] includes
  ind3 <- (bg1 <= bg2) & (ed1 >= ed2) & (ed2-bg2)/(ed1-bg1) > thresh
  ## [bg1,ed1] is included
  ind4 <- (bg1 >= bg2) & (ed1 <= ed2) & (ed1-bg1)/(ed2-bg2) > thresh
  ## distance of breakpoints less than d
  ind5 <- (abs(bg1 - bg2) <= d) & (abs(ed1 - ed2) <= d)
  ind <- which((ind0 | ind1 | ind2 | ind3 | ind4) & ind5)
  ind
}

# mat1 <- rbind(c(11,20),c(11,21),c(12,30))
# mat2 <- rbind(c(15,20),c(16,21),c(22,30))
# apply(mat1, 1, get.reciprocal.overlaps, mat=mat2, thresh=0.5)


### Get reciprocal overlap of thresh % (default=0.5)
### and compares origin for 1kG for data release
get.reciprocal.overlaps.origin <- function(v, mat, thresh=0.5, origin)  {
  bg1 <- as.numeric(as.character(v[1]))
  ed1 <- as.numeric(as.character(v[2]))
  orig <- as.character(v[3])
  bg2 <- mat[,1]
  ed2 <- mat[,2]
  ind0 <- which(bg1==bg2 & ed1==ed2)
  ## [bg1,ed1] shifted left
  ind1 <- which((bg1 <= bg2) & (ed1 <= ed2) & ((ed1-bg2)/(ed1-bg1) > thresh) & ((ed1-bg2)/(ed2-bg2) > thresh))
  ## [bg1,ed1] shifted right
  ind2 <- which((bg1 >= bg2) & (ed1 >= ed2) & ((ed2-bg1)/(ed1-bg1) > thresh) & ((ed2-bg1)/(ed2-bg2) > thresh))
  ## [bg1,ed1] includes
  ind3 <- which((bg1 <= bg2) & (ed1 >= ed2) & (ed2-bg2)/(ed1-bg1) > thresh)
  ## [bg1,ed1] is included
  ind4 <- which((bg1 >= bg2) & (ed1 <= ed2) & (ed1-bg1)/(ed2-bg2) > thresh)
  ind01 <- union(ind0, ind1)
  ind02 <- union(ind01, ind2)
  ind03 <- union(ind02, ind3)
  ind <- union(ind03, ind4)
  indorigin <- which(orig==origin[ind])
  find <- ind[indorigin]
  find
}

#### Takes index list from get.reciprocal.overlaps and converts into a matrix
#### which indicates which entries overlap
## ind = list or matrix of overlap indices

fill.overlap.matrix <- function(ind)  {     # changed to incorporate ind as matrix 2010-01-05
  ## ind can be a list or a matrix
  if (is.list(ind) | (is.integer(ind) & !is.matrix(ind)))  {
    n <- length(ind)
  }  else  {
    n <- ncol(ind)         # important to use col because of output from get.reciprocal.overlaps
  }
  ## fill matrix using the overlap indices
  m <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n)  {
    if (is.list(ind) | (is.integer(ind) & !is.matrix(ind)))  {
      indcol <- ind[[i]]
    }  else  {
      indcol <- ind[,i]
    }
    m[i,indcol] <- 1
  }
  m
}
  

## seems to give the same answers or not quite?
get.overlaps3.old <- function(v, mat)  {
  bg1 <- v[1]
  ed1 <- v[2]
  bg2 <- mat[,1]
  ed2 <- mat[,2]
  ind1 <- which((bg2 >= bg1) & (bg2 <= ed1))
  ind2 <- which((ed2 >= bg1) & (ed2 <= ed1))
  ind3 <- which((bg2 <= bg1) & (ed2 >= bg1))
  ind12 <- union(ind1, ind2)
  ind <- union(ind12, ind3)
  ind
}

get.overlaps3 <- function(v, mat)  {
  bg1 <- as.numeric(as.character(v[1]))
  ed1 <- as.numeric(as.character(v[2]))
  bg2 <- as.numeric(as.character(mat[,1]))
  ed2 <- as.numeric(as.character(mat[,2]))
  ind1 <- which((bg2 >= bg1) & (bg2 <= ed1))
  ind2 <- which((ed2 >= bg1) & (ed2 <= ed1))
  ind3 <- which((bg2 <= bg1) & (ed2 >= bg1))
  ind4 <- which((bg2 >= ed1) & (ed2 <= ed1))
  ind12 <- union(ind1, ind2)
  ind13 <- union(ind12, ind3)
  ind <- union(ind13, ind4)
  ind
}

#apply(mat1, 1, get.overlaps3, mat=mat2)

### get percentage overlap of v with mat, what sum and percentage of mat is captured
## v = vector
## mat = matrix
get.percentage.overlap <- function(v, mat)  {
  bg <- as.numeric(as.character(v[1]))
  ed <- as.numeric(as.character(v[2]))
  bgv <- as.numeric(as.character(mat[,1]))
  edv <- as.numeric(as.character(mat[,2]))
  ## slot v into mat, separate for bg and ed
  v1 <- sort(c(bg, bgv))
  v2 <- sort(c(ed, edv))
  indbg <- which(v1==bg)[1]
  inded <- which(v2==ed)[1]
  if (indbg<inded)  {
    ## collate pieces
    smat <- mat[indbg:(inded-1),]  
    mreg <- sum(smat[,2] - smat[,1] + 1)               ## centre region
    if (indbg>1)  {  
      lreg <- mat[indbg-1,2] - bg + 1                  ## left border region that is overlapped
      tlreg <- mat[indbg-1,2] - mat[indbg-1,1] + 1     ## total left border region
    }  else  {
      lreg <- tlreg <- 0
    }
    if (inded <= nrow(mat))  {
      rreg <- ed - mat[inded,1] + 1                    ## right border region that is overlapped
      trreg <- mat[inded,2] - mat[inded,1] + 1         ## total left border region
    }  else  {
      rreg <- trreg <- 0
    }
    ## if no overlap
    if (lreg<0)  {
      lreg <- tlreg <- 0
    }
    if (rreg<0)  {
      rreg <- trreg <- 0
    }
    oreg <- lreg + mreg + rreg
    treg <- mreg + tlreg + trreg
    p <- round(100 * oreg / treg, 2)
  }  else  {
    oreg <- p <- 0
  }
  ## total regions
  list(r=oreg, p=p)
}

### get percentage overlap of mat1 with mat2, what sum and percentage of mat2 is captured
## mat1 = matrix
## mat2 = matrix
get.percentage.overlap.chr <- function(mat1, mat2)  {
  all <- NULL
  for (i in 1:24)  {
    print(i)
    indchr1 <- which(mat1[,1]==i)
    indchr2 <- which(mat2[,1]==i)
    p <- apply(mat1[indchr1,2:3], 1, function(v) get.percentage.overlap(v, mat2[indchr2,2:3]))
    all <- c(all, p)
  }
  all
}

### get region of mat that overlaps v
## v = vector
## mat = matrix
get.overlapping.region <- function(v, mat)  {
  bg <- as.numeric(as.character(v[1]))
  ed <- as.numeric(as.character(v[2]))
  bgv <- as.numeric(as.character(mat[,1]))
  edv <- as.numeric(as.character(mat[,2]))
  ## slot v into mat, separate for bg and ed
  v1 <- sort(c(bg, bgv))
  v2 <- sort(c(ed, edv))
  indbg <- which(v1==bg)[1]
  inded <- which(v2==ed)[1]
  if (indbg<inded)  {
    ## collate pieces
    smat <- mat[indbg:(inded-1),]                    ## centre region
    if (indbg>1)  {  
      lreg <- c(bg, mat[indbg-1,2])                  ## left border region that is overlapped
    }  else  {
      lreg <- c(0,0)
    }
    if (inded <= nrow(mat))  {
      rreg <- c(mat[inded,1], ed)                    ## right border region that is overlapped
    }  else  {
      rreg <- c(0,0)
    }
    ## if no overlap
    if (lreg[1]>lreg[2])  {
      lreg <- c(0,0)
    }
    if (rreg[1]>rreg[2])  {
      rreg <- c(0,0)
    }
    areg <- rbind(lreg, smat, rreg)
  }  else  {
    areg <- c(0,0)
  }
  ## remove zero regions
  if (!is.vector(areg))  {
    indrm <- apply(areg, 1, function(v) all(v==c(0,0)))
    areg <- areg[!indrm,]
  }
  areg
}


### get percentage overlap of mat1 with mat2, what sum and percentage of mat2 is captured
## mat1 = matrix
## mat2 = matrix
get.overlapping.region.chr <- function(mat1, mat2)  {
  all <- NULL
  for (i in 1:24)  {
    print(i)
    indchr1 <- which(mat1[,1]==i)
    indchr2 <- which(mat2[,1]==i)
    p <- apply(mat1[indchr1,2:3], 1, function(v) get.overlapping.region(v, mat2[indchr2,2:3]))
    all <- c(all, p)
  }
  all
}


get.overlaps3.thresh <- function(v, mat, thresh)  {
  bg1 <- v[1]
  ed1 <- v[2]
  bg2 <- mat[,1]
  ed2 <- mat[,2]
  ind1 <- which((bg2 >= bg1) & (bg2 <= ed1))
  ind2 <- which((ed2 >= bg1) & (ed2 <= ed1))
  ind3 <- which((bg2 <= bg1) & (ed2 >= bg1))
  ind4 <- which(abs((bg2-bg1)-(ed2-ed1)) <= thresh)      ## additional threshold for size
  ind12 <- union(ind1, ind2)
  ind13 <- union(ind12, ind3)
  ind <- intersect(ind13, ind4)
  ind
}

# mat <- rbind(c(10,20),c(11,21),c(12,30))
# mat2 <- rbind(c(15,20),c(16,21),c(22,30))
# apply(mat1, 1, get.overlaps3.thresh, mat=mat2, thresh=2)


### within one matrix
get.overlaps4 <- function(mat)  {
  bg <- mat[,1]
  ed <- mat[,2]
  over <- NULL
  for (i in 1:length(bg))  {
    ind1 <- which((bg >= bg[i]) & (bg <= ed[i]))
    ind2 <- which((ed >= bg[i]) & (ed <= ed[i]))
    ind3 <- which((bg <= bg[i]) & (ed >= bg[i]))
    ind12 <- union(ind1, ind2)
    ind <- union(ind12, ind3)
    if (length(ind)==0)
      ind <- NA
    over[i] <- list(ind)
  }
  over
}

### within one matrix
### store result in matrix instead of list as in get.overlaps4
get.overlaps5 <- function(mat)  {
  n <- nrow(mat)
  bg <- mat[,1]
  ed <- mat[,2]
  over <- matrix(0, nrow=n, ncol=n)
  diag(over) <- 1
  for (i in 1:length(bg))  {
    ind1 <- which((bg >= bg[i]) & (bg <= ed[i]))
    ind2 <- which((ed >= bg[i]) & (ed <= ed[i]))
    ind3 <- which((bg <= bg[i]) & (ed >= bg[i]))
    ind12 <- union(ind1, ind2)
    ind <- union(ind12, ind3)
    if (length(ind)!=0)  {
      over[i,ind] <- 1
    }
  }
  over
}

### get reads from mat that are within v
get.including.reads <- function(v, mat,thresh)  {
  bg1 <- v[1]
  ed1 <- v[2]
  bg2 <- mat[,1]
  ed2 <- mat[,2]
  #ind <- which((bg2 >= bg1) & (ed2 <= ed1))
  ind <- which(abs(bg2-bg1)<thresh & abs(ed2-ed1)<thresh)
  ind
}

#### Gets connected components using matrix multiplication
## mat = matrix that contains which entries overlap

get.conn.comp <- function(mat)  {
  mat <- t(apply(mat, 1, as.numeric))
  m <- diag(x = 1, nrow=nrow(mat), ncol=ncol(mat))
  n <- nrow(mat)
  if (any((m-mat)!=0))  {
    nmat <- mat%*%mat
    nmat[which(nmat>0)] <- 1
    mat2 <- mat
    ## repeat matrix multiplication until matrix remains unchanged
    while(any((nmat-mat2)!=0))  {
      mat2 <- nmat
      nmat <- nmat%*%mat
      nmat[which(nmat>0)] <- 1
    }
  }
  else  {
    nmat <- m
  }
  fmat <- nmat+t(nmat)
  fmat[fmat!=0] <- 1
  fmat
}


###### New version of merging SV calls without matrix multiplication ######
###### but with connected blocks by position ###### 

#### Gets blocks clustered by position using physical RP depth 
## dat = matrix with chr, start, end

get.blocks <- function(dat)  {
  ## sort by start and end positions, add code 1 for start and -1 for end (=flag), add ID for RP (=del)
  bg <- data.frame(dat[,1:2], rep(1, nrow(dat)), 1:nrow(dat))
  ed <- data.frame(dat[,c(1,3)], rep(-1, nrow(dat)), 1:nrow(dat))
  names(bg) <- names(ed) <- c("chr","pos", "flag", "del")
  pos <- rbind(bg, ed)
  ord <- order(pos[,1], pos[,2])
  spos <- pos[ord,]
  ## cumulative sum of physical RP depth
  depth <- cumsum(spos[,3])
  dpos <- cbind(spos, depth)
  dpos
}

#### Gets clusters as a list of row indices
## dpos = output from get.blocks, i.e. matrix with chr, pos, flag, del

cluster.blocks <- function(dpos)  {
  depth <- dpos[,5]
  ## uses physical RP depth==0 to distinguish clusters, ind2 is the end index
  ind2 <- which(depth==0)
  k <- length(ind2)-1
  ## ind1 is the start index
  if (k>0)  {                         ## edited on 12-01-2010
    ind1 <- c(1, ind2[1:k]+1)
  }  else  {
    ind1 <- 1
  }
  ## combine indices for each cluster in a list (can also happen to be a matrix)
  if (length(ind1)!=1)  {
    cluster <- apply(cbind(ind1, ind2), 1, function(v) v[1]:v[2])
  }  else  {
    cluster <- list(ind1:ind2)
  }
  cluster
}

#### Merges calls that have reciprocal overlaps using connected components
## m = matrix with chr, start, end

#get.reciprocal.conn.comp <- function(m, thresh, d=500)  {
get.reciprocal.conn.comp <- function(m, thresh)  {
  ## get indices for overlaps
  if (thresh==-1)  {
    indlist <- apply(m[,-1], 1, get.overlaps3, mat=m[,-1])    ## added on 2/3/2010
  }  else  {
    #indlist <- apply(m[,-1], 1, get.reciprocal.overlaps, mat=m[,-1], thresh, d)
    indlist <- apply(m[,-1], 1, get.reciprocal.overlaps, mat=m[,-1], thresh)
  }
  if (is.matrix(indlist) | (!is.matrix(indlist) & length(indlist) != length(unlist(indlist))))  {
    ## get overlap matrix
    indover.mat <- fill.overlap.matrix(indlist)
    ## get connected components
    overlaps <- get.conn.comp(indover.mat)
    clust <- NULL
    ## merge overlapping calls
    for (i in 1:nrow(overlaps))  {
      ind <- which(overlaps[i,]>0)
      if (length(ind)>1)  {
        clustlines <- merge.lines.basic(m, ind, pcols=c(2,2,3,3), outer=FALSE)
      }  else  {
        clustlines <- c(m[ind,1], m[ind,2], m[ind,2], m[ind,3], m[ind,3])
      }
      clust <- rbind(clust, clustlines)
    }
    ## remove duplicate entries
    v <- unique(unlist(apply(clust, 1, paste, collapse="\t")))
    mat <- vecstr2matstr(v, "\t")
  }  else  {
    mat <- cbind(m[,1], m[,2], m[,2], m[,3], m[,3])
  }
  mat
}

get.reciprocal.conn.comp.ann <- function(m, thresh, pcols=c(2,2,3,3), d=500)  {
  ## get indices for overlaps
  if (thresh==-1)  {
    indlist <- apply(m[,-1], 1, get.overlaps3, mat=m[,-1])    ## added on 2/3/2010
  }  else  {
    indlist <- apply(m[,-1], 1, get.reciprocal.overlaps, mat=m[,-1], thresh, d)
  }
  if (is.matrix(indlist) | (!is.matrix(indlist) & length(indlist) != length(unlist(indlist))))  {
    ## get overlap matrix
    indover.mat <- fill.overlap.matrix(indlist)
    ## get connected components
    overlaps <- get.conn.comp(indover.mat)
    clust <- NULL
    ## merge overlapping calls
    for (i in 1:nrow(overlaps))  {
      ind <- which(overlaps[i,]>0)
      if (length(ind)>1)  {
        clustlines <- merge.lines.basic(m[min(pcols):max(pcols)], ind, pcols=pcols, outer=FALSE)
      }  else  {
        clustlines <- c(m[ind,1], m[ind,2], m[ind,2], m[ind,3], m[ind,3])
      }
      clust <- rbind(clust, clustlines)
    }
    ## remove duplicate entries
    v <- unique(unlist(apply(clust, 1, paste, collapse="\t")))
    mat <- vecstr2matstr(v, "\t")
  }  else  {
    mat <- cbind(m[,1], m[,2], m[,2], m[,3], m[,3])
  }
  mat
}

#### Pairs start and end within a cluster, uses output from get.blocks and cluster.blocks 
## dpos = output from get.blocks
## cluster = output from cluster.blocks

pair.sv.calls <- function(dpos, cluster)  {
  mdat <- merge(dpos[cluster,], dpos[cluster,], by.x=4, by.y=4)
  indrm <- which(mdat[,3]==mdat[,7] | mdat[,3]>mdat[,7])
  mdat[-indrm,]
}

#### Re-formats singleton calls that are not within a cluster 
## m = output from get.blocks

pos2del <- function(m)  {
  c(m[1,1], m[1,2], m[1,2], m[2,2], m[2,2])
}

##### Define region for clusters that have a reciprocal overlap
get.reciprocal.overlap.within.cluster <- function(dpos, cluster, thresh)  {
  ## pair up start and end
  sdat <- pair.sv.calls(dpos, cluster)
  m <- sdat[,c(2,3,7)]
  ## get reciprocal overlap
  if (nrow(m)!=0)  {
    mclust <- get.reciprocal.conn.comp(m, thresh)
  }  else  {
    mclust <- rep(NA, 5)
  }
  mclust
}

#### Merges clusters using reciprocal overlap and connected components
#### and recombines with unclustered singletons
## dpos = output from get.blocks
## cluster = output from cluster.blocks

#merge.cluster.blocks <- function(dpos, cluster, thresh=0.5, d=500) {      ## thresh=-1 for any bp overlap
merge.cluster.blocks <- function(dpos, cluster, thresh=0.5) {              ## thresh=-1 for any bp overlap
  if (is.matrix(cluster))  {
    n <- apply(cluster, 2, length)
  }  else  {
    n <- sapply(cluster, length)
  }
  ## merge clusters
  inds <- which(n>2)                 ## inds <- inds[5520:5530]     inds <- 9891
##   amclust <- arclust <- NULL
##   for (i in inds)  {
##     print(i)
##     ## pair up start and end
##     sdat <- pair.sv.calls(dpos, cluster[[i]])
##     m <- sdat[,c(2,3,7)]
##     ## get reciprocal overlap
##     if (nrow(m)!=0)  {
##       #mclust <- get.reciprocal.conn.comp(m, thresh, d)
##       mclust <- get.reciprocal.conn.comp(m, thresh)
##       print(mclust)
##     }  else  {
##       mclust <- NULL
##     }
##     amclust <- rbind(amclust, mclust)  
##   }
  if (length(inds) > 0)  {
    amclust <- t(sapply(inds, function(v) get.reciprocal.overlap.within.cluster(dpos, cluster[[v]], thresh)))
  }  else  {
    amclust <- NULL                          ## added 27/2/2012
  }
  ## prepare rest of calls
  indr <- which(n==2)
  ## check numbers
  if (is.matrix(cluster))  {
    test <- (length(inds) + length(indr)) == ncol(cluster)
  }  else  {
    test <- (length(inds) + length(indr)) == length(cluster)
  }
  if (test==FALSE)
    print(paste("all clusters considered =", test))
  ## collect rest
  arclust <- NULL
  if (length(indr)!=0)  {
    if (is.matrix(cluster))  {
      arclust <- t(apply(cluster[,indr], 2, function(v) pos2del(dpos[v,])))
    }  else  {
      arclust <- t(sapply(cluster[indr], function(v) pos2del(dpos[v,])))
    }
  }
  ## combine with rest of calls
  allclust <- rbind(amclust, arclust)
  allclust <- t(apply(allclust, 1, as.numeric))
  ord <- order(allclust[,1], allclust[,2], allclust[,3], allclust[,4])
  allclust <- allclust[ord,]
  ## remove NAs
  indna <- is.na(allclust[,1])
  allclust <- allclust[!indna,]
  allclust
}


### Function that converts output from the cluster Perl script into
### matrix for cluster positions and into a list for the types such
### as FF=1, FR=2, RF=4, RR=8

pairs2tablist <- function(dat)  {
  typ <- sapply(dat, function(s) {
    types.str <- strsplit(s,"types ")[[1]][2]
    strsplit(types.str," ")[[1]] })
  tab <- sapply(dat, function(s) {
    table.str <- strsplit(s,"types")[[1]][1]
    strsplit(table.str," ")[[1]]  })
  mat <- matrix(unlist(tab), ncol=13, nrow=length(dat), byrow=T)
  list(mat=mat, typ=typ)
}

editPairs <- function(dat)  {
  tablist <- pairs2tablist(dat)
  typ <- tablist$typ
  mat <- tablist$mat
  ind <- c(1:3,5,7:9,11,13)
  mat <- mat[,ind]
  mat <- t(apply(mat, 1, gsub, pattern="(", replacement="", extended=F))
  mat <- t(apply(mat, 1, gsub, pattern=")", replacement="", extended=F))
  mat1 <- t(apply(mat[,-ncol(mat)], 1, as.integer))
  mat2 <- as.real(mat[,ncol(mat)])
  mat <- data.frame(mat1,mat2)
  list(typ=typ, mat=mat)
}


### Use the following commands if more than one chromosome

chrpairs2tablist <- function(dat)  {
  typ <- sapply(dat, function(s) {
    types.str <- strsplit(s,"types ")[[1]][2]
    strsplit(types.str," ")[[1]] })
  tab <- sapply(dat, function(s) {
    table.str <- strsplit(s,"types")[[1]][1]
    strsplit(table.str," ")[[1]]  })
  mat <- data.frame(matrix(unlist(tab), ncol=14, nrow=length(dat), byrow=T))
  list(mat=mat, typ=typ)
}

editChrPairs <- function(dat)  {
  tablist <- chrpairs2tablist(dat)
  typ <- tablist$typ
  mat <- tablist$mat
  ind <- c(1:4,6,8:10,12,14)
  mat <- mat[,ind]
  mat <- t(apply(mat, 1, gsub, pattern="(", replacement="", extended=F))
  mat <- data.frame(t(apply(mat, 1, gsub, pattern=")", replacement="", extended=F)))
  list(typ=typ, mat=mat)
}


### Function that converts output from the cluster Perl script into
### matrix for cluster positions and into a list for singlets such
### as 64 and 130 and lib origin

singlets2tablist <- function(dat)  {
  lib <- sapply(dat, function(s) {
    libs.str <- strsplit(s,"libs ")[[1]][2]
    strsplit(libs.str," ")[[1]] })
  typ <- sapply(dat, function(s) {
    str <- strsplit(s,"types ")[[1]][2]
    types.str <- strsplit(str,"libs")[[1]][1]
    strsplit(types.str, " ")})
  tab <- sapply(dat, function(s) {
    table.str <- strsplit(s,"types")[[1]][1]
    strsplit(table.str," ")[[1]]  })
  mat <- matrix(as.numeric(unlist(tab)), ncol=4, nrow=length(dat), byrow=T)
  list(mat=mat, typ=typ, lib=lib)
}

### get types from list typ from function singlets2tablist and put them in a matrix
### using types 1,4,8,32,64,130

get.type <- function(tlist, type=c(1,4,8,32,64,130))  {
  tlist <- lapply(tlist, gsub, pattern="\\(", replacement="", perl=T)
  tlist <- lapply(tlist, gsub, pattern="\\)", replacement="", perl=T)
  tlist <- lapply(tlist, function(v) matrix(as.numeric(v), byrow=T, ncol=2, nrow=length(v)/2))
  n <- length(tlist)
  m <- length(type)
  type.mat <- matrix(0, ncol=m, nrow=n)
  for (i in 1:n)  {
    for (j in 1:m)  {
      ind <- which(tlist[[i]][,1]==type[j])
      if (length(ind)!=0)  {
        type.mat[i,j] <- tlist[[i]][ind,2]
      }
    }
  }
  type.mat
}

### get libraries from list lib from function singlets2tablist and put them in a matrix
### using libraries from 1:16

get.lib <- function(tlist,num=16)  {
  tlist <- lapply(tlist, gsub, pattern="\\(", replacement="", perl=T)
  tlist <- lapply(tlist, gsub, pattern="\\)", replacement="", perl=T)
  tlist <- lapply(tlist, function(v) matrix(as.numeric(v), byrow=T, ncol=2, nrow=length(v)/2))
  n <- length(tlist)
  type.mat <- matrix(0, ncol=num, nrow=n)
  for (i in 1:n)  {
    for (j in 1:num)  {
      ind <- which(tlist[[i]][,1]==j)
      if (length(ind)!=0)  {
        type.mat[i,j] <- tlist[[i]][ind,2]
      }
    }
  }
  type.mat
}

### Convert perl singlet cluster into matrix format
# dat = tab delimited string from readLines, output from perl script 
#       scripts/project_1000_genomes/startCluster.pl

singletcluster2matrix <- function(dat)  {
  clustlist <- singlets2tablist(dat)
  clust <- clustlist$mat
  tmat <- get.type(clustlist$typ)
  libmat <- get.lib(clustlist$lib)
  fclust <- cbind(clust, tmat, libmat)
  fclust
}



### Compare clusters by indicating whether the cluster start/end falls
### within a threshold
               
cluster.diff <- function(x)  {
  d <- abs(as.real(x[1:(length(x)-1)]) - as.real(x[2:length(x)]))
  d
}

compare.clusters <- function(x1, x2, thresh)  {
  ord1 <- order(x1[,1])                  # start of cluster
  ord2 <- order(x2[,1])
  sx1 <- x1[ord1,]
  sx2 <- x2[ord2,]
  m1 <- cbind(sx1, rep(1,nrow(x1)))
  m2 <- cbind(sx2, rep(2,nrow(x2)))
  m <- rbind(m1, m2)
  ord <- order(m[,1])
  m <- m[ord,]
  d <- c(cluster.diff(m[,1]), NA)        # calculate difference of cluster start
  ind <- d < thresh
  #indn <- which(d < thresh)             # maybe better without adding TRUE to next element
  #ind[indn+1] <- TRUE;
  #if (is.na(ind[length(ind)]))
  #  ind[length(ind)] <- FALSE            # FALSE for last element
  mat <- data.frame(m, ind)
  mat
}


### Missed clusters within ELAND/MAQ

cluster.analysis <- function(clust)  {
  clind <- which(clust[,11])
  sclind <- which(clust[clind,10] != clust[clind+1,10])
  sclindlen <- which(clust[clind,10] != clust[clind+1,10] & clust[clind+1,11]==FALSE)
  sclind1 <- which(clust[clind,10] == clust[clind+1,10] & clust[clind,10]==1)
  sclind2 <- which(clust[clind,10] == clust[clind+1,10] & clust[clind,10]==2)

  clustfind <- sort(union(clind[sclind],clind[sclind]+1))
  find1 <- sort(union(clind[sclind1],clind[sclind1]+1))
  find2 <- sort(union(clind[sclind2],clind[sclind2]+1))
  find <- sort(c(find1,find2))

  cluster.intersection <- clust[clustfind,]
  tomerge <- clust[find,]
  n <- length(sclindlen)
  n1 <- length(find1)
  n2 <- length(find2)
  list(cluster.intersection=cluster.intersection, tomerge=tomerge, n=n, n1=n1, n2=n2)
}

### Get dimension statistics
get.dim <- function(dat1, dat2, clust)  {
  m <- cluster.analysis(clust)
  n <- m$n
  n1 <- m$n1
  n2 <- m$n2
  dimens <- c(length(dat1), length(dat2), n, n1, n2)
  dimens
}

### Compare ELAND with MAQ, ELAND with ELAND and MAQ with MAQ
compare.all.clusters <- function(dat1, dat2, mat1, mat2)  {
  dimat <- matrix(0, nrow=5, ncol=5)
  thresh <- c(1,10,100,200,300)
  for (i in 1:length(thresh))  {
    clust <- compare.clusters(mat1, mat2, thresh=thresh[i])
    dimat[i,] <- get.dim(dat1, dat2, clust)
  }
  list(clust=clust, dimat=dimat)
}

### Compare ELAND and MAQ with other data sources
# percentage in respect to the 1st matrix
get.all.clusters <- function(mat1, mat2, full="no")  {
  n1 <- n2 <- NULL
  thresh <- c(10,50,100,200)
  for (i in 1:length(thresh))  {
    clust <- compare.clusters(mat1, mat2, thresh=thresh[i])
    ind <- which(clust[,11])
    ind <- sort(union(ind,ind+1))
    if (full=="yes")  {
      print(clust[ind,1:10])
    }
    n1[i] <- round(length(which(clust[ind,10]==1)),2)                        # not quite accurate, better to put in smaller matrix as mat1
    n2[i] <- round(length(which(clust[ind,10]==1))/nrow(mat1),2)*100
    n <- rbind(n1,n2)
  }
  print(n)
}


### Merge clusters from different sources, subfunction that merges two
### types below a distance threshold

merge.types <- function(typ1, typ2)  {
  typ1 <- unlist(strsplit(as.character(typ1), " "))
  typ1 <- gsub("(", "", typ1, extended=F)
  typ1 <- gsub(")", "", typ1, extended=F)
  typ2 <- unlist(strsplit(as.character(typ2), " "))
  typ2 <- gsub("(", "", typ2, extended=F)
  typ2 <- gsub(")", "", typ2, extended=F)
  typ <- c(typ1,typ2)
  even <- seq(2, length(typ), 2)
  odd <- seq(1, length(typ), 2)
  typmat <- cbind(typ[odd],typ[even])

  mat <- matrix(0, ncol=2, nrow=6)
  mat[,1] <- c(1,2,4,8,18,"NA")
  for (i in 1:nrow(typmat))  {
    ind <- which(as.character(typmat[i,1])==as.character(mat[,1]))
    mat[ind,2] <- as.numeric(as.character(mat[ind,2])) + as.numeric(as.character(typmat[i,2]))
    if (length(ind)==0)
      print(paste("Type =", typmat[i,1]))
  }
  # remove zero entries
  ind <- as.numeric(as.character(mat[,2])) != 0
  nmat <- as.matrix(mat[ind,])
  if (ncol(nmat)==1)
    nmat <- t(nmat)
  # print(nrow(nmat))
  
  v <- NULL
  for (i in 1:nrow(nmat))  {
    v <- paste(v, nmat[i,1], " (", nmat[i,2], ") ", sep="")
  }
  n <- nchar(v)
  v <- substr(v, 1, n-1)
  v
}

merge.types2 <- function(typs, ind)  {
  typs <- as.character(typs[ind])
  typ <- unlist(lapply(typs, strsplit, " "))
  ind1 <- which(typ=="1")
  ind2 <- which(typ=="2")
  ind4 <- which(typ=="4")
  ind8 <- which(typ=="8")
  ind18 <- which(typ=="18")
  typ <- gsub("(", "", typ, extended=F)
  typ <- as.numeric(as.character(gsub(")", "", typ, extended=F)))
  st <- c(1,2,4,8,18)
  ss <- c(sum(typ[ind1+1]), sum(typ[ind2+1]), sum(typ[ind4+1]), sum(typ[ind8+1]), sum(typ[ind18+1]))
  ind <- which(ss==0)
  if (length(ind)!=0)  {
    st <- st[-ind]
    ss <- ss[-ind]
  } 
  v <- paste(st, " (", ss, ")", sep="", collapse=" ")
  v
}

### If types and libs are stored in a matrix
merge.types.libs <- function(dat, ind)  {
  ndat <- t(apply(as.matrix(dat[ind,]), 1, function(v) as.numeric(as.character(v))))
  fdat <- apply(ndat, 2, sum)
  fdat
}

### Merge clusters from different sources, subfunction that merges two
### lines below a distance threshold, uses merge.types function

merge.lines.old <- function(mat, i)  {
  bg1 <- min(as.numeric(as.character(mat[i,1])), as.numeric(as.character(mat[i+1,1])))
  ed1 <- max(as.numeric(as.character(mat[i,3])), as.numeric(as.character(mat[i+1,3])))
  bg2 <- min(as.numeric(as.character(mat[i,5])), as.numeric(as.character(mat[i+1,5])))
  ed2 <- max(as.numeric(as.character(mat[i,7])), as.numeric(as.character(mat[i+1,7])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  svd <- bg2-ed1-1
  n <- as.numeric(as.character(mat[i,8])) + as.numeric(as.character(mat[i+1,8]))
  avlen <- round((as.numeric(as.character(mat[i,8])) * as.numeric(as.character(mat[i,9])) + as.numeric(as.character(mat[i+1,8])) * as.numeric(as.character(mat[i+1,9]))) / n, 1)
  type <- merge.types(mat[i,10], mat[i+1,10])
    #paste(as.character(mat[i,10]), as.character(mat[i+1,10]), sep=",")
  #orig1 <- as.numeric(as.character(unlist(strsplit(as.character(mat[i,10]),","))))
  #orig2 <- as.numeric(as.character(unlist(strsplit(as.character(mat[i+1,10]),","))))
  orig1 <- as.numeric(as.character(unlist(strsplit(as.character(mat[i,11]),","))))
  orig2 <- as.numeric(as.character(unlist(strsplit(as.character(mat[i+1,11]),","))))
  ored <- sort(union(orig1,orig2))
  orig <- paste(ored, sep=",",collapse=",")
  newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, n, avlen, type, orig)
  newline
}

#### Merges overlapping calls by determining inner and outer boundaries of merged calls
## mat = matrix of chr, start, end OR chr, outer start, inner start, inner end, outer end
## outer = if TRUE  then chr, outer start, outer end are returned,
##         if FALSE then chr, outer start, inner start, inner end, outer end are returned

merge.lines.basic <- function(mat, ind, pcols=c(2,2,3,3), outer=TRUE)  {
  chr <- min(as.numeric(as.character(mat[ind,1])))
  if (outer==TRUE)  {
    bg <- min(as.numeric(as.character(mat[ind,pcols[1]])))
    ed <- max(as.numeric(as.character(mat[ind,pcols[4]])))
    newline <- c(chr, bg, ed)
  }  else  {
    bg.outer <- min(as.numeric(as.character(mat[ind,pcols[1]])))
    bg.inner <- max(as.numeric(as.character(mat[ind,pcols[2]])))
    ed.inner <- min(as.numeric(as.character(mat[ind,pcols[3]])))
    ed.outer <- max(as.numeric(as.character(mat[ind,pcols[4]])))
    newline <- c(chr, bg.outer, bg.inner, ed.inner, ed.outer)
  }
  newline
}

### To merge the calls from all centres for the Trio-CEU
# mat = chr1, start2, end3, rpdel4, rddel5, rddup6, origin7

merge.lines.basic.ann.old <- function(mat, ind, pcols=c(2,2,3,3), tcols=c(4:6), ocols=7, outer=TRUE, secmerge=FALSE)  {
  chr <- min(as.numeric(as.character(mat[ind,1])))
  type <- apply(mat[ind,tcols], 2, sum)
  if (secmerge == FALSE) {
    origin <- paste(sort(unique(as.numeric(as.character(mat[ind,ocols])))), collapse=",",sep="")
  }  else  {
    origin <- paste(sort(unique(as.numeric(as.character(unlist(vecstr2matstr(mat[ind,ocols],",")))))), collapse=",",sep="")
  }
  if (outer==TRUE)  {
    bg <- min(as.numeric(as.character(mat[ind,pcols[1]])))
    ed <- max(as.numeric(as.character(mat[ind,pcols[2]])))
    newline <- c(chr, bg, ed, type, origin)
  }  else  {
    bg.outer <- min(as.numeric(as.character(mat[ind,pcols[1]])))
    bg.inner <- max(as.numeric(as.character(mat[ind,pcols[2]])))
    ed.inner <- min(as.numeric(as.character(mat[ind,pcols[3]])))
    ed.outer <- max(as.numeric(as.character(mat[ind,pcols[4]])))
    newline <- c(chr, bg.outer, bg.inner, ed.inner, ed.outer, type, origin)
  }
  newline
}

### To merge the calls from all centres for the Trio-CEU and keep track of origin for intersection and union
# mat = chr1, start2, end3, rpdel4, rddel5, rddup6, origin7

merge.lines.basic.ann <- function(m,ind,pcols=c(2,2,3,3),tcols=c(4:6),ocols=7,outer=TRUE,secmerge=FALSE,track=FALSE,exclude=NULL,num=TRUE)  {
  chr <- min(as.numeric(as.character(m[ind,1])))
  if (!is.null(tcols)) {
    type <- apply(m[ind,tcols], 2, sum)
  }  else  {
    type <- NULL
  }
  if (secmerge == FALSE) {
    if (num==TRUE)  {  ## sort numeric
      origin <- paste(sort(unique(as.numeric(as.character(m[ind,ocols])))), collapse=",",sep="")
    }  else  {         ## sort string
      origin <- paste(sort(unique(as.character(m[ind,ocols]))), collapse=",",sep="")
    }    
  }  else  {
    if (num==TRUE)  {  ## sort numeric
      origin <- paste(sort(unique(as.numeric(as.character(unlist(vecstr2matstr(m[ind,ocols],",")))))), collapse=",",sep="")
    }  else  {         ## sort string
      origin <- paste(sort(unique(as.character(unlist(vecstr2matstr(m[ind,ocols],","))))), collapse=",",sep="")
    }    
  }
  if (outer==TRUE)  {    ## supply only 2-tuple for pcols
    ## find begin and end
    if (is.null(exclude))  {
      bg <- min(as.numeric(as.character(m[ind,pcols[1]])))
      ed <- max(as.numeric(as.character(m[ind,pcols[2]])))
    }  else  {    ## exclude contains vector of rows to remove (e.g. RD)
      indrm <- as.vector(unlist(sapply(exclude, function(v) which(v==m[ind,ocols]))))
      if (length(indrm)!=0 & length(m[ind,1][-indrm])!=0)  {
        bg <- min(as.numeric(as.character(m[ind,pcols[1]][-indrm])))
        ed <- max(as.numeric(as.character(m[ind,pcols[2]][-indrm])))
      }  else  {
        bg <- min(as.numeric(as.character(m[ind,pcols[1]])))
        ed <- max(as.numeric(as.character(m[ind,pcols[2]])))
      }
    }
    ## keep track of intersection and union origin
    if (track==TRUE)  {  
      indm <- which(m[ind,pcols[1]]==bg)
      obg <- paste(m[ind,ocols][indm], collapse=",", sep="")
      indm <- which(m[ind,pcols[2]]==ed)
      oed <- paste(m[ind,ocols][indm], collapse=",", sep="")
      newline <- c(chr, bg, ed, type, origin, obg, oed)
    }  else  {
      newline <- c(chr, bg, ed, type, origin)
    }      
  }  else  {             ## supply 4-tuple for pcols
    ## find begin and end interval 
    if (is.null(exclude))  {
      bg.outer <- min(as.numeric(as.character(m[ind,pcols[1]])))
      bg.inner <- max(as.numeric(as.character(m[ind,pcols[2]])))
      ed.inner <- min(as.numeric(as.character(m[ind,pcols[3]])))
      ed.outer <- max(as.numeric(as.character(m[ind,pcols[4]])))
    }  else  {
      indrm <- as.vector(unlist(sapply(exclude, function(v) which(v==m[ind,ocols]))))
      if (length(indrm)!=0 & length(m[ind,1][-indrm])!=0)  {
        bg.outer <- min(as.numeric(as.character(m[ind,pcols[1]][-indrm])))
        bg.inner <- max(as.numeric(as.character(m[ind,pcols[2]][-indrm])))
        ed.inner <- min(as.numeric(as.character(m[ind,pcols[3]][-indrm])))
        ed.outer <- max(as.numeric(as.character(m[ind,pcols[4]][-indrm])))
      }  else  {
        bg.outer <- min(as.numeric(as.character(m[ind,pcols[1]])))
        bg.inner <- max(as.numeric(as.character(m[ind,pcols[2]])))
        ed.inner <- min(as.numeric(as.character(m[ind,pcols[3]])))
        ed.outer <- max(as.numeric(as.character(m[ind,pcols[4]])))
      }
    }
    ## keep track of intersection and union origin
    if (track==TRUE)  {
      indm <- which(m[ind,pcols[1]]==bg.outer)
      obg.outer <- paste(sort(unique(as.numeric(as.character(m[ind,ocols][indm])))), collapse=",", sep="")
      indm <- which(m[ind,pcols[2]]==bg.inner)
      obg.inner <- paste(sort(unique(as.numeric(as.character(m[ind,ocols][indm])))), collapse=",", sep="")
      indm <- which(m[ind,pcols[3]]==ed.inner)
      oed.inner <- paste(sort(unique(as.numeric(as.character(m[ind,ocols][indm])))), collapse=",", sep="")
      indm <- which(m[ind,pcols[4]]==ed.outer)
      oed.outer <- paste(sort(unique(as.numeric(as.character(m[ind,ocols][indm])))), collapse=",", sep="")
      newline <- c(chr, bg.outer, bg.inner, ed.inner, ed.outer, type, origin, obg.outer, obg.inner, oed.inner, oed.outer)
    }  else  {
      newline <- c(chr, bg.outer, bg.inner, ed.inner, ed.outer, type, origin)
    }
  }
  newline
}

merge.lines <- function(mat, i, j)  {
  bg1 <- min(as.numeric(as.character(mat[i,1])), as.numeric(as.character(mat[j,1])))
  ed1 <- max(as.numeric(as.character(mat[i,3])), as.numeric(as.character(mat[j,3])))
  bg2 <- min(as.numeric(as.character(mat[i,5])), as.numeric(as.character(mat[j,5])))
  ed2 <- max(as.numeric(as.character(mat[i,7])), as.numeric(as.character(mat[j,7])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  svd <- bg2-ed1-1
  n <- as.numeric(as.character(mat[i,8])) + as.numeric(as.character(mat[j,8]))
  avlen <- round((as.numeric(as.character(mat[i,8])) * as.numeric(as.character(mat[i,9])) + as.numeric(as.character(mat[j,8])) * as.numeric(as.character(mat[j,9]))) / n, 1)
  type <- merge.types(mat[i,10], mat[j,10])
  orig1 <- as.numeric(as.character(unlist(strsplit(as.character(mat[i,11]),","))))
  orig2 <- as.numeric(as.character(unlist(strsplit(as.character(mat[j,11]),","))))
  ored <- sort(union(orig1,orig2))
  orig <- paste(ored, sep=",",collapse=",")
  newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, n, avlen, type, orig)
  newline
}

merge.lines2 <- function(mat, ind)  {
  bg1 <- min(as.numeric(as.character(mat[ind,1])))
  ed1 <- max(as.numeric(as.character(mat[ind,3])))
  bg2 <- min(as.numeric(as.character(mat[ind,5])))
  ed2 <- max(as.numeric(as.character(mat[ind,7])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  svd <- bg2-ed1-1
  n <- sum(as.numeric(as.character(mat[ind,8])))
  avlen <- round(sum((as.numeric(as.character(mat[ind,8])) * as.numeric(as.character(mat[ind,9])))) / n, 1)
  type <- merge.types2(mat[,10], ind)
  ## origin
  ored <- unlist(strsplit(as.character(mat[ind,11]),","))
  ored <- sort(unique(as.numeric(ored)))
  orig <- paste(ored, sep=",",collapse=",")
  ## collate
  newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, n, avlen, type, orig)
  newline
}

## No origin is only difference to merge.lines2
merge.lines3 <- function(mat, ind)  {
  bg1 <- min(as.numeric(as.character(mat[ind,1])))
  ed1 <- max(as.numeric(as.character(mat[ind,3])))
  bg2 <- min(as.numeric(as.character(mat[ind,5])))
  ed2 <- max(as.numeric(as.character(mat[ind,7])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  svd <- bg2-ed1-1
  n <- sum(as.numeric(as.character(mat[ind,8])))
  avlen <- round(sum(as.numeric(as.character(mat[ind,8])) * as.numeric(as.character(mat[ind,9]))) / n, 1)
  type <- merge.types2(as.character(mat[,10]), ind) 
  newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, n, avlen, type)
  newline
}

## Similar to merge.lines2, but types and libs are stored in matrix
merge.lines4 <- function(mat, ind, lib=TRUE)  {
  bg1 <- min(as.numeric(as.character(mat[ind,1])))
  ed1 <- max(as.numeric(as.character(mat[ind,3])))
  bg2 <- min(as.numeric(as.character(mat[ind,5])))
  ed2 <- max(as.numeric(as.character(mat[ind,7])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  svd <- bg2-ed1-1
  n <- sum(as.numeric(as.character(mat[ind,8])))
  avlen <- round(sum((as.numeric(as.character(mat[ind,8])) * as.numeric(as.character(mat[ind,9])))) / n, 1)
  type <- merge.types.libs(mat[,10:15], ind)
  if (lib==TRUE)  {
    ## origin 
    orig <- merge.types.libs(mat[,16:31], ind)
    ## collate
    newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, n, avlen, type, orig)
  }  else  {
    newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, n, avlen, type)
  }
  newline
}

## Similar to merge.lines2, but types and libs are stored in matrix
## for mouse data - chr17, rmove chr
merge.lines4a <- function(mat, ind)  {
  bg1 <- min(as.numeric(as.character(mat[ind,1])))
  ed1 <- max(as.numeric(as.character(mat[ind,3])))
  bg2 <- min(as.numeric(as.character(mat[ind,5])))
  ed2 <- max(as.numeric(as.character(mat[ind,7])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  svd <- bg2-ed1-1
  n <- sum(as.numeric(as.character(mat[ind,8])))
  avlen <- round(sum((as.numeric(as.character(mat[ind,8])) * as.numeric(as.character(mat[ind,9])))) / n, 1)
  type <- merge.types.libs(mat[,10:17], ind)
  newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, n, avlen, type)
  newline
}

## Similar to merge.lines2, but types and libs are stored in matrix - for 1000 genomes Trio-CEU (July 2008 freeze)
## Use also for Melanesian
merge.lines5 <- function(mat, ind, tcol=c(12:19), lcol=20, lib=TRUE, origin=FALSE, sample.id=FALSE)  {
  bg1 <- min(as.numeric(as.character(mat$start1[ind])))
  ed1 <- max(as.numeric(as.character(mat$end1[ind])))
  bg2 <- min(as.numeric(as.character(mat$start2[ind])))
  ed2 <- max(as.numeric(as.character(mat$end2[ind])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  sv <- bg2-ed1
  n <- sum(as.numeric(as.character(mat$n[ind])))
  avlen <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avlen[ind])))) / n, 1)
  avqual <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avqual[ind])))) / n, 1)
  ## type
  if (!is.null(tcol))  {
    type <- merge.types.libs(mat[,tcol], ind)
  } 
  ## origin
  if (lib)  {
    orig <- merge.types.libs(mat[,lcol:ncol(mat)], ind)
    newline <- c(bg1, cd1, ed1, sv, bg2, cd2, ed2, n, avlen, avqual, type, orig)
  }  else if (origin)  {
    v <- mat[ind,lcol]
    if (length(grep(",", v))==0)  {
      if (sample.id==FALSE)  {
        orig <- paste(sort(unique(as.numeric(as.character(v)))), collapse=",")
      }  else  {
        orig <- paste(sort(unique(as.character(v))), collapse=",")
      }
    }  else  {
      if (sample.id==FALSE)  {
        orig <- paste(sort(unique(as.numeric(as.character(unlist(vecstr2matstr(v,",")))))), collapse=",",sep="")
      }  else  {
        orig <- paste(sort(unique(as.character(unlist(vecstr2matstr(v,","))))), collapse=",",sep="")
      }
    }
    newline <- c(bg1, cd1, ed1, sv, bg2, cd2, ed2, n, avlen, avqual, type, orig)
  } else if (is.null(tcol))  {
    newline <- c(bg1, cd1, ed1, sv, bg2, cd2, ed2, n, avlen, avqual, orig)
  } else  {
    newline <- c(bg1, cd1, ed1, sv, bg2, cd2, ed2, n, avlen, avqual, type)
  }
  newline
}

merge.lines5.simple <- function(mat, ind, pcol=c(2,8), ocol=7, num=TRUE)  {
  ## get intersection of deletion
  bg <- max(as.numeric(as.character(mat[ind,pcol[1]])))
  ed <- min(as.numeric(as.character(mat[ind,pcol[2]])))
  n <- sum(as.numeric(as.character(mat$n[ind])))
  avlen <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avlen[ind])))) / n, 1)
  avqual <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avqual[ind])))) / n, 1)
  ## origin
  v <- mat[ind,ocol]
  if (num==TRUE)  {
    if (length(grep(",", v))==0)  {
      orig <- paste(sort(unique(as.numeric(as.character(v)))), collapse=",")
    }  else  {
      orig <- paste(sort(unique(as.numeric(as.character(unlist(vecstr2liststr(v,",")))))), collapse=",",sep="")
    }
  }  else  {
    if (length(grep(",", v))==0)  {
      orig <- paste(sort(unique(as.character(v))), collapse=",")
    }  else  {
      orig <- paste(sort(unique(as.character(unlist(vecstr2liststr(v,","))))), collapse=",",sep="")
    }
  }
  newline <- c(bg, ed, n, avlen, avqual, orig)
  as.character(newline)
}

merge.lines.main <- function(mat, ind)  {
  bg1 <- min(as.numeric(as.character(mat$start1[ind])))
  ed1 <- max(as.numeric(as.character(mat$end1[ind])))
  bg2 <- min(as.numeric(as.character(mat$start2[ind])))
  ed2 <- max(as.numeric(as.character(mat$end2[ind])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  sv <- bg2-ed1
  n <- sum(as.numeric(as.character(mat$n[ind])))
  avlen <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avlen[ind])))) / n, 1)
  avqual <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avqual[ind])))) / n, 1)
  newline <- c(bg1, cd1, ed1, sv, bg2, cd2, ed2, n, avlen, avqual)
  newline
}

merge.lines.simple <- function(mat, ind)  {
  bg <- max(as.numeric(as.character(mat$start[ind])))
  ed <- min(as.numeric(as.character(mat$end[ind])))
  n <- sum(as.numeric(as.character(mat$n[ind])))
  avlen <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avlen[ind])))) / n, 1)
  avqual <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avqual[ind])))) / n, 1)
  newline <- c(bg, ed, n, avlen, avqual)
  newline
}

## Similar to merge.lines5, but types and libs are stored in matrix - for 1000 genomes Trio-CEU (July 2008 freeze)
## Use for inversions !!!
## c("chr","start1","bpl1","bpl2","end1","sv","start2","bpr1","bpr2","end2","n","avlen","avqual","t1","t2","t4","t8","t18","t32")
merge.lines.inv <- function(mat, ind, tcol=c(14:19), lib=TRUE)  {
  bg1 <- min(as.numeric(as.character(mat$start1[ind])))
  ed1 <- max(as.numeric(as.character(mat$end1[ind])))
  bg2 <- min(as.numeric(as.character(mat$start2[ind])))
  ed2 <- max(as.numeric(as.character(mat$end2[ind])))
  bpl1 <- max(as.numeric(as.character(mat$bpl1[ind])))
  bpl2 <- min(as.numeric(as.character(mat$bpl2[ind])))
  bpr1 <- max(as.numeric(as.character(mat$bpr1[ind])))
  bpr2 <- min(as.numeric(as.character(mat$bpr2[ind])))
  svd <- bg2-ed1
  n <- sum(as.numeric(as.character(mat$n[ind])))
  avlen <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avlen[ind])))) / n, 1)
  avqual <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avqual[ind])))) / n, 1)
  ## type
  type <- merge.types.libs(mat[,tcol], ind)
  ## origin
  if (lib)  {
    orig <- merge.types.libs(mat[,20:ncol(mat)], ind)
    newline <- c(bg1, bpl1, bpl2, ed1, svd, bg2, bpr1, bpr2, ed2, n, avlen, avqual, type, orig)
  }  else  {
    newline <- c(bg1, bpl1, bpl2, ed1, svd, bg2, bpr1, bpr2, ed2, n, avlen, avqual, type)
  }
  newline
}

## Similar to merge.lines2, but types and libs are stored in matrix - for 1000 genomes Trio-CEU (July 2008 freeze)
## Insertions
merge.lines6 <- function(mat, ind)  {
  bg1 <- min(as.numeric(as.character(mat$start1[ind])))
  ed1 <- max(as.numeric(as.character(mat$end1[ind])))
  bg2 <- min(as.numeric(as.character(mat$start2[ind])))
  ed2 <- max(as.numeric(as.character(mat$end2[ind])))
  cd1 <- ed1-bg1+1
  cd2 <- ed2-bg2+1
  svd <- bg2-ed1-1
  sn1 <- sum(as.numeric(as.character(mat$n1[ind])))
  sn2 <- sum(as.numeric(as.character(mat$n2[ind])))
  n1 <- as.numeric(as.character(mat$n1[ind]))
  n2 <- as.numeric(as.character(mat$n2[ind]))
  m <- apply(cbind(n1,n2), 1, sum)
  n <- sum(n1,n2)
  avqual <- round(sum(m * as.numeric(as.character(mat$avqual[ind]))) / n, 1)  
  ## type
  type <- merge.types.libs(mat[,12:19], ind)
  ## origin 
  orig <- merge.types.libs(mat[,20:25], ind)
  ## collate
  newline <- c(bg1, cd1, ed1, svd, bg2, cd2, ed2, sn1, sn2, avqual, type, orig)
  newline
}

## Similar to merge.lines2, but types and libs are stored in matrix - for 1000 genomes Trio-CEU (July 2008 freeze)
## Translocations
merge.lines7 <- function(mat, ind)  {
  bg <- min(as.numeric(as.character(mat$start[ind])))
  ed <- max(as.numeric(as.character(mat$end[ind])))
  cd <- ed-bg+1
  n <- sum(as.numeric(as.character(mat$n[ind])))
  avqual <- round(sum((as.numeric(as.character(mat$n[ind])) * as.numeric(as.character(mat$avqual[ind])))) / n, 1)
  ## type
  type <- merge.types.libs(mat[,7:14], ind)
  ## origin 
  orig <- merge.types.libs(mat[,15:20], ind)
  ## collate
  newline <- c(bg, cd, ed, n, avqual, type, orig)
  newline
}

## Merging lines for translocation clusters
merge.lines.trans <- function(mat)  {
  bg1 <- min(as.numeric(as.character(mat$start1)))
  ed1 <- max(as.numeric(as.character(mat$start1)))
  bg2 <- min(as.numeric(as.character(mat$start2)))
  ed2 <- max(as.numeric(as.character(mat$start2)))
  cd1 <- ed1-bg1
  cd2 <- ed2-bg2
  #n <- length(ind)
  n <- nrow(mat)
  qu <- round(mean(mat$quality1),1)
  newline <- c(mat$chr1[1], bg1, ed1, cd1, mat$chr2[1], bg2, ed2, cd2, n, qu)
  newline
}

## Merging lines for translocation clusters
merge.lines.trans.lib <- function(mat)  {
  bg1 <- min(as.numeric(as.character(mat$start1)))
  ed1 <- max(as.numeric(as.character(mat$end1)))
  bg2 <- min(as.numeric(as.character(mat$start2)))
  ed2 <- max(as.numeric(as.character(mat$end2)))
  cd1 <- ed1-bg1
  cd2 <- ed2-bg2
  n <- sum(mat$n)
  qu <- round(mean(mat$quality),1)
  newline <- c(mat$chr1[1], bg1, ed1, cd1, mat$chr2[1], bg2, ed2, cd2, n, qu)
  newline
}

### Merge clusters from different sources, uses function merge.lines
### and therefore also merge.types

merge.clusters <- function(mat, thresh)  {
  mat <- as.matrix(mat)
  d1.left <- cluster.diff(mat[,1])
  d2.left <- cluster.diff(mat[,3])
  d1.right <- cluster.diff(mat[,5])
  d2.right <- cluster.diff(mat[,7])
  ind1 <- d1.left < thresh | d2.left < thresh
  ind2 <- d1.right < thresh | d2.right < thresh
  ind <- ind1 & ind2
  fclust <- matrix(0, nrow=nrow(mat), ncol=ncol(mat))
  ## check 1st line
  if (!ind[1])  {                               
    fclust[1,] <- mat[1,]
    j <- 1
    k <- 2
  }  else  {
    j <- 0
    k <- 1
  }
  ## for loop from 1st line if TRUE or from 2nd line if 1st line is FALSE
  for (i in k:(nrow(mat)-1))  {
    j <- j+1
    if (ind[i])  {     ## if TRUE
      if (ind[i-1])
        j <- j-1
      ## copy new cluster into next line of mat
      # print(paste("i =",i))
      fclust[j,] <- merge.lines(mat, i, i+1)
      mat[i+1,] <- fclust[j,]
    }  else  {         ## if FALSE
      if (!ind[i-1])  {        ## if previous ind == FALSE, meaning no previous cluster, starting new cluster
        fclust[j,] <- mat[i,]
      }  else 
      j <- j-1
    } 
  }
  ## last row
  if (!ind[nrow(mat)-1])  {   
    j <- j+1
    fclust[j,] <- mat[nrow(mat),]
  }
  ## edit, change to integer and remove rows with zeros
  indi <- c(10,11)
  fmat <- apply(fclust[,-indi], 2, as.integer)
  ind0 <- apply(fmat, 1, function(v) any(v!=0))
  fmat <- data.frame(fmat, fclust[,10:11])
  fmat <- fmat[ind0,]
  names(fmat) <- c("bg1", "int1", "ed1", "sv", "bg2", "int2", "ed2", "n", "avlen", "type", "origin")
  fmat
}


### Merge overlapping clusters
# m = chr, start, end
merge.clusters.basic <- function(m,cols=c(2,3),thresh=NULL,pcols=c(2,2,3,3),tcols=c(4:6),ocols=7,outer=TRUE,
                                 ann=FALSE,secmerge=FALSE,num=TRUE, d=500)
{
  if (is.vector(m))  {
    mat <- t(as.matrix(as.character(m[,c(1,pcols,tcols,ocols)])))
  }  else  {
    mp <- data.frame(t(apply(m[,cols], 1, as.character)))
    mp <- data.frame(t(apply(mp, 1, as.numeric)))
    if (is.null(thresh)) {
      ## 1bp overlap
      indlist <- apply(mp, 1, get.overlaps3, mat=mp)  
    }  else  {
      ## e.g. 50% overlap
      indlist <- apply(mp, 1, get.reciprocal.overlaps, mat=mp, thresh=thresh, d)
    }
    if (is.list(indlist))  {
      if (ann==FALSE)  {
        clust <- lapply(indlist, merge.lines.basic, mat=m, pcols=pcols, outer=outer)
      }  else  {
        clust <- lapply(indlist, merge.lines.basic.ann, mat=m, pcols=pcols,tcols=tcols,ocols=ocols,outer=outer,secmerge=secmerge,num=num)
      }    
      v <- unique(unlist(lapply(clust, paste, collapse="\t")))
      mat <- vecstr2matstr(v, "\t")
    }  else  {
      mat <- m
    }
    ord <- order(as.numeric(as.character(mat[,2])))
    mat <- t(apply(mat[ord,], 1, as.character))
  }
  mat
}

## for merging calls across centres
merge.clusters.conn.comp <- function(m,cols=c(2,3),thresh=NULL,pcols=c(2,2,3,3),tcols=c(4:6),ocols=7,
                                     outer=TRUE,ann=FALSE,secmerge=FALSE,track=FALSE,exclude=NULL,num=TRUE, d=500)  {
  if (is.vector(m))  {
    if (track==FALSE)  {
      mat <- t(as.matrix(as.character(m[,c(1,pcols,tcols,ocols)])))
    }  else  {
      mat <- t(as.matrix(as.character(m[,c(1,pcols,tcols,rep(ocols,5))])))
    }
  }  else  {
    if (nrow(m) > 1)  {
      mp <- data.frame(t(apply(m[,cols], 1, as.character)))
      mp <- data.frame(t(apply(mp, 1, as.numeric)))
      if (is.null(thresh)) {
        ## 1bp overlap
        indlist <- apply(mp, 1, get.overlaps3, mat=mp)  
        indover.mat <- fill.overlap.matrix(indlist)
      }  else  {
        ## e.g. 50% overlap
        indlist <- apply(mp, 1, get.reciprocal.overlaps, mat=mp, thresh=thresh, d)
        indover.mat <- fill.overlap.matrix(indlist)
      }
      ## get clusters via connected components
      overlaps <- get.conn.comp(indover.mat)
      ## merge lines that overlap and include rest
      clust <- NULL
      for (i in 1:nrow(overlaps))  {
        ind <- which(overlaps[i,]>0)
        ##if (length(ind)!=0)  {
        if (ann==FALSE)  {
          clustlines <- merge.lines.basic(m, ind, pcols=pcols, outer=outer)
          clust <- rbind(clust, clustlines)
        }  else  {
          clustlines <- merge.lines.basic.ann(m,ind,pcols=pcols,tcols=tcols,ocols=ocols,
                                              outer=outer,secmerge=secmerge,track=track,exclude=exclude,num=num)
          clust <- rbind(clust, clustlines)
        }    
      }
      v <- unique(unlist(apply(clust, 1, paste, collapse="\t")))
      mat <- vecstr2matstr(v, "\t")
      if (is.matrix(mat)) {
        ord <- order(as.numeric(as.character(mat[,2])))
        mat <- t(apply(as.matrix(mat[ord,]), 1, as.character))
      }  else  {
        mat <- as.character(mat)
      }
    }  else  {
      if (track==FALSE & ann==TRUE)  {
        mat <- as.character(m[,c(1,pcols,tcols,ocols)])
      }  else if (track==TRUE & ann==TRUE) {
        mat <- as.character(m[,c(1,pcols,tcols,rep(ocols,5))])
      }  else if (track==FALSE) {
        mat <- as.character(m[,c(1,pcols)])
      }
    }
  }
  mat
}
 
## for merging clusters from scripts/project_1000_genomes/pairCluster1000Chr.pl output
merge.clusters.conn.comp2 <- function(m, cols=c(2,8), thresh=0.5, tcol=c(12:19), lcol=20, lib=TRUE, origin=FALSE,
                                      sample.id=FALSE, num=FALSE, simple=FALSE, d=500)  {
  if (nrow(m) > 1)  {
    mp <- data.frame(t(apply(m[,cols], 1, as.character)))
    mp <- data.frame(t(apply(mp, 1, as.numeric)))
    if (is.null(thresh)) {
      ## 1bp overlap
      indlist <- apply(mp, 1, get.overlaps3, mat=mp)  
      indover.mat <- fill.overlap.matrix(indlist)
    }  else  {
      ## e.g. 50% overlap
      indlist <- apply(mp, 1, get.reciprocal.overlaps, mat=mp, thresh=thresh, d)
      indover.mat <- fill.overlap.matrix(indlist)
    }
    ## get clusters via connected components
    overlaps <- get.conn.comp(indover.mat)
    ## merge lines that overlap and include rest
    clust <- NULL
    for (i in 1:nrow(overlaps))  {
      ind <- which(overlaps[i,]>0)
      if (length(ind)>0)  {
        if (simple==FALSE)  {
          clustlines <- merge.lines5(m, ind, tcol=tcol, lcol=lcol, lib=lib, origin=origin, sample.id)
        }  else  {
          clustlines <- merge.lines5.simple(m, ind, pcol=cols, ocol=lcol, num)
        }
        clust <- rbind(clust, clustlines)
      } 
    }     
    v <- unique(unlist(apply(clust, 1, paste, collapse="\t")))
    mat <- vecstr2matstr(v, "\t")
    chr <- rep(m[1,1], nrow(mat))    ## add chr
    mat <- data.frame(chr, mat)           ## combine chr with dat
    ord <- order(as.numeric(as.character(mat[,2])))
    #mat <- t(apply(mat[ord,], 1, as.numeric))
    mat <- mat[ord,]
  }  else  {
    mat <- t(as.matrix(c(as.character(m[1,1]), as.character(m[,-1]))))
  }
  mat
}

### Loop through chr
merge.clusters.basic.chr <- function(dat,cols=c(2,3),thresh=NULL,pcols=c(2,2,3,3),tcols=c(4:6),ocols=7,
                                     outer=TRUE,ann=FALSE,conn=FALSE,secmerge=FALSE,track=FALSE,exclude=NULL,num=TRUE)  {
  m <- NULL
  for (i in 1:24)  {
    print(i)
    indchr <- which(dat[,1]==i)
    if (length(indchr) > 1)  {
      if (conn==FALSE)  {
##         mm <- merge.clusters.basic(dat[indchr,],cols=cols,thresh=thresh,pcols=pcols,tcols=tcols,ocols=ocols,
##                                    outer=outer,ann=ann, secmerge=secmerge, track=track, exclude=exclude, num=num)
        mm <- merge.clusters.basic(dat[indchr,],cols=cols,thresh=thresh,pcols=pcols,tcols=tcols,ocols=ocols,
                                   outer=outer,ann=ann, secmerge=secmerge, num=num)
        m <- rbind(m, mm)
      }  else  {
        mm <- merge.clusters.conn.comp(dat[indchr,],cols=cols,thresh=thresh,pcols=pcols,tcols=tcols,ocols=ocols,
                                       outer=outer,ann=ann,secmerge=secmerge,track=track,exclude=exclude, num=num)
        m <- rbind(m, mm)
      }
    } else if (length(indchr) == 1)  {
      if (track==FALSE & ann==TRUE)  {
        m <- rbind(m, c(as.character(dat[indchr,c(1,pcols,tcols)]), as.character(dat[indchr,ocols])))
      }  else if (track==TRUE & ann==TRUE) {
        m <- rbind(m, c(as.character(dat[indchr,c(1,pcols,tcols)]), as.character(dat[indchr,rep(ocols,5)])))
      }  else if (ann==FALSE & outer==FALSE) {
        m <- rbind(m, as.character(dat[indchr,c(1,pcols)]))
      }  else if (ann==FALSE & outer==TRUE) {
        m <- rbind(m, as.character(dat[indchr,]))
      }
    }
  }
  m
}

### Loop through chr for clusters generated by scripts/project_1000_genomes/pairCluster1000Chr.pl
merge.clusters.conn.comp.by.chr <- function(dat, cols=c(2,8), thresh=0.5, lib=TRUE)  {
  m <- NULL
  for (i in 1:24)  {
    indchr <- which(dat[,1]==i)
    if (length(indchr)!=0)  {
        mm <- merge.clusters.conn.comp2(dat[indchr,],cols=cols,thresh=thresh,lib=lib)
        m <- rbind(m, mm)
    }
  }
  m
}


### Merge clusters from different sources, uses the overlap2 function, the function merge.lines
### and therefore also merge.types
# m = dataframe in format from pairCluster.pl script
merge.clusters2 <- function(m)  {
  pos <- cbind(as.numeric(as.character(m[,1])), as.numeric(as.character(m[,7])))
  indlist <- apply(pos, 1, get.overlaps3, mat=pos)
  clust <- lapply(indlist, merge.lines2, mat=m)
  v <- unique(unlist(lapply(clust, paste, collapse="\t")))
  v
}

## No origin is only difference to merge.lines2
merge.clusters3 <- function(m, c1=1, c2=7)  {
  pos <- cbind(as.numeric(as.character(m[,c1])), as.numeric(as.character(m[,c2])))
  indlist <- apply(pos, 1, get.overlaps3, mat=pos)
  clust <- lapply(indlist, merge.lines3, mat=m)
  v <- unique(unlist(lapply(clust, paste, collapse="\t")))
  v
}

## Similar to merge.lines2, but types and libs are stored in matrix
merge.clusters4 <- function(m, c1=1, c2=7, lib=TRUE)  {
  pos <- cbind(as.numeric(as.character(m[,c1])), as.numeric(as.character(m[,c2])))
  indlist <- apply(pos, 1, get.overlaps3, mat=pos)
  ## merge lines that overlap
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  #clust <- lapply(indlist[indl], merge.lines4, mat=m, lib)
  clust <- lapply(indlist[indl], merge.lines4, mat=m, lib)
  v <- unique(unlist(lapply(clust, paste, collapse="\t")))
  clust <- vecstr2matstr(v, "\t")
  clust <- as.data.frame(clust)
  names(clust) <- names(m)
  ## recombine with lines that were not merged
  fmat <- rbind(m[-indl,], clust)
  ord <- order(as.numeric(as.character(fmat[,c1])))
  fmat <- fmat[ord,]
  fmat
}

## Similar to merge.lines2, but types and libs are stored in matrix
## for mouse data
merge.clusters4a <- function(m, c1=1, c2=7)  {
  pos <- cbind(as.numeric(as.character(m[,c1])), as.numeric(as.character(m[,c2])))
  indlist <- apply(pos, 1, get.overlaps3, mat=pos)
  ## merge lines that overlap
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    clust <- lapply(indlist[indl], merge.lines4a, mat=m)
    v <- unique(unlist(lapply(clust, paste, collapse="\t")))
    clust <- vecstr2matstr(v, "\t")
    clust <- as.data.frame(clust)
    names(clust) <- names(m)
    ## recombine with lines that were not merged
    fmat <- rbind(m[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,c1])))
    fmat <- fmat[ord,]
  }
  else {
    fmat <- m
  }
  fmat
}

## Similar to merge.lines2, but types and libs are stored in matrix - for 1000 genomes (July 2008 freeze)
## use also for Melanesian and Trio-YRI and inversions (invs=TRUE)
merge.clusters5 <- function(m, c1=2, c2=8, thresh=NULL, lib=TRUE, recip=TRUE, tcol=c(12:19), lcol=20, invs=FALSE, origin=FALSE, d=500)  {
  pos <- cbind(as.numeric(as.character(m[,c1])), as.numeric(as.character(m[,c2])))
  if (recip==FALSE)  {
    if (is.null(thresh))  {
      indlist <- apply(pos, 1, get.overlaps3, mat=pos)
    }  else  {
      indlist <- apply(pos, 1, get.overlaps3.thresh, mat=pos, thresh=thresh)
    }
  }  else  {
    indlist <- apply(pos, 1, get.reciprocal.overlaps, mat=pos, thresh=thresh, d)
  }
  ## merge lines that overlap
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    if (invs==FALSE)  {
      clust <- lapply(indlist[indl], merge.lines5, mat=m, tcol=tcol, lcol=lcol, lib=lib, origin=origin)
    }  else  {
      clust <- lapply(indlist[indl], merge.lines.inv, mat=m, tcol=tcol, lib=lib)
    }
    v <- unique(unlist(lapply(clust, paste, collapse=" ")))
    clust <- vecstr2matstr(v, " ")
    clust <- as.data.frame(clust)
    chr <- rep(m$chr[1], nrow(clust))
    clust <- data.frame(chr, clust)
    names(clust) <- names(m)
    ## recombine with lines that were not merged
    fmat <- rbind(m[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,c1])))
    fmat <- fmat[ord,]
  }  else  {
    fmat <- m
  }
  fmat
}

## Similar to merge.lines2, but no types and no libs - for 1000 genomes (Aug 2010 freeze)
merge.clusters.main <- function(m, c1=2, c2=8)  {
  pos <- cbind(as.numeric(as.character(m[,c1])), as.numeric(as.character(m[,c2])))
  indlist <- apply(pos, 1, get.overlaps3, mat=pos)
  ## merge lines that overlap
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    clust <- lapply(indlist[indl], merge.lines.main, mat=m)
    v <- unique(unlist(lapply(clust, paste, collapse=" ")))
    clust <- vecstr2matstr(v, " ")
    clust <- as.data.frame(clust)
    chr <- rep(m$chr[1], nrow(clust))
    clust <- data.frame(chr, clust)
    names(clust) <- names(m)
    ## recombine with lines that were not merged
    fmat <- rbind(m[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,c1])))
    fmat <- fmat[ord,]
  }  else  {
    fmat <- m
  }  
  fmat
}

## Similar to merge.lines2, but types and libs are stored in matrix - for 1000 genomes (July 2008 freeze)
## Insertions
merge.clusters6 <- function(m, c1=2, c2=8)  {
  pos <- cbind(as.numeric(as.character(m[,c1])), as.numeric(as.character(m[,c2])))
  indlist <- apply(pos, 1, get.overlaps3, mat=pos)
  ## merge lines that overlap
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    clust <- lapply(indlist[indl], merge.lines6, mat=m)
    v <- unique(unlist(lapply(clust, paste, collapse=" ")))
    clust <- vecstr2matstr(v, " ")
    clust <- as.data.frame(clust)
    chr <- rep(m$chr[1], nrow(clust))
    clust <- data.frame(chr, clust)
    names(clust) <- names(m)
    ## recombine with lines that were not merged
    fmat <- rbind(m[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,c1])))
    fmat <- fmat[ord,]
  }  else  {
    fmat <- m
  }
  fmat
}
## Similar to merge.lines2, but types and libs are stored in matrix - for 1000 genomes (July 2008 freeze)
## Translocations
merge.clusters7 <- function(m, c1=2, c2=8)  {
  pos <- cbind(as.numeric(as.character(m[,c1])), as.numeric(as.character(m[,c2])))
  indlist <- apply(pos, 1, get.overlaps3, mat=pos)
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    ## merge lines that overlap
    clust <- lapply(indlist[indl], merge.lines7, mat=m)
    v <- unique(unlist(lapply(clust, paste, collapse=" ")))
    clust <- vecstr2matstr(v, " ")
    clust <- as.data.frame(clust)
    chr <- rep(m$chr[1], nrow(clust))
    clust <- data.frame(chr, clust)
    names(clust) <- names(m)
    ## recombine with lines that were not merged
    fmat <- rbind(m[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,c1])))
    fmat <- fmat[ord,]
  }  else  {
    fmat <- m
  }
  fmat
}

### Remove duplicates (identical start positions on both chrs)
### and retain the pair with the highest quality

rm.clusters.trans.0.old <- function(mat)  {
  s1 <- mat$start1
  s2 <- mat$start2
  retain <- NULL
  indlist <- rm.identical.pos(s1, s2)
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    i <- 1
    while (i <= length(indl))  {
      inddup <- indlist[indl][[i]]
      (indmax <- which(mat$quality1[inddup]==max(mat$quality1[inddup]))[1])
      (retain <- rbind(retain,mat[inddup,][indmax,]))
      i <- i+length(inddup)
    }
    rmat <- rbind(mat[-indl,], retain)
  }  else  {
    rmat <- mat
  }
  rmat
}

rm.clusters.trans.0.old <- function(mat)  {
  s1 <- mat$start1
  s2 <- mat$start2
  retain <- NULL
  indlist <- rm.identical.pos(s1, s2)
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    qualind <- sapply(indlist[indl], function(v) which.max(mat$quality1[v]))
    i <- 1
    while (i <= length(indl))  {
      inddup <- indlist[indl][[i]]      
      (retain <- rbind(retain,mat[inddup,][qualind[i],]))
      i <- i+length(inddup)
    }
    rmat <- rbind(mat[-indl,], retain)
  }  else  {
    rmat <- mat
  }
  rmat
}

rm.clusters.trans.0 <- function(mat)  {
  s1 <- mat$start1
  s2 <- mat$start2
  retain <- NULL
  indlist <- rm.identical.pos(s1, s2)
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    qualind <- unique(sapply(indlist[indl], function(v) v[which.max(mat$quality1[v])]))
    retain <- mat[qualind,]
    rmat <- rbind(mat[-indl,], retain)
  }  else  {
    rmat <- mat
  }
  rmat
}

### Merge clusters of start positions on chrA combined with associated clusters on chrB
# start1 = start positions on chrA
# start2 = start positions on chrB
# thresh.low = lower threshold for distance between start positions (0 means same positions can occur several times)
# thresh.up = upper threshold for distance between start positions (8*max(sd) = 250 for g1k)

merge.clusters.trans <- function(mat, thresh.low=0, thresh.up=250)  {
  start1 <- mat$start1
  start2 <- mat$start2
  indlist <- get.overlap.0(start1, start2, thresh.low, thresh.up)
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    ## merge lines that overlap and skip those that don't
    clust <- lapply(indlist[indl], merge.lines.trans, mat=mat)
    v <- unique(unlist(lapply(clust, paste, collapse="\t")))
    clust <- vecstr2matstr(v, "\t")
    clust <- as.data.frame(clust)
  }
  clust
}


merge.clusters.trans2.old <- function(mat, thresh.low=0, thresh.up=250)  {
  start1 <- mat$start1
  start2 <- mat$start2
  indlist <- get.overlap.0(start1, start2, thresh.low, thresh.up)
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  clust <- NULL
  if (length(indl)!=0)  {
    overlist <- indlist[indl]
    k <- 1
    finoverlist <- overlist[1]
    for (i in 1:(length(overlist)-1))  {
      interind <- intersect(finoverlist[[k]], overlist[[i+1]])
      if (length(interind)!=0)  {
        finoverlist[k] <- list(union(finoverlist[[k]], overlist[[i+1]]))
      }  else  {
        k <- k+1
        finoverlist[k] <- overlist[i+1]
      }
    }
    ## merge lines that overlap
    clust <- lapply(finoverlist, merge.lines.trans, mat=mat)
    v <- unique(unlist(lapply(clust, paste, collapse="\t")))
    clust <- vecstr2matstr(v, "\t")
    clust <- as.data.frame(clust)
    names(clust) <- c("chr1","start1","end1","intvl1","chr2","start2","end2","intvl2","n","quality")
    ## recombine with lines that were not merged
    intvl <- rep(0, nrow(mat))
    n <- rep(1, nrow(mat))
    m <- data.frame(mat$chr1,mat$start1,mat$start1,intvl,mat$chr2,mat$start2,mat$start2,intvl,n,mat$quality2)
    names(m) <- c("chr1","start1","end1","intvl1","chr2","start2","end2","intvl2","n","quality")
    fmat <- rbind(m[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,1])),as.numeric(as.character(fmat[,5])),as.numeric(as.character(fmat[,2])),as.numeric(as.character(fmat[,6])))
    fmat <- fmat[ord,]
  }  else  {
    intvl <- rep(0, nrow(mat))
    n <- rep(1, nrow(mat))
    fmat <- data.frame(mat$chr1,mat$start1,mat$start1,intvl,mat$chr2,mat$start2,mat$start2,intvl,n,mat$quality2)
    names(fmat) <- c("chr1","start1","end1","intvl1","chr2","start2","end2","intvl2","n","quality")
  }
  fmat
}

merge.clusters.trans2 <- function(mat, thresh.low=0, thresh.up=250)  {
  start1 <- mat$start1
  start2 <- mat$start2
  indover <- get.overlap.0.mat(start1, start2, thresh.low, thresh.up)
  ## get clusters via connected components
  overlaps <- get.conn.comp(indover)
  ## merge lines that overlap and include rest
  clust <- NULL
  for (i in 1:nrow(overlaps))  {
    ind <- which(overlaps[i,]>0)
    clustlines <- merge.lines.trans(mat[ind,])
    clust <- rbind(clust, clustlines)
  }
  v <- unique(unlist(apply(clust, 1, paste, collapse="\t")))
  clust <- vecstr2matstr(v, "\t")
  clust <- as.data.frame(clust)
  names(clust) <- c("chr1","start1","end1","intvl1","chr2","start2","end2","intvl2","n","quality")
  ord <- order(as.numeric(as.numeric(as.character(clust[,2])),as.numeric(as.character(clust[,6]))))
  clust[ord,]
}


merge.clusters.trans.lib.old <- function(mat, thresh.low=0, thresh.up=250)  {
  start1 <- mat$start1
  start2 <- mat$start2
  indlist <- get.overlap.0(start1, start2, thresh.low, thresh.up)
  #indlist <- apply(mat[,2:3], 1, get.overlaps3, mat=mat[,2:3])
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  if (length(indl)!=0)  {
    ## merge lines that overlap
    clust <- lapply(indlist[indl], merge.lines.trans.lib, mat=mat)
    v <- unique(unlist(lapply(clust, paste, collapse="\t")))
    clust <- vecstr2matstr(v, "\t")
    clust <- as.data.frame(clust)
    names(clust) <- c("chr1","start1","end1","intvl1","chr2","start2","end2","intvl2","n","quality")
    ## recombine with lines that were not merged
    fmat <- rbind(mat[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,1])),as.numeric(as.character(fmat[,5])),as.numeric(as.character(fmat[,2])),as.numeric(as.character(fmat[,6])))
    fmat <- fmat[ord,]
  }  else  {
    fmat <- mat
  }
  fmat
}

merge.clusters.trans.lib.old2 <- function(mat, thresh.low=0, thresh.up=250)  {
  start1 <- mat$start1
  start2 <- mat$start2
  indlist <- get.overlap.0(start1, start2, thresh.low, thresh.up)
  indlen <- unlist(lapply(indlist, length))
  indl <- which(indlen > 1)
  overlist <- indlist[indl]
  clust <- NULL
  if (length(overlist)!=0)  {
    ## find union of overlaps
    k <- 1
    finoverlist <- overlist[1]
    for (i in 1:(length(overlist)-1))  {
      interind <- intersect(finoverlist[[k]], overlist[[i+1]])
      if (length(interind)!=0)  {
        finoverlist[k] <- list(union(finoverlist[[k]], overlist[[i+1]]))
      }  else  {
        k <- k+1
        finoverlist[k] <- overlist[i+1]
      }
    }
    ## merge lines that overlap
    clust <- lapply(finoverlist, merge.lines.trans.lib, mat=mat)
    v <- unique(unlist(lapply(clust, paste, collapse="\t")))
    clust <- vecstr2matstr(v, "\t")
    clust <- as.data.frame(clust)
    names(clust) <- c("chr1","start1","end1","intvl1","chr2","start2","end2","intvl2","n","quality")
    ## recombine with lines that were not merged
    fmat <- rbind(mat[-indl,], clust)
    ord <- order(as.numeric(as.character(fmat[,1])),as.numeric(as.character(fmat[,5])),as.numeric(as.character(fmat[,2])),as.numeric(as.character(fmat[,6])))
    fmat <- fmat[ord,]
  }  else  {
    fmat <- mat
  }
  fmat
}


merge.clusters.trans.lib <- function(mat, thresh.low=0, thresh.up=250)  {
  start1 <- mat$start1
  start2 <- mat$start2
  indover <- get.overlap.0.mat(start1, start2, thresh.low, thresh.up)
  ## get clusters via connected components
  overlaps <- get.conn.comp(indover)
  ## merge lines that overlap and include rest
  clust <- NULL
  for (i in 1:nrow(overlaps))  {
    ind <- which(overlaps[i,]>0)
    clustlines <- merge.lines.trans.lib(mat[ind,])
    clust <- rbind(clust, clustlines)
  }
  v <- unique(unlist(apply(clust, 1, paste, collapse="\t")))
  clust <- vecstr2matstr(v, "\t")
  clust <- as.data.frame(clust)
  names(clust) <- c("chr1","start1","end1","intvl1","chr2","start2","end2","intvl2","n","quality")
  ord <- order(as.numeric(as.character(clust[,1])),as.numeric(as.character(clust[,5])),as.numeric(as.character(clust[,2])),as.numeric(as.character(clust[,6])))
  clust[ord,]
}

merge.read.clusters <- function(mat, thresh1=0.75, thresh2=250, p.over=TRUE, tcol=c(12:19), d=500)  {
  mmat <- cbind(mat$start1, mat$end2)
  if(p.over)  {
    indover <- apply(mmat, 1, get.reciprocal.overlaps, mat=mmat, thresh=thresh1, d)
    indover.mat <- fill.overlap.matrix(indover)
  } else  {
    indover.mat <- get.overlap.0.mat(mmat[,1], mmat[,2], thresh.low=0, thresh.up=thresh2) 
  }
  ## get clusters via connected components
  overlaps <- get.conn.comp(indover.mat)
  ## merge lines that overlap and include rest
  clust <- NULL
  for (i in 1:nrow(overlaps))  {
    ind <- which(overlaps[i,]>0)
    #if (length(ind)!=0)  {
      clustlines <- merge.lines5(mat, ind, lib=FALSE, tcol=tcol)
      clust <- rbind(clust, clustlines)
    #}
    #clust <- rbind(clust, as.matrix(mat[i,-1]))
  }
  v <- unique(unlist(apply(clust, 1, paste, collapse="\t")))
  clust <- vecstr2matstr(v, "\t")
  chr <- rep(mat[1,1], nrow(clust))
  clust <- data.frame(chr, clust)
  names(clust) <- names(mat)
  clust
}


### Merge libraries to individuals from 1000 genomes project
# d = path name for input files
# fname = path and file name for output files

merge.individs <- function(d, fname, c2=8)  {       ## should be replaced by merge.libs !!!!!!!
  alldat <- NULL
  ## combine libraries
  for (i in 1:length(d))  {
    print(i)
    dat <- read.delim(d[i], sep=" ")
    alldat <- rbind(alldat, dat) 
  }  
  ## loop through each chr for merging
  clust <- NULL
  for (i in 1:24)  {
    ## get chr
    ind <- which(alldat[,1]==i & alldat[,9]>=2)     ## have at least 2 RPs per lib
    ## order by position
    if (length(ind)!=0)  {
      ord <- order(alldat[ind,2])
      mdat <- alldat[ind,][ord,]
      ## merge clusters
      clust <- rbind(clust, merge.clusters5(mdat, c1=2, c2=c2))
    }
  }
  write.table(clust, fname, sep=" ",quote=F, row.names=F)
}

merge.libs <- function(d, fname, c2=8, rpnum=1, thresh=NULL, lib=TRUE, recip=TRUE)  {          ## replaces merge.individs !!!!!
  alldat <- NULL
  ## combine libraries
  for (i in 1:length(d))  {
    print(i)
    dat <- read.delim(d[i], sep=" ")
    alldat <- rbind(alldat, dat) 
  }  
  ## loop through each chr for merging
  clust <- NULL
  for (i in 1:24)  {
    ## get chr
    ind <- which(alldat[,1]==i & alldat[,9]>=rpnum)                      ## have at least rpnum RPs per lib
    ## order by position
    if (length(ind)!=0)  {
      ord <- order(alldat[ind,2])
      mdat <- alldat[ind,][ord,]
      ## merge clusters
      clust <- rbind(clust, merge.clusters5(mdat, c1=2, c2=c2, thresh, lib, recip))
    }
  }
  write.table(clust, fname, sep=" ",quote=F, row.names=F)
}

## uses connected components
merge.libs.conn.comp <- function(d,fname,cols=c(2,8),rpnum=1,thresh=0.5,lib=TRUE,col.names=NULL)  {  ## replaces merge.individs !!!!!
  alldat <- NULL
  ## combine libraries
  for (i in 1:length(d))  {
    print(i)
    dat <- read.delim(d[i], header=F, sep=" ")
    alldat <- rbind(alldat, dat) 
  }  
  ## loop through each chr for merging
  clust <- NULL
  for (i in 1:24)  {
    ## get chr
    ind <- which(alldat[,1]==i & alldat[,9]>=rpnum)                      ## have at least rpnum RPs per lib
    ## order by position
    if (length(ind)!=0)  {
      ord <- order(alldat[ind,2])
      mdat <- alldat[ind,][ord,]
      names(mdat) <- col.names
      ## merge clusters
      clust <- rbind(clust, merge.clusters.conn.comp2(mdat, cols=cols, thresh=thresh, lib=lib))
    }
  }
  write.table(clust, fname, sep=" ",quote=F, row.names=F, col.names=col.names)
}


merge.libs.trans <- function(d, fname, rp.thresh=2, thresh.low=0, thresh.up=250)  {
  alldat <- NULL
  ## combine libraries
  for (i in 1:length(d))  {
    print(i)
    dat <- read.delim(d[i])
    alldat <- rbind(alldat, dat) 
  }
  fclust <- NULL
  for (i in 1:24)  {     ## loop for 1st mate
    print(i)
    for (j in 1:24)  {   ## loop for 2nd mate
      ## get chr
      ind <- which(alldat[,1]==i & alldat[,5]==j)                      ## have at least rpnum RPs per lib
      ## order by position
      if (length(ind)!=0)  {
        ord <- order(alldat[ind,2], alldat[ind,6])
        mdat <- alldat[ind,][ord,]
        ## merge clusters
        clust <- merge.clusters.trans.lib(mdat, thresh.low, thresh.up)
        fclust <- rbind(fclust, clust)
      }
    }
  }
  ## filter for RP number >= rp.thresh
  indnum <- which(as.numeric(as.character(fclust[,9])) >= rp.thresh)
  fclust <- fclust[indnum,]
  write.table(fclust, fname, sep="\t",quote=F, row.names=F)
}


### Merge individuals from 1000 genomes project 
# d = path name for input files
# fname = path and file name for output files
# summary2 =  summary statistics for libraries
# i = index for library

filter.individs <- function(d, fname, summary2, i, num=2)  {
  invs <- read.delim(d, sep=" ")
  ## filter for anomalous insert sizes depending on library
  mu <- min(summary2$ins_mean[i:(i+1)])
  sde <- min(summary2$ins_sd[i:(i+1)])
  ## filter for anomalous FR depending on library
  thresh1 <- mu + 4*sde - 35
  thresh2 <- 2*(mu + 4*sde)
  ind2 <- which(invs[,5]>=thresh1 & invs[,3]<=thresh2 & invs[,7]<=thresh2 & invs[,3]>0 & invs[,7]>0 & invs$n>=num)
  invsed <- invs[ind2,]
  invsed <- data.frame(invsed)
  write.table(invsed, fname, sep=" ",quote=F, row.names=F)
}

## for Trio-YRI
# d = directory
# fname = file name for output
# mu = minimum of library means
# sde = minimum of library sd
# num = minumum number of RPs

filter.individs2 <- function(d, fname, thresh1, thresh2, num=2)  {
  invs <- read.delim(d, sep=" ")
  ## filter for anomalous FR depending on library
  ind2 <- which(invs[,5]>=thresh1 & invs[,3]<=thresh2 & invs[,7]<=thresh2 & invs[,3]>0 & invs[,7]>0 & invs[,9]>=num)
  invsed <- invs[ind2,]
  invsed <- data.frame(invsed)
  write.table(invsed, fname, sep=" ",quote=F, row.names=F)
}

filter.invs.individs <- function(d, fname, summary2, i, num=2)  {
  invs <- read.delim(d, sep=" ")
  ## filter for anomalous insert sizes depending on library
  mu <- min(summary2$ins_mean[i:(i+1)])
  sde <- min(summary2$ins_sd[i:(i+1)])
  ## filter for anomalous FR depending on library
  thresh1 <- mu + 4*sde - 35
  thresh2 <- 2*(mu + 4*sde)
  #ind2 <- which(invs[,5]>=0 & invs[,3]<=thresh2 & invs[,7]<=thresh2 & invs[,3]>0 & invs[,7]>0 & invs$n>=num & invs$t1>1 & invs$t8>1)
  ind2 <- which(invs[,5]>=0 & invs[,3]<=thresh2 & invs[,7]<=thresh2 & invs[,3]>=-35 & invs[,7]>=-35 & invs$n>=num)
  invsed <- invs[ind2,]
  invsed <- data.frame(invsed)
  write.table(invsed, fname, sep=" ",quote=F, row.names=F)
}


### Extract read pair sequence from MAQ and ELAND alignments

extract.reads <- function(mat, loc, type) {
  if (type=="eland")  {
    indneg <- mat[,14] < 0
    bg <- mat[,7]
    ed <- mat[,21]
    bg[indneg] <- bg[indneg] + 35         # add 35 to bg if offset < 0
    ed[!indneg] <- ed[!indneg] + 35       # add 35 to ed if offset > 0
    ind <- ed - bg < 0                    # swap bg and ed if ed < bg
    nbg <- bg
    ned <- ed
    nbg[ind] <- ed[ind]
    ned[ind] <- bg[ind]
    bg <- nbg
    ed <- ned
    # swap read pair sequence for those where bg and ed were swapped
    rp1 <- as.character(mat[,3])                        
    rp2 <- as.character(mat[,17])
    rp2[ind] <- as.character(mat[ind,3]) 
    rp1[ind] <- as.character(mat[ind,17]) 
  }
  if (type=="maq")  { 
    bg <- mat[,2]
    ed <- bg + mat[,4]
  }
  fmat <- NULL
  for (i in 1:nrow(mat))  {
    ind <- any((bg[i] >= loc[,2] & bg[i] <= loc[,3])  |  (ed[i] >= loc[,2] & ed[i] <= loc[,3]))
    if (ind)  {
      if (type=="eland")  {       
        s <- paste(as.character(rp1[i]), as.character(rp2[i]), sep="")
      }
      if (type=="maq")  {       
        s <- paste(as.character(mat[i,5]), as.character(mat[i,8]), sep="")
      }      
      v <- c(bg[i], ed[i], s)
      fmat <- rbind(fmat, v)
    }
  }
  fmat12 <- apply(fmat[,-3], 2, as.integer)
  fmat <- data.frame(fmat12, fmat[,3])
  names(fmat) <- c("start","end","seq")
  fmat
}


### Compile all reads from each cluster and the appropriate chrX sequence

compile.cluster.reads <- function(mat, loc, fname)   {
  for (i in 1:nrow(loc))  {
    ind <- which(mat[,1] >= loc[i,2] & mat[,1] <= loc[i,3])
    h <- s <- f <- NULL
    for (j in ind)  {
      h <- paste(">", mat[j,1], mat[j,2])
      s <- as.character(mat[j,3])
      f <- rbind(f, h, s)
    }
    if (length(f) != 0)
      write(f, fname[i])
  }
}


##### Split chrX sequences in separate files

split.fasta <- function(fasta, fname)  {
  n <- length(fasta)/2
  h <- get.header(fasta, table=F)
  s <- get.sequence(fasta, table=F)
  for (i in 1:n)  {
    f <- rbind(h[i], s[i])
    write(f, fname[i])
  }
}


##### Get type information from paired clusters

get.eland.maq.types <- function(num, loc, aclustt, pos20)  {
  bg <- loc[num,2]
  ed <- loc[num,3]

  ind <- aclustt$start1 >= bg & aclustt$end2 <= ed
  print(aclustt[ind,])

  ind <- pos20[,4] >= bg & pos20[,5] <= ed
  print(pos20[ind,])
}


###### General overlap function ######
# table: start, end,
# add: length, id

### not working yet

find.overlaps <- function(mat1, mat2)  {
  id1 <- paste("a", c(1:nrow(mat1)), sep="")
  id2 <- paste("b", c(1:nrow(mat2)), sep="")

  #len1 <- mat1[,2] - mat1[,1] + 1
  #len2 <- mat2[,2] - mat2[,1] + 1
  
  nmat1 <- data.frame(mat1, id1)
  nmat2 <-  data.frame(mat2, id2)
  names(nmat1) <- names(nmat2) <- c("bg","ed","id")
  
  mat <- rbind(nmat1, nmat2)
  ord <- order(mat[,1],mat[,2])
  rmat <- mat[ord,]

  n <- nrow(rmat)
  ind <- rmat[1:(n-1),2] - rmat[2:n,1] >= 0     # new start < previous end
  # ind2 <- rmat[1:(n-1),2] - rmat[2:n,2] >= 0     # new end < previous end
  #ind <- ind1 | ind2
  ind <- c(ind,NA)
  print(ind)
  fmat <- data.frame(rmat, ind)

  fclust <- NULL
  ## check 1st line
  if (!ind[1])  {                               
    fclust[1] <- list(fmat[1,3])    # collect ids
    j <- 1
    k <- 2
  }  else  {
    j <- 0
    k <- 1
  }
  ## for loop from 1st line if TRUE or from 2nd line if 1st line is FALSE
  for (i in k:(nrow(fmat)-1))  {
    j <- j+1
    if (ind[i])  {     ## if TRUE
      if (ind[i-1])
        j <- j-1
      ## copy new cluster into next line of fmat
      # print(paste("i =",i))
      fclust[j,] <- merge.lines(fmat, i, i+1)
      fmat[i+1,] <- fclust[j,]
    }  else  {         ## if FALSE
      if (!ind[i-1])  {        ## if previous ind == FALSE, meaning no previous cluster, starting new cluster
        fclust[j,] <- fmat[i,]
      }  else 
      j <- j-1
    } 
  }
  ## last row
  if (!ind[nrow(fmat)-1])  {   
    j <- j+1
    fclust[j,] <- fmat[nrow(fmat),]
  }
  ## edit, change to integer and remove rows with zeros
  indi <- c(10,11)
  ffmat <- apply(fclust[,-indi], 2, as.integer)
  ind0 <- apply(ffmat, 1, function(v) any(v!=0))
  ffmat <- data.frame(ffmat, fclust[,10:11])
  ffmat <- ffmat[ind0,]
  names(ffmat) <- c("bg1", "int1", "ed1", "sv", "bg2", "int2", "ed2", "n", "avlen", "type", "origin")
  ffmat
}



### Cluster positions by using a threshold
### Output are two lists: one contains the vectors of the row number
### and the other vectors of the position
# pos = position
# thresh = threshold for clusters

get.cluster <- function(pos, thresh=350, ord=NULL)  {
  i <- 1
  j <- 0
  k <- -1
  n <- length(pos)
  dpos <- pos[2:n] - pos[1:(n-1)]
 
  clustpos <- NULL
  clustrow <- NULL
  while (i < n)  {
    d <- dpos[i]
    if (d <= thresh)  {
      ## if cluster threshold is observed
      if ((i - k) == 1)  {
        ## add to same cluster
        if (is.null(ord))  {
          clustrow[j] <- list(c(clustrow[[j]], i+1))
        }  else  {
          clustrow[j] <- list(c(clustrow[[j]], ord[i+1]))
        }
        clustpos[j] <- list(c(clustpos[[j]], pos[i+1]))
      }  else  {
        ## start new cluster
        j <- j+1
        if (is.null(ord))  {
          clustrow[j] <- list(c(i, i+1))
        }  else  {
          clustrow[j] <- list(c(ord[i], ord[i+1]))
        }
        clustpos[j] <- list(c(pos[i], pos[i+1]))
      }
      k <- i
      i <- i+1
    }
    else  {
      i <- i+1
    }
  }
  list(clustrow=clustrow, clustpos=clustpos)
}


### collate clusters in matrix, by using the begin and end positions of cluster

get.cluster.matrix <- function(pos, lib, thresh=350, ord=NULL)  {
  i <- 1
  j <- 0
  k <- -1
  n <- length(pos)
  dpos <- pos[2:n] - pos[1:(n-1)]
 
  clustpos <- data.frame(0, nrow=length(pos), ncol=4)
  while (i < n)  {
    print(i)
    d <- dpos[i]
    if (d <= thresh)  {
      ## if cluster threshold is observed
      if ((i - k) == 1)  {
        ## add to same cluster, use min and max position
        clustpos[j,1:2] <- c(clustpos[j,1], pos[i+1])
        clustpos[j,3] <- clustpos[j,3]+1
        id <- as.numeric(unlist(strsplit(clustpos[j,4],",")))
        clustpos[j,4] <- paste(sort(unique(c(id, lib[i+1]))), collapse=",")
      }  else  {
        ## start new cluster
        j <- j+1
        clustpos[j,1:2] <- c(pos[i], pos[i+1])
        clustpos[j,3] <- 1
        clustpos[j,4] <- paste(sort(unique(c(lib[i], lib[i+1]))), collapse=",")
      }
      k <- i
      i <- i+1
    }
    else  {
      i <- i+1
    }
  }
  clustpos
}


### Links the clusters for the start positions to the clusters of the
### end positions by using the row numbers
### Output are two lists: one contains the row numbers and the other the size

find.corresp.cluster <- function(row1, row2)  {
  r <- NULL
  size1 <- NULL
  size2 <- NULL
  intsize <- NULL
  for (i in 1:length(row1))  {
    ## Intersection of row1 and row2, using row1 as the reference
    allrows <- lapply(row2, intersect, row1[[i]])
    len <- unlist(lapply(allrows, length))
    ind <- which(len!=0)
    ## Intersecting row numbers
    r[i] <- list(ind)
    ## Size of intersection
    intsize[i] <- list(len[ind])
    ## Size of the corresponding linked cluster
    size1[i] <- list(unlist(lapply(row1[i], length)))
    size2[i] <- list(unlist(lapply(row2[ind], length)))
  }
  list(r=r, intsize=intsize, size1=size1, size2=size2)
}

## returns only the row numbers and the size of intersection
find.corresp.cluster.basic <- function(row1, row2)  {
  r <- NULL
  intsize <- NULL
  for (i in 1:length(row1))  {
    ## Intersection of row1 and row2, using row1 as the reference
    allrows <- lapply(row2, intersect, row1[[i]])
    len <- unlist(lapply(allrows, length))
    ind <- which(len!=0)
    ## Intersecting row numbers
    r[i] <- list(ind)
    ## Size of intersection
    intsize[i] <- list(len[ind])
  }
  list(r=r, intsize=intsize)
}

## returns only the row numbers and the size of intersection
find.corresp.cluster.basic2 <- function(row1, row2)  {
  r <- NULL
  intsize <- NULL
  ## Use only clusters > 1
  csize <- unlist(lapply(row1, length))
  indsize <- which(csize > 1)
  for (i in 1:length(row1[indsize]))  {
    ## Intersection of row1 and row2, using row1 as the reference
    allrows <- lapply(row2, intersect, row1[indsize][[i]])
    len <- unlist(lapply(allrows, length))
    ind <- which(len!=0)
    ## Intersecting row numbers
    r[i] <- list(ind)
    ## Size of intersection
    intsize[i] <- list(len[ind])
  }
  list(r=r, intsize=intsize)
}


### Get the link between all start and end clusters and save them to an RData file

### for one chromosome
get.corresp.clusters <- function(start, ende, fn)  {
  ## get relevant chromosome start and ends
  bg <- sort(start);                         print(any((bg-start)!=0))
  ed <- ende
  ord <- order(ed)
  ed <- ed[ord]
  ## cluster starts and ends
  cbg <- get.cluster(bg) 
  ced <- get.cluster(ed, ord=ord) 
  ## cluster of start and end positions, contains the row numbers
  bgrow <- cbg$clustrow
  edrow <- ced$clustrow
  ## Find the corresponding clusters
  bgcorrclust <- find.corresp.cluster.basic2(bgrow, edrow)
  ## Save corresponding cluster data
  nfn <- paste(fn, ".RData",sep="")
  save(bgcorrclust, file=nfn)
}

### Get the corresponding clusters for start and end positions for all chromosomes
# chr = vector of chromosomes
# start = vector of start positions
# end = vector of end positions
# fn = file name
get.all.corresp.clusters <- function(chr, start, ende, fn)  {
  chrlevels <- c(as.character(c(1:22)),"X","Y")
  chrlen <- length(chrlevels)
  for (i in 1:chrlen)  {
    ## get relevant chromosome start and ends
    ind <- which(chr==chrlevels[i])
    bg <- sort(start[ind]);                         print(paste("i =", i, ",  ", any((bg-start[ind])!=0)))
    ed <- ende[ind]
    ord <- order(ed)
    ed <- ed[ord]
    ## cluster starts and ends
    cbg <- get.cluster(bg) 
    ced <- get.cluster(ed, ord=ord) 
    ## cluster of start and end positions, contains the row numbers
    bgrow <- cbg$clustrow
    edrow <- ced$clustrow
    ## Find the corresponding clusters
    bgcorrclust <- find.corresp.cluster(bgrow, edrow)
    ## Save corresponding cluster data
    nfn <- paste(fn, i, ".RData",sep="")
    save(bgcorrclust, file=nfn)
  }
}

### difference to get.all.corresp.clusters that here it is one chr at a time
get.all.corresp.clusters.chr <- function(chr, bg, ed, fn)  {
  ## get relevant chromosome start and ends
  if (any((sort(bg)-bg)!=0))   print("error")
  ord <- order(ed)
  ed <- ed[ord]
  ## cluster starts and ends
  cbg <- get.cluster(bg) 
  ced <- get.cluster(ed, ord=ord) 
  ## cluster of start and end positions, contains the row numbers
  bgrow <- cbg$clustrow
  edrow <- ced$clustrow
  ## Find the corresponding clusters
  bgcorrclust <- find.corresp.cluster.basic(bgrow, edrow)
  ## Save corresponding cluster data
  save(bgcorrclust, file=fn)
}


### Collate all clusters in the final format 

get.all.clusters.chr <- function(chrnum, bg, ed, bgrow, edrow, bgcorrclust, flag, lib)   {
  ## list of row numbers (e.g. end positions) that are linked to the other clusters (e.g. start positions)
  rbg <- bgcorrclust$r
  ## remove start clusters that are not linked to any end clusters
  indf <- which(unlist(lapply(rbg, function(v) length(v)!=0)))
  ## list of cluster sizes of the intersecting clusters
  intsbg <- bgcorrclust$intsize
  iclust <- length(unlist(intsbg))
  ## iterate through clusters
  clust <- matrix(0, nrow=iclust, ncol=13)
  m <- length(rbg)
  k <- 1
  for (i in 1:m)  {
    ## get cluster indices
    ind <- rbg[[i]]
    if (length(ind)!=0)  {
      ## get row indices for identified clusters
      indrowbg <- bgrow[[i]]
      indrowedlist <- edrow[ind]
      intsbgv <- unlist(intsbg[i])
      ## iterate through clusters that are not 1:1
      for (j in 1:length(ind))  {
        ## get row indices for identified clusters
        indrowed <- indrowedlist[j]
        ## get corresponding rows, i.e. pairs
        introws <- intersect(indrowbg, unlist(indrowed));
        ## get cluster boundaries
        start1 <- min(bg[introws])
        end1 <- max(bg[introws])
        start2 <- min(ed[introws])
        end2 <- max(ed[introws])
        ## get intervals
        intvl1 <- end1 - start1 + 1
        intvl2 <- end2 - start2 + 1
        sv <- start2 - end1 - 1
        ## get cluster size
        n <- intsbgv[j]
        ## get average length
        avlen <- round(mean(ed[introws]-bg[introws]),1)
        ## get type
        types <- factor(flag[introws], levels=c(1,2,4,8,18,130))
        tt <- as.data.frame(table(types))
        indtt <- which(tt[,2]!=0)
        nametype <- tt[indtt,1]
        f <- tt[indtt,2]
        type <- paste(nametype, " (", f, ")", sep="", collapse=" ")
        ## get library code
        ids <- sort(unique(lib[introws]))
        id <- paste(ids, collapse=",")
        ## collate columns
        clust[k,] <- c(i, chrnum, start1, intvl1, end1, sv, start2, intvl2, end2, n, avlen, type, id)
        k <- k+1
      }
    }
  }
  clust <- data.frame(clust)
  names(clust) <- c("clust","chr","start1","intvl1","end1","sv","start2","intvl2","end2","n","avlen","type","id")
  clust
}



### Collate all clusters in the final matrix (!) format, no libraries are given

get.all.clusters.mformat <- function(chrnum, bg, ed, bgrow, edrow, bgcorrclust, flag)   {
  ## list of row numbers (e.g. end positions) that are linked to the other clusters (e.g. start positions)
  rbg <- bgcorrclust$r
  ## remove start clusters that are not linked to any end clusters
  indf <- which(unlist(lapply(rbg, function(v) length(v)!=0)))
  ## list of cluster sizes of the intersecting clusters
  intsbg <- bgcorrclust$intsize
  iclust <- length(unlist(intsbg))
  ## iterate through clusters
  #clust <- matrix(0, nrow=iclust, ncol=17)
  clust <- NULL
  m <- length(rbg)
  k <- 1
  for (i in 1:m)  {
    ## get cluster indices
    ind <- rbg[[i]]
    if (length(ind)!=0)  {
      ## get row indices for identified clusters
      indrowbg <- bgrow[[i]]
      indrowedlist <- edrow[ind]
      intsbgv <- unlist(intsbg[i])
      ## iterate through clusters that are not 1:1
      for (j in 1:length(ind))  {
        ## get row indices for identified clusters
        indrowed <- indrowedlist[j]
        ## get corresponding rows, i.e. pairs
        introws <- intersect(indrowbg, unlist(indrowed));
        if (length(introws)!=0)  {
          ## get cluster boundaries
          start1 <- min(bg[introws])
          end1 <- max(bg[introws])
          start2 <- min(ed[introws])
          end2 <- max(ed[introws])
          ## get intervals
          intvl1 <- end1 - start1 + 1
          intvl2 <- end2 - start2 + 1
          sv <- start2 - end1 - 1
          ## get cluster size
          n <- intsbgv[j]
          ## get average length
          avlen <- round(mean(ed[introws]-bg[introws]),1)
          ## get type
          types <- factor(flag[introws], levels=c(1,2,4,8,18,130))
          tt <- as.data.frame(table(types))
          type <- tt[,2]
          ## collate columns
          clust <- rbind(clust, c(k, chrnum, start1, intvl1, end1, sv, start2, intvl2, end2, n, avlen, type))
          k <- k+1
        }
      }
    }
  }
  clust <- data.frame(clust)
  typname <- paste("t", c(1,2,4,8,18,130), sep="")
  names(clust) <- c("clust","chr","start1","intvl1","end1","sv","start2","intvl2","end2","n","avlen",typname)
  clust
}


### Collate all clusters for all chr in the final matrix (!) format 

get.all.clusters.chr.mformat <- function(chrnum, bg, ed, bgrow, edrow, bgcorrclust, flag, lib)   {
  ## list of row numbers (e.g. end positions) that are linked to the other clusters (e.g. start positions)
  rbg <- bgcorrclust$r
  ## remove start clusters that are not linked to any end clusters
  indf <- which(unlist(lapply(rbg, function(v) length(v)!=0)))
  ## list of cluster sizes of the intersecting clusters
  intsbg <- bgcorrclust$intsize
  iclust <- length(unlist(intsbg))
  ## iterate through clusters
  clust <- matrix(0, nrow=iclust, ncol=33)
  m <- length(rbg)
  k <- 1
  for (i in 1:m)  {
    ## get cluster indices
    ind <- rbg[[i]]
    if (length(ind)!=0)  {
      ## get row indices for identified clusters
      indrowbg <- bgrow[[i]]
      indrowedlist <- edrow[ind]
      intsbgv <- unlist(intsbg[i])
      ## iterate through clusters that are not 1:1
      for (j in 1:length(ind))  {
        ## get row indices for identified clusters
        indrowed <- indrowedlist[j]
        ## get corresponding rows, i.e. pairs
        introws <- intersect(indrowbg, unlist(indrowed));
        ## get cluster boundaries
        start1 <- min(bg[introws])
        end1 <- max(bg[introws])
        start2 <- min(ed[introws])
        end2 <- max(ed[introws])
        ## get intervals
        intvl1 <- end1 - start1 + 1
        intvl2 <- end2 - start2 + 1
        sv <- start2 - end1 - 1
        ## get cluster size
        n <- intsbgv[j]
        ## get average length
        avlen <- round(mean(ed[introws]-bg[introws]),1)
        ## get type
        types <- factor(flag[introws], levels=c(1,2,4,8,18,130))
        tt <- as.data.frame(table(types))
        type <- tt[,2]
        ## get library code
        ids <- factor(lib[introws], levels=c(1:16))
        idt <- as.data.frame(table(ids))
        id <- idt[,2]
        ## collate columns
        clust[k,] <- c(i, chrnum, start1, intvl1, end1, sv, start2, intvl2, end2, n, avlen, type, id)
        k <- k+1
      }
    }
  }
  clust <- data.frame(clust)
  typname <- paste("t", c(1,2,4,8,18,130), sep="")
  libname <- paste("lib", c(1:16), sep="")
  names(clust) <- c("clust","chr","start1","intvl1","end1","sv","start2","intvl2","end2","n","avlen",typname, libname)
  clust
}


### Allocate same cluster number to overlapping clusters
### not used !!!
allocate.cluster.num <- function(clust)  {
  cstart <- as.numeric(as.character(clust$start1))
  cend <- as.numeric(as.character(clust$end2))
  
  indover <- get.overlaps4(cbind(cstart, cend))
  ## index for clusters
  inds <- which(unlist(lapply(indover, length))>1)       ## 173 clusters have overlaps in chr1
  
  overlaps <- NULL
  for (j in 1:length(indover))  {
    ## find intersection for clusters
    if (is.element(j, inds))  {
      overlaps[j] <- list(which(lapply(indover, function(v) which(length(intersect(v,indover[[j]]))!=0))!=0))
      
    }
    ## if single element allocate j
    else  {
      if (indover[[j]]==j)  {
        overlaps[j] <- list(j)
      }  else  {
        print(paste("error in", j))
      }
    }
  }
  k <- length(overlaps)
  cnum <- 0
  i <- 1
  nclust <- NULL
  while (i <= k)  {
    cnum <- cnum + 1
    ind <- overlaps[[i]]
    cn <- length(ind)
    nclust <- rbind(nclust, data.frame(rep(cnum,cn), clust[ind,]))
    i <- max(ind)+1
  }
  names(nclust) <- c("clustnum", names(clust))
  nclust
}


### using connected components
allocate.cluster.num2 <- function(clust)  {
  cstart <- as.numeric(as.character(clust$start1))
  cend <- as.numeric(as.character(clust$end2))
  ## get connected components, store in matrix
  indover <- get.overlaps5(cbind(cstart, cend))
  ## get clusters via connected components
  overlaps <- get.conn.comp(indover)
  ## allocate cluster numbers
  k <- nrow(overlaps)
  cnum <- 0
  i <- 1
  nclust <- NULL
  while (i <= k)  {
    cnum <- cnum + 1
    ind <- which(overlaps[i,]>0)
    cn <- length(ind)
    nclust <- rbind(nclust, data.frame(rep(cnum,cn), clust[ind,]))
    i <- max(ind)+1
  }
  names(nclust) <- c("clustnum", names(clust))
  nclust
}


### Collate all clusters, add the same cluster number for overlapping clusters and merge clusters
# chr = integer that indicates chr
# bg = start positions
# ed = end positions
# pflag = paired flag
# lib = library id
# fn = file name that contains the corresponding clusters

collate.all.clusters.chr <- function(chr, bg, ed, pflag, lib, fn)  {
  fclust <- NULL
  nclust <- NULL
  mclust <- NULL
  ## get relevant chromosome start and ends, pflag and library id
  print(paste("i =", i, ",  ", any((bg-sort(bg))!=0)))
  ord <- order(ed)
  ed <- ed[ord]
  ## cluster starts and ends
  cbg <- get.cluster(bg) 
  ced <- get.cluster(ed, ord=ord) 
  ## cluster of start and end positions, contains the row numbers
  bgrow <- cbg$clustrow
  edrow <- ced$clustrow
  ## get the corresponding clusters
  load(fn)
  ## get clusters
  clust <- get.all.clusters.chr(chr, bg, ed, bgrow, edrow, bgcorrclust, pflag, lib)
  fclust <- rbind(fclust, clust)   
  ## allocate cluster numbers indicating overlapping clusters
  clust2 <- allocate.cluster.num2(clust)
  nclust <- rbind(nclust, clust2)  
  ## merge clusters
  clust3 <- merge.clusters2(clust2[,4:14]) 
  chrnum <- rep(chr, length(clust3))
  fclust3 <- paste(chrnum, clust3, sep="\t")
  mclust <- c(mclust, fclust3)
  fmclust <- vecstr2matstr(mclust, "\t")
  fmclust <- data.frame(fmclust)
  names(fmclust) <- c("chr","start1","intvl1","end1","sv","start2","intvl2","end2","n","avlen","type","id")
  list(fclust=fclust, nclust=nclust, fmclust=fmclust)
}



### Filter for anomalous insert sizes in the Trio-CEU of the 1000
### genomes project depending on the minimum mean and sd of the insert
### sizes for each library combination
# dat = vector containing SV size and library columns
# mu = vector of insert means
# sde = vector of insert SDs

filter.sv <- function(dat, mu, sde)  {
  sv <- dat[1] + 36
  lib <- dat[-1]
  ind <- which(lib!=0)
  max.mu <- max(mu[ind])
  max.sd <- max(sde[ind])
  sv >= max.mu + 4*max.sd
}


### Combine 2 rows of singleton in order to find plus clusters matched
### to minus clusters. Matrix needs to be prepared so that two rows
### are in one line.
# v =  two rows in one line

combine.2.rows.of.singletons <- function(v)  {
  chr <- v[1]
  fst.row <- as.numeric(v[2:21])
  snd.row <- as.numeric(v[23:42])
  start1 <- fst.row[1]
  intvl1 <- fst.row[2] + 35
  end1  <- fst.row[3] + 35
  sv <- snd.row[1] - end1 + 1
  start2 <- snd.row[1]
  intvl2 <- snd.row[2] + 35
  end2 <- snd.row[3] + 35
  n1 <- fst.row[4]
  n2 <- snd.row[4]
  avqual <- round((fst.row[4]*fst.row[5] + snd.row[4]*snd.row[5]) / (fst.row[4]+snd.row[4]),1)
  types <- fst.row[6:13] + snd.row[6:13]
  libs <- fst.row[14:19] + snd.row[14:19]
  orient <- c(fst.row[20], snd.row[20])
  n.row <- c(chr,start1,intvl1,end1,sv,start2,intvl2,end2,n1,n2,avqual,types,libs,orient)
  n.row
}

combine.2.rows.of.singletons.basic <- function(v)  {
  chr <- v[1]
  fst.row <- as.numeric(v[2:15])
  snd.row <- as.numeric(v[17:30])
  start1 <- fst.row[1]
  intvl1 <- fst.row[2] + 35
  end1  <- fst.row[3] + 35
  sv <- snd.row[1] - end1 + 1
  start2 <- snd.row[1]
  intvl2 <- snd.row[2] + 35
  end2 <- snd.row[3] + 35
  n1 <- fst.row[4]
  n2 <- snd.row[4]
  avqual <- round((fst.row[4]*fst.row[5] + snd.row[4]*snd.row[5]) / (fst.row[4]+snd.row[4]),1)
  types <- fst.row[6:13] + snd.row[6:13]
  orient <- c(fst.row[14], snd.row[14])
  n.row <- unlist(c(chr,start1,intvl1,end1,sv,start2,intvl2,end2,n1,n2,avqual,types,orient))
  n.row
}

### Reformat singletons so that plus and minus clusters can be combined
# datp = clusters of singletons on the positive strand
# datm = clusters of singletons on the negative strand

reformat.singletons <- function(datp, datm)  {
  n <- ncol(datp)
  id <- rep(1, nrow(datp))
  datp <- data.frame(datp, id)
  id <- rep(2, nrow(datm))
  datm <- data.frame(datm, id)
  dat <- rbind(datp, datm)
  ## replace X by 23 and Y by 24
  chr <- as.character(dat[,1])
  indx <- which(chr=="X")
  indy <- which(chr=="Y")
  if (length(indx)!=0)  {
    chr[indx] <- "23"
  }
  if (length(indy)!=0)  {
    chr[indy] <- "24"
  }
  dat[,1] <- chr
  ord <- order(dat[,1],dat[,2])
  dat <- dat[ord,]
  ## loop through chr
  f.dat <- NULL
  allchr <- levels(as.factor(dat[,1]))
  for (j in allchr)  {
    indchr <- which(dat[,1]==j); print(paste(j, length(indchr)))
    n <- length(indchr)
    if (n > 1)  {
      comb.dat <- data.frame(dat[indchr,][1:(n-1),], dat[indchr,][2:n,])
      chr.dat <- t(apply(comb.dat, 1, combine.2.rows.of.singletons.basic))       
      f.dat <- rbind(f.dat, chr.dat)
    } # else  {
      #f.dat <- rbind(f.dat, dat[indchr,])
    #}
  }
  f.dat <- data.frame(f.dat)
  f.dat
}
  

##### Edit annotation #####

### convert numeric chr into string chr

chr2string <- function(dat)  {
  chr <- paste("chr",as.character(dat[,1]), sep="")
  chr[chr=="chr23"] <- "chrX"
  chr[chr=="chr24"] <- "chrY"
  dat <- data.frame(chr, dat[,-1], row.names=1:nrow(dat))
  dat
}
#### Converts string chr into numeric chr and returns set ordered by chr and start position
## mat = matrix with chr in 1st column

edit.chr <- function(mat, chrn=1, posn=2, ord=TRUE)  {
  chr <- as.character(mat[,chrn])
  chr <- gsub("chr","", chr)
  chr <- gsub("Chr","", chr)
  chr[chr=="X"] <- 23
  chr[chr=="Y"] <- 24
  mat[,chrn] <- as.numeric(chr)
  if(ord)  {
    ## order by chr and start position
    ord <- order(as.numeric(as.character(mat[,chrn])), as.numeric(as.character(mat[,posn])))
    mat <- data.frame(mat[ord,], stringsAsFactors=F)
  }
  mat
}

### Edit format of clusters in libraries in 1000 genomes project
# dat = Trio-CEU data frame

edit.format.lib <- function(dat)  {
  ndat <- dat[,c(1,4,6,9:11,20:25)]
  names(ndat) <- c("chr","start","end","n","avlen","avqual", "NA12878_1", "NA12878_2", "NA12891_1", "NA12891_2", "NA12892_1", "NA12892_2")
  ndat
}

### Edit format of clusters in libraries in 1000 genomes project
# dat = Trio-YRI data frame

edit.format.lib.yri <- function(dat, libname)  {
  ndat <- dat[,c(1,4,6,9:11,20:ncol(dat))]
  names(ndat) <- c("chr","start","end","n","avlen","avqual", libname)
  ndat
}

### Edit format of clusters in individuals in 1000 genomes project
# dat = Trio-CEU data frame

edit.format.individ <- function(dat)  {
  ndat <- dat[,c(1,4,6,9:11)]
  ## combine libraries for individuals
  ndat <- data.frame(ndat, dat[,20]+dat[,21], dat[,22]+dat[,23], dat[,24]+dat[,25])
  names(ndat) <- c("chr","start","end","n","avlen","avqual", "NA12878", "NA12891", "NA12892")
  ndat
}


### Edit format of clusters in the Trio-CEU in 1000 genomes project
# dat = Trio-CEU data frame

edit.format.trio <- function(dat)  {
  ndat <- dat[,c(1,4,6,9:11)]
  ## combine libraries for individuals
  ndat <- data.frame(ndat, dat[,20]+dat[,21], dat[,22]+dat[,23], dat[,24]+dat[,25])
  names(ndat) <- c("chr","start","end","n","avlen","avqual", "NA12878", "NA12891", "NA12892")
  ndat
}

### Edit format of clusters in the Trio-YRI in 1000 genomes project
# dat = Trio-YRI data frame
# cols1 = columns of 1st individual NA19238
# cols2 = columns of 1st individual NA19239
# cols3 = columns of 1st individual NA19240
# libnames = c("NA19238","NA19239","NA19240")

edit.format.trio.yri <- function(dat, ncols1, ncols2, ncols3, libnames)  {
  ndat <- dat[,c(1,4,6,9:11)]
  ## combine libraries for individuals
  ndat <- data.frame(ndat, apply(dat[,ncols1],1,sum), apply(dat[,ncols2],1,sum), apply(dat[,ncols3],1,sum))
  names(ndat) <- c("chr","start","end","n","avlen","avqual", libnames)
  ndat
}


### Get annotation for overlaps
# mat1 = reference matrix
# mat2 = matrix searched for overlaps

get.ann.chr <-  function(mat1, mat2)  {
  ind <- NULL
  ## loop through chr
  for (i in 1:24)  {
    print(i)
    indd <- which(as.numeric(as.character(mat1[,1]))==i)
    indann <- which(as.numeric(as.character(mat2[,1]))==i)
    ind[i] <- list(apply(mat1[indd,2:3], 1, get.overlaps3, mat2[indann,2:3]))
  }
  ind
}

get.ann.chr.thresh <-  function(mat1, mat2, thresh=0.5)  {
  ind <- NULL
  ## loop through chr
  for (i in 1:24)  {
    indd <- which(as.numeric(as.character(mat1[,1]))==i)
    indann <- which(as.numeric(as.character(mat2[,1]))==i)
    ind[i] <- list(apply(mat1[indd,2:3], 1, get.reciprocal.overlaps, mat2[indann,2:3], thresh=thresh))
  }
  ind
}

### get reciprocal overlap - need to make sure that row.names=c(1:nrow(mat))
# mat1 = matrix overlaps are searched for
# mat2 = matrix overlaps are searched against, so mat1 in mat2
# thresh = 0.5 (at least 50% overlap)
get.reciprocal.overlaps.chr <-  function(mat1, mat2, thresh=0.5)  {
  ind <- NULL
  ## loop through chr
  for (i in 1:24)  {
    #print(i)
    indd <- which(as.numeric(as.character(mat1[,1]))==i)
    indann <- which(as.numeric(as.character(mat2[,1]))==i)
    if (length(indd)!=0 & length(indann)!=0)  {     ## added on 24-11-2011 in case not all chroms are used
      over <- apply(mat1[indd,2:3], 1, get.reciprocal.overlaps, mat2[indann,2:3], thresh=thresh)
      if (is.matrix(over))  {           ## edited on 30/6/2010, because if there is just one element that has a hit against more than one,
        over <- over[1,1]               ## the returned class is a matrix and it does not have a name which causes problems below
      }
      ## if (is.list(over))  {             ## edited on 28/6/2012, if it is already a list, I don't want a list of a list
##         ind[i] <- list(unlist(over))
##       }
      #else  {
        ind[i] <- list(over)
      #}
    }
  }
  v <- rep(0, nrow(mat1))
  indname <- as.numeric(names(which(unlist(sapply(ind, sapply, length))>0)))
  v[indname] <- 1
  v
}

### get any bp overlap - need to make sure that row.names=c(1:nrow(mat))
# mat1 = matrix overlaps are searched for
# mat2 = matrix overlaps are searched against, so mat1 in mat2
get.any.bp.overlap.chr <-  function(mat1, mat2)  {
  ind <- NULL
  ## loop through chr
  for (i in 1:24)  {
    print(i)
    indd <- which(as.numeric(as.character(mat1[,1]))==i)
    indann <- which(as.numeric(as.character(mat2[,1]))==i)
    if (length(indd)!=0 & length(indann)!=0)  {
      over <- apply(mat1[indd,2:3], 1, get.overlaps3, mat2[indann,2:3])
      if (is.matrix(over))  {           ## edited on 30/6/2010, because if there is just one element that has a hit against more than one,
        over <- over[1,1]               ## the returned class is a matrix and it does not have a name which causes problems below
      }
      ind[i] <- list(over)
    }
  }
  v <- rep(0, nrow(mat1))
  indname <- as.numeric(names(which(unlist(sapply(ind, sapply, length))>0)))
  v[indname] <- 1
  v
}

### get reciprocal overlap - need to make sure that row.names=c(1:nrow(mat))
# mat1 = matrix overlaps are searched for
# mat2 = matrix overlaps are searched against, so mat1 in mat2
# thresh = 0.5 (at least 50% overlap)
get.reciprocal.overlaps.chr.p <-  function(mat1, mat2, thresh=0.5)  {
  ind <- NULL
  ## loop through chr
  for (i in 1:24)  {
    #print(i)
    indd <- which(as.numeric(as.character(mat1[,1]))==i)
    indann <- which(as.numeric(as.character(mat2[,1]))==i)
    listover <- apply(mat1[indd,2:3], 1, get.reciprocal.overlaps.p, mat2[indann,2:3], thresh=thresh)
    over <- listover$ind
    if (is.matrix(over))  {      
      over <- over[1,1]          
    }
    ind[i] <- list(over)
  }
  v <- rep(0, nrow(mat1))
  p <- rep(0, nrow(mat1))
  indname <- as.numeric(names(which(unlist(sapply(ind, sapply, length))>0)))
  v[indname] <- 1
  p[indname] <- 1   ## not finished
  v
}

### get reciprocal overlap - need to make sure that row.names=c(1:nrow(mat))
# mat1 = matrix overlaps are searched for
# mat2 = matrix overlaps are searched against, so mat1 in mat2
# thresh = 0.5 (at least 50% overlap)
get.reciprocal.overlaps.chr.old <-  function(mat1, mat2, thresh=0.5, d=500)  {
  ind <- NULL
  v <- rep(0, nrow(mat1))
  ## loop through chr
  for (i in 1:24)  {
    #print(i)
    indd <- which(as.numeric(as.character(mat1[,1]))==i)
    indann <- which(as.numeric(as.character(mat2[,1]))==i)
    over <- apply(mat1[indd,2:3], 1, get.reciprocal.overlaps, mat2[indann,2:3], thresh=thresh, d)
    if (is.matrix(over))  {                       ## edited on 30/6/2010, because if there is just one element that
      v[indd][as.numeric(colnames(over))] <- 1    ## has a hit against more than one, the returned class is a matrix
    }  else if (is.list(over))  {                 ## and it does not have a name which causes problems below
      indname <- as.numeric(names(which(unlist(sapply(over, length))>0)))
      v[indname] <- 1                             ## edited on 10/3/2011 because vector (=1:1 relationship) was not accounted for
    }  else  {
      v[over==1] <- 1 
    }
  }
  v
}


### Compile annotation
# indf = list of list of indices for all annotations (1st level) and all chr (2nd level)
# which.ann = index of annotation required
# ann = annotation data vector

get.annotation <- function(dels, indf, which.ann, ann)  {
  s.ann <- vector(length=nrow(dels))
  ## specify annotation
  ind <- indf[[which.ann]]
  ## loop through chr
  for (i in 1:24)  {
    indchr <- ind[[i]]
    indav <- which(sapply(indchr, length)!=0)
    rn <- as.numeric(names(indchr[indav]))
    if (length(indav)!=0)  {
      for (j in 1:length(indav))  {
        s.ann.ind <- indchr[[indav[j]]]
        s.ann[rn[j]] <- paste(ann[s.ann.ind], collapse=",")
      }
    }
  }
  ind <- which(s.ann=="FALSE")
  s.ann[ind] <- NA
  s.ann
}

### Function to add the overlapping deletions from other sources, returns vector
# dat = RP data
# dat2 = overlapping data from other sources

get.ann.overlap <- function(dat, dat2, colnum)  {
  ann <- NULL
  k <- 1
  test <- NULL
  ## loop through chr
  for (i in 1:24)  {
    #print(i)
    ind1 <- which(dat[,1]==i); #test <- sum(test,length(ind1)); print(test)
    ind2 <- which(dat2[,1]==i)
    ord <- order(dat2[ind2,2])
    dat2.ord <- dat2[ind2,][ord,]
    if (length(ind1)!=0 & length(ind2)!=0)  {
      ## get 1bp overlaps between dat and dat2
      over <- apply(dat[ind1,c(2,3)], 1, get.overlaps3, mat=dat2.ord[,2:3])   ## start and end position in col 2 and 3
      n <- length(over)
      m <- sapply(over, length)
      ## 1:1 relationship between dat and dat2
      if (!is.matrix(over) & all(m==1) & n!=0)  {
        #print("vector")
        for (j in 1:length(over))  {       ## vector if not list !!!
          ann[k] <- paste(unique(dat2.ord[over[j],colnum]), collapse=",")
          k <- k+1
        }
      }
      ## 1:n relationship between dat and dat2 with equal n
      else if (is.matrix(over) & n!=0)  {
        #print("matrix")
        for (j in 1:ncol(over))  {       ## matrix if not list !!!
          ann[k] <- paste(unique(dat2.ord[over[,j],colnum]), collapse=",")
          k <- k+1
        }
      }
      ## 1:n relationship between dat and dat2 with unequal n
      else if (any(m!=1) & n!=0)  {    ## list !!!
        #print("list")
        ## loop through rows in dat for a given chr
        for (j in 1:n)  {
          if (length(over[[j]])!=0)  {
            ann[k] <- paste(unique(dat2.ord[over[[j]],colnum]), collapse=",")
          }  else  {
            #ann[k] <- NA
            is.na(ann) <- k
          }
          k <- k+1
        }
      }
      ## no overlap
      else  {
        ny <- length(ind1)
        #ann <- c(ann, rep(NA, ny))
        na.vec <- NULL
        is.na(na.vec) <- 1:ny
        ann <- c(ann, na.vec)
        k <- k + ny
      }
    }
    ## no data for that particular chr
    else if (length(ind1)!=0 & length(ind2)==0)  { 
      ny <- length(ind1)
      #ann <- c(ann, rep(NA, ny))
      na.vec <- NULL
      is.na(na.vec) <- 1:ny
      ann <- c(ann, na.vec)
      k <- k + ny
    }
    #print(paste("ann =", length(ann)))
  }
  ann
}



##### Plot all reads within 1kb of a deletion

# del.bg = deletion start
# del.ed = deletion end
# position = vector of 1st position of read
# insert = vector of insert size
# type = vector of paired flag
# sing.qual = vector of single mapping quality 
# qual = vector of mapping quality 

plot.reads2 <- function(del.bg, del.ed, position, insert, type, sing.qual, qual, thresh=250, del=TRUE, mt=NULL, chr=NULL)   {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  mycols <- c("darkgrey","darkred","red","orange","yellow2","green","blue","darkgreen","black")
  ## generate plot 
  plot.bg <- (del.bg - 1000)/s
  plot.end <- (del.ed + 1000)/s
  intvl <- plot.end - plot.bg
  m <- max(c(qual, sing.qual))
  if (is.null(mt))  {
    mt <- paste("chr",chr,":",del.bg,"-",del.ed," ", sep="")
  }
  plot(c(plot.bg - intvl/20, plot.end + intvl/3), c(-2,m), ty="n", xlab="position", ylab="quality scores", 
       cex.axis=1.3, cex.lab=1.3, cex.main=1.4, main=mt, las=1)
  x <- plot.end + intvl/10
  legend(x, m, legend=c("normal","FF","FR","RF","RR","diff chr","single","S-W","deletion"),
         lty=1, lwd=2, cex=0.8, col=mycols)
  ## find types
  ind1 <- which(type==1)
  ind2 <- which(type==2 | (type==18 & insert > thresh))
  ind4 <- which(type==4)
  ind8 <- which(type==8)
  ind18 <- which(type==18 & insert <= thresh)
  ind32 <- which(type==32)      # different chr
  ind64 <- which(type==64)      # one read not mapped
  ind130 <- which(type==130)    # S-W alignment
  ## plot read pairs for all types
  indlist <- list(ind18, ind1, ind2, ind4, ind8, ind32, ind64, ind130)
  for (j in 1:length(indlist))  {
    ind <- indlist[[j]]
    bg <- as.integer(position[ind])
    ed <- bg + as.integer(insert[ind])
    qu <- qual[ind]
    squ <- sing.qual[ind]
    n <- length(bg)
    if (n != 0)  {
      for (k in 1:n)  {
        if (j==1 | j==2 | j==3 | j==4 | j==5 | j==8)  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((ed[k]-35)/s, ed[k]/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((bg[k]+35)/s, (ed[k]-35)/s), c(qu[k],qu[k]), lwd=0.5, col=mycols[j])
        }  else  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(squ[k],squ[k]), lwd=2, col=mycols[j])
        }
      }
    }
  }
  if (del==TRUE) {
    ## deletion
    lines(c(del.bg/s, del.ed/s), c(-2,-2), col="black", lwd=2)
  }
}


###### Function for filtering the directory structure for plotting reads around deletions ######

# fd = path of reads
# id = id of individual

filter.dir <- function(fd, id)  {
  d <- dir(fd, full.names=TRUE)
  f <- dir(fd)
  ind <- grep(id, d)
  d <- d[ind]
  f <- f[ind]
  ind <- grep("cluster", d)
  d <- d[ind]
  f <- f[ind]
  ind <- grep("_1.txt", d)
  d1 <- d[ind]
  f <- f[ind]
  ind <- grep("_2.txt", d)
  d2 <- d[ind]
  fnum <- gsub("cluster","", f)
  idname <- paste("_", id, "_1.txt", sep="")
  fnum <- as.integer(gsub(idname, "", fnum))
  ord <- order(fnum)
  fnum <- fnum[ord]
  d1 <- d1[ord]
  d2 <- d2[ord]
  list(d1=d1, d2=d2)
}

# fd = path of reads
# id = id of individual

filter.rd.rp.dir <- function(fd, id, rp=TRUE)  {
  d <- dir(fd, full.names=TRUE)
  f <- dir(fd)
  ind <- grep(id, d)
  d <- d[ind]
  f <- f[ind]
  ind <- grep("cluster", d)
  d <- d[ind]
  f <- f[ind]
  if (is.null(rp))  {
    indrp <- grep(".txt", d)
    d <- d[indrp]
    f <- f[indrp]
    ind <- grep("_1.txt", d)
    d1 <- d[ind]
    f <- f[ind]
    ind <- grep("_2.txt", d)
    d2 <- d[ind]
  }  else  {
    if (rp==TRUE)  {
      indrp <- grep("rp.txt", d)
      d <- d[indrp]
      f <- f[indrp]
    }  else  {
      indrp <- grep("rd.txt", d)
      d <- d[indrp]
      f <- f[indrp]
    }
    ind <- grep("_1_", d)
    d1 <- d[ind]
    f <- f[ind]
    ind <- grep("_2_", d)
    d2 <- d[ind]
  }
  fnum <- vecstr2matstr(f, "_")[,1]
  fnum <- gsub("cluster","", fnum)
  ord <- order(as.numeric(fnum))
  fnum <- fnum[ord]
  d1 <- d1[ord]
  d2 <- d2[ord]
  list(d1=d1, d2=d2)
}


##### Prepares reads from directories for plotting by libraries

# fn = file name of plot
# d = directory structure, taking into account 1st and 2nd lib, output from filter.dir
# dels = start and end position of deletions
# thresh = threshold for anomalous pairs (vector of size 2 for 1st and 2nd lib)
# chr = chr vector

get.plot.reads.lib <- function(fn, d, dels, thresh, chr="20")  {
  pdf(fn, height=3, width=4, pointsize=7, onefile=T)
  #for (i in 1:(length(d[[1]])-1))  {
  for (i in 1:length(d[[1]]))  {
    print(i)
    ## 1st lib
    if (file.info(d[[1]][i])$size > 0)  {
      clust1 <- read.delim(d[[1]][i], header=FALSE)
      lib <- rep(1, nrow(clust1))
      clust1 <- data.frame(clust1, lib)
    }  else  {
      clust1 <- c(V1=0,V2=0,V3="+",V4=0,V5=0,V6=0,V7=0,V8=0,lib=0)
    }
    ## 2nd lib
    if (file.info(d[[2]][i])$size > 0)  {
      clust2 <- read.delim(d[[2]][i], header=FALSE)
      lib <- rep(2, nrow(clust2))
      clust2 <- data.frame(clust2, lib)
    }  else  {
      clust2 <- c(V1=0,V2=0,V3="-",V4=0,V5=0,V6=0,V7=0,V8=0,lib=0)
    }
    ## combine and plot
    clust <- rbind(clust1, clust2)
    #clust <- clust1  # remove again
    plot.reads.lib(dels[i,1],dels[i,2],clust[,2],clust[,4],clust[,5],clust[,7],clust[,8],clust[,9],thresh,chr[i])
  }
  dev.off()
}


##### Prepares reads from directories for plotting in PNG format by libraries

# fn = file name of plot
# d = directory structure, taking into account 1st and 2nd lib, output from filter.dir
# dels = start and end position of deletions
# thresh = threshold for anomalous pairs (vector of size 2 for 1st and 2nd lib)

get.plot.reads.lib.png <- function(fn, d, dels, thresh)  {
  for (i in 1:length(d[[1]]))  {
    png(fn[i])
    print(i)
    ## 1st lib
    if (file.info(d[[1]][i])$size > 0)  {
      clust1 <- read.delim(d[[1]][i], header=FALSE)
      lib <- rep(1, nrow(clust1))
      clust1 <- data.frame(clust1, lib)
    }  else  {
      clust1 <- c(V1=0,V2=0,V3="+",V4=0,V5=0,V6=0,V7=0,V8=0,lib=0)
    }
    ## 2nd lib
    if (file.info(d[[2]][i])$size > 0)  {
      clust2 <- read.delim(d[[2]][i], header=FALSE)
      lib <- rep(2, nrow(clust2))
      clust2 <- data.frame(clust2, lib)
    }  else  {
      clust2 <- c(V1=0,V2=0,V3="-",V4=0,V5=0,V6=0,V7=0,V8=0,lib=0)
    }
    ## combine and plot
    clust <- rbind(clust1, clust2)
    plot.reads.lib(dels[i,1],dels[i,2],clust[,2],clust[,4],clust[,5],clust[,7],clust[,8],clust[,9],thresh)
    dev.off()
  }
}


##### Plot all reads within 1kb of a deletion, library specific

# del.bg = deletion start
# del.ed = deletion end
# position = vector of 1st position of read
# insert = vector of insert size
# type = vector of paired flag
# sing.qual = vector of single mapping quality 
# qual = vector of mapping quality 

plot.reads.lib <- function(del.bg, del.ed, position, insert, type, sing.qual, qual, lib, thresh, chr="20")   {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  mycols <- c("darkgrey","darkred","orange","yellow2","green","blue","darkgreen","red")
  ## generate plot 
  plot.bg <- (del.bg - 1000)/s
  plot.end <- (del.ed + 1000)/s
  intvl <- plot.end - plot.bg
  m <- max(c(qual, sing.qual))
  plot(c(plot.bg - intvl/20, plot.end + intvl/3), c(-2,m), ty="n", xlab="position", ylab="quality scores", 
       cex.axis=1.3, cex.lab=1.3, cex.main=1.4, main=paste("chr",chr, ":",del.bg,"-",del.ed," ", sep=""), las=1)
  x <- plot.end + intvl/10
  legend(x, m, legend=c("normal","FF","FR","RF","RR","diff chr","single","S-W","deletion"),
         lty=1, lwd=2, cex=0.8, col=c("darkgrey","darkred","red","orange","yellow2","green","blue","darkgreen","black"))
  ## find types
  ind1 <- which(type==1)
  ind2 <- c(which(lib==1 & (type==2 | (type==18 & insert > thresh[1]))),   # lib specific thresholds
            which(lib==2 & (type==2 | (type==18 & insert > thresh[2]))))
  ind4 <- which(type==4)
  ind8 <- which(type==8)
  ind18 <- c(which(lib==1 & type==18 & insert <= thresh[1]),               # lib specific thresholds
             which(lib==2 & type==18 & insert <= thresh[2]))
  ind32 <- which(type==32)      # different chr
  ind64 <- which(type==64)      # one read not mapped
  ind130 <- which(type==130)    # S-W alignment
  ## plot read pairs for all types
  indlist <- list(ind18, ind1, ind4, ind8, ind32, ind64, ind130, ind2)
  for (j in 1:length(indlist))  {
    ind <- indlist[[j]]
    bg <- as.integer(position[ind])
    ed <- bg + as.integer(insert[ind])
    qu <- qual[ind]
    squ <- sing.qual[ind]
    n <- length(bg)
    if (n != 0)  {
      for (k in 1:n)  {
        if (j==1 | j==2 | j==3 | j==4 | j==5 | j==8)  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((ed[k]-35)/s, ed[k]/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((bg[k]+35)/s, (ed[k]-35)/s), c(qu[k],qu[k]), lwd=0.1, col=mycols[j])
        }  else  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(squ[k],squ[k]), lwd=2, col=mycols[j])
        }
      }
    }
  }
  ## deletion
  lines(c(del.bg/s, del.ed/s), c(-2,-2), col="black", lwd=2)
}




##### Plot all reads within 1kb of a deletion

# del.bg = deletion start
# del.ed = deletion end
# position = vector of start positions of reads
# insert = vector of insert sizes
# type = vector of paired flags
# sing.qual = vector of single mapping qualities 
# qual = vector of mapping qualities
# thresh = cut-off for anomalous read pairs
# chr = chromosome

plot.reads <- function(del.bg,del.ed,position,insert,type,sing.qual,qual,thresh,chr,primer=NULL,leg.cex=0.8,h=NULL,marg=1000,v=FALSE,left=TRUE)   {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  mycols <- c("darkgrey","darkred","orange","yellow2","green","blue","darkgreen","red")
  ## generate plot 
  plot.bg <- (del.bg - marg)/s
  plot.end <- (del.ed + marg)/s
  intvl <- plot.end - plot.bg
  #m <- max(c(qual, sing.qual), na.rm=TRUE)
  m <- 100
  header <- paste("chr",chr,":",del.bg,"-",del.ed," ", sep="")
  if (!is.null(h))
    header <- paste(h, header)
  plot(c(plot.bg - intvl/20, plot.end + intvl/3), c(-2,m), ty="n", xlab="position", ylab="quality scores", 
       cex.axis=1.3, cex.lab=1.3, cex.main=1.4, main=header, las=1, ylim=c(0,100))
  x <- plot.end + intvl/10
  if (is.null(primer))  {
    legend(x, m, legend=c("normal","FF","FR","RF","RR","diff chr","single","S-W","deletion"),
           lty=1, lwd=2, cex=leg.cex, col=c("darkgrey","darkred","red","orange","yellow2","green","blue","darkgreen","black"))
  }  else  {
    legend(x, m, legend=c("normal","FF","FR","RF","RR","diff chr","single","S-W","deletion","primer"),
           lty=c(rep(1,9),4), lwd=2, cex=leg.cex,
           col=c("darkgrey","darkred","red","orange","yellow2","green","blue","darkgreen","black","goldenrod3"))
  }
  ## find types
  ind1 <- which(type==1)
  ind2 <- which(type==2 | (type==18 & abs(insert) > thresh))  # threshold for anomalous pairs
  ind4 <- which(type==4)
  ind8 <- which(type==8)
  ind18 <- which(type==18 & insert <= thresh)            # threshold for normal pairs
  ind32 <- which(type==32)      # different chr
  ind64 <- which(type==64)      # one read not mapped
  ind130 <- which(type==130)    # S-W alignment
  ## plot read pairs for all types
  indlist <- list(ind18, ind1, ind4, ind8, ind32, ind64, ind130, ind2)
  for (j in 1:length(indlist))  {
    ind <- indlist[[j]]
    bg <- as.integer(position[ind])
    ## for negative inserts bg = bg+35
    bg[insert[ind]<0] <- bg[insert[ind]<0]+35
    ed <- bg + as.integer(insert[ind])
    qu <- as.numeric(as.character(qual[ind]))
    squ <- as.numeric(as.character(sing.qual[ind]))
    n <- length(bg)
    if (n != 0)  {
      for (k in 1:n)  {
        if (j==1 | j==2 | j==3 | j==4 | j==5 | j==8)  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((ed[k]-35)/s, ed[k]/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((bg[k]+35)/s, (ed[k]-35)/s), c(qu[k],qu[k]), lwd=0.1, col=mycols[j])
        }  else  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(squ[k],squ[k]), lwd=2, col=mycols[j])
        }
      }
    }
  }
  ## deletion
  lines(c(del.bg/s, del.ed/s), c(-2,-2), col="black", lwd=2)
  if (v==TRUE)  {
    if (left==TRUE)  {
      abline(v=del.bg/s, col="black", lty=4)
    } else {
     abline(v=del.ed/s, col="black", lty=4)
    } 
  } 
  if (!is.null(primer))  {
    abline(v=primer[1], col="goldenrod3", lty=4, lwd=1.2)
    abline(v=primer[2], col="goldenrod3", lty=4, lwd=1.2)
  }    
}

plot.reads.jitter <- function(del.bg,del.ed,position,insert,type,sing.qual,qual,thresh,chr,primer=NULL,leg.cex=0.8,h=NULL,marg=1000,v=FALSE,left=TRUE)   {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  mycols <- c("darkgrey","darkred","orange","yellow2","green","blue","darkgreen","red")
  ## generate plot 
  plot.bg <- (del.bg - marg)/s
  plot.end <- (del.ed + marg)/s
  intvl <- plot.end - plot.bg
  #m <- max(c(qual, sing.qual), na.rm=TRUE)
  #m <- 100
  m <- 65
  header <- paste("chr",chr,":",del.bg,"-",del.ed," ", sep="")
  if (!is.null(h))
    header <- paste(h, header)
  ## plot deletion size
  plot(c(plot.bg - intvl/20, plot.end + intvl/3), c(-2,m), ty="n", xlab="position", ylab="quality scores", 
       cex.axis=1.3, cex.lab=1.3, cex.main=1.4, main=header, las=1, ylim=c(0,100))
  x <- plot.end + intvl/10
  if (is.null(primer))  {
    legend(x, m, legend=c("normal","FF","FR","RF","RR","diff chr","single","S-W","deletion"),
           lty=1, lwd=2, cex=leg.cex, col=c("darkgrey","darkred","red","orange","yellow2","green","blue","darkgreen","black"))
  }  else  {
    legend(x, m, legend=c("normal","FF","FR","RF","RR","diff chr","single","S-W","deletion","primer"),
           lty=c(rep(1,9),4), lwd=2, cex=leg.cex,
           col=c("darkgrey","darkred","red","orange","yellow2","green","blue","darkgreen","black","goldenrod3"))
  }
  ## find types
  ind1 <- which(type==1)
  ind2 <- which(type==2 | (type==18 & abs(insert) > thresh))  # threshold for anomalous pairs
  ind4 <- which(type==4)
  ind8 <- which(type==8)
  ind18 <- which(type==18 & insert <= thresh)            # threshold for normal pairs
  ind32 <- which(type==32)      # different chr
  ind64 <- which(type==64)      # one read not mapped
  ind130 <- which(type==130)    # S-W alignment
  ## plot read pairs for all types
  indlist <- list(ind18, ind1, ind4, ind8, ind32, ind64, ind130, ind2)
  for (j in 1:length(indlist))  {
    ind <- indlist[[j]]
    bg <- as.integer(position[ind])
    ## for negative inserts bg = bg+35
    bg[insert[ind]<0] <- bg[insert[ind]<0]+35
    ed <- bg + as.integer(insert[ind])
    qu <- as.numeric(as.character(qual[ind]))
    squ <- as.numeric(as.character(sing.qual[ind]))
    n <- length(bg)
    if (n != 0)  {
      for (k in 1:n)  {
        if (j==1 | j==2 | j==3 | j==4 | j==5 | j==8)  {
          y <- jitter(qu[k], 2)
          lines(c(bg[k]/s, (bg[k]+35)/s), c(y,y), lwd=2, col=mycols[j])
          lines(c((ed[k]-35)/s, ed[k]/s), c(y,y), lwd=2, col=mycols[j])
          lines(c((bg[k]+35)/s, (ed[k]-35)/s), c(y,y), lwd=0.1, col=mycols[j])
        }  else  {
          y <- jitter(squ[k], 1)
          lines(c(bg[k]/s, (bg[k]+35)/s), c(y,y), lwd=2, col=mycols[j])
        }
      }
    }
  }
  ## deletion
  lines(c(del.bg/s, del.ed/s), c(-2,-2), col="black", lwd=2)
  if (v==TRUE)  {
    if (left==TRUE)  {
      abline(v=del.bg/s, col="black", lty=4)
    } else {
     abline(v=del.ed/s, col="black", lty=4)
    } 
  } 
  if (!is.null(primer))  {
    abline(v=primer[1], col="goldenrod3", lty=4, lwd=1.2)
    abline(v=primer[2], col="goldenrod3", lty=4, lwd=1.2)
  }    
}

plot.reads.paper <- function(del.bg,del.ed,position,insert,type,sing.qual,qual,thresh,chr,leg.cex=0.8,h=NULL,marg=1000)   {
  par(mai=c(0.6,0.6,0.3,0.1))
  #par(mai=c(1.2,1.2,1,0.5))     ## for png
  s <- 1000000    ## scaling factor
  mycols <- c("darkgrey","red")
  ## generate plot 
  plot.bg <- (del.bg - marg)/s
  plot.end <- (del.ed + marg)/s
  intvl <- plot.end - plot.bg
  m <- 100
  header <- paste("chr",chr,":",del.bg,"-",del.ed," ", sep="")
  if (!is.null(h))
    header <- paste(h, header, sep="")
  plot(c(plot.bg - intvl/20, plot.end + intvl/3), c(-2,m), ty="n", xlab="position", ylab="quality scores", 
       cex.axis=1.3, cex.lab=1.3, cex.main=1.4, main=header, las=1, ylim=c(0,100))
  x <- plot.end + intvl/50
  legend(x, m, legend=c("normal RP","extend RP","deletion"),
         lty=1, lwd=2, cex=leg.cex, col=c("darkgrey","red","black"))
  ## find deletion RPs
  ind2 <- which(abs(insert) > thresh)  # threshold for anomalous pairs
  ind18 <- which(insert <= thresh)            # threshold for normal pairs
  ## plot read pairs for all types
  indlist <- list(ind18, ind2)
  for (j in 1:length(indlist))  {
    ind <- indlist[[j]]
    bg <- as.integer(position[ind])
    ## for negative inserts bg = bg+35
    bg[insert[ind]<0] <- bg[insert[ind]<0]+35
    ed <- bg + as.integer(insert[ind])
    qu <- as.numeric(as.character(qual[ind]))
    squ <- as.numeric(as.character(sing.qual[ind]))
    n <- length(bg)
    if (n != 0)  {
      for (k in 1:n)  {
        if (bg[k]/s >= plot.bg & bg[k]/s <= plot.end & ed[k]/s >= plot.bg & ed[k]/s <= plot.end)  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((ed[k]-35)/s, ed[k]/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((bg[k]+35)/s, (ed[k]-35)/s), c(qu[k],qu[k]), lwd=0.1, col=mycols[j])
        }
      }
    }
  }
  ## deletion
  lines(c(del.bg/s, del.ed/s), c(-2,-2), col="black", lwd=2)
}


##### Plot all reads within 1kb of a deletion, library specific

# del.bg = deletion start
# del.ed = deletion end
# position = vector of start positions of reads
# insert = vector of insert sizes
# type = vector of paired flags
# sing.qual = vector of single mapping qualities 
# qual = vector of mapping qualities
# lib = vector indicating libraries
# thresh = cut-off for anomalous read pairs (two cutoffs for two libraries here)
# chr = chromosome

plot.reads.lib <- function(del.bg, del.ed, position, insert, type, sing.qual, qual, lib, thresh, chr)   {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  mycols <- c("darkgrey","darkred","orange","yellow2","green","blue","darkgreen","red")
  ## generate plot 
  plot.bg <- (del.bg - 1000)/s
  plot.end <- (del.ed + 1000)/s
  intvl <- plot.end - plot.bg
  m <- max(c(qual, sing.qual))
  plot(c(plot.bg - intvl/20, plot.end + intvl/3), c(-2,m), ty="n", xlab="position", ylab="quality scores", 
       cex.axis=1.3, cex.lab=1.3, cex.main=1.4, main=paste("chr",chr,":",del.bg,"-",del.ed," ", sep=""), las=1)
  x <- plot.end + intvl/10
  legend(x, m, legend=c("normal","FF","FR","RF","RR","diff chr","single","S-W","deletion"),
         lty=1, lwd=2, cex=0.8, col=c("darkgrey","darkred","red","orange","yellow2","green","blue","darkgreen","black"))
  ## find types
  ind1 <- which(type==1)
  ind2 <- c(which(lib==1 & (type==2 | (type==18 & insert > thresh[1]))),   # lib specific thresholds
            which(lib==2 & (type==2 | (type==18 & insert > thresh[2]))))
  ind4 <- which(type==4)
  ind8 <- which(type==8)
  ind18 <- c(which(lib==1 & type==18 & insert <= thresh[1]),               # lib specific thresholds
             which(lib==2 & type==18 & insert <= thresh[2]))
  ind32 <- which(type==32)      # different chr
  ind64 <- which(type==64)      # one read not mapped
  ind130 <- which(type==130)    # S-W alignment
  ## plot read pairs for all types
  indlist <- list(ind18, ind1, ind4, ind8, ind32, ind64, ind130, ind2)
  for (j in 1:length(indlist))  {
    ind <- indlist[[j]]
    bg <- as.integer(position[ind])
    ed <- bg + as.integer(insert[ind])
    qu <- qual[ind]
    squ <- sing.qual[ind]
    n <- length(bg)
    if (n != 0)  {
      for (k in 1:n)  {
        if (j==1 | j==2 | j==3 | j==4 | j==5 | j==8)  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((ed[k]-35)/s, ed[k]/s), c(qu[k],qu[k]), lwd=2, col=mycols[j])
          lines(c((bg[k]+35)/s, (ed[k]-35)/s), c(qu[k],qu[k]), lwd=0.1, col=mycols[j])
        }  else  {
          lines(c(bg[k]/s, (bg[k]+35)/s), c(squ[k],squ[k]), lwd=2, col=mycols[j])
        }
      }
    }
  }
  ## deletion
  lines(c(del.bg/s, del.ed/s), c(-2,-2), col="black", lwd=2)
}

### Plot deletion cluster
# x1 = start of deletion
# x2 = end of deletion
# y = depth

plot.cluster.basic1 <- function(x1, x2, y, read=100)  {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  ## generate plot 
  bg <- min(x1, x2)/s 
  ed <- max(x1, x2)/s 
  plot(c(bg, ed), c(-0.5, max(y)+0.5), ty="n", xlab="position", ylab="",
       cex.axis=1.3, cex.lab=1.3, main="", las=1, yaxt="n")
  n <- length(x1)
  for (k in 1:n)  {
    yk <- jitter(y[k])
    lines(c(x1[k]/s, (x1[k]+read)/s), c(yk,yk), lwd=3)
    lines(c((x2[k]-read)/s, x2[k]/s), c(yk,yk), lwd=3)
    lines(c((x1[k]+read)/s, (x2[k]-read)/s), c(yk,yk), lwd=0.3)
  }
}

plot.cluster.basic2 <- function(x1, x2, y)  {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  ## generate plot 
  bg <- min(x1, x2)/s 
  ed <- max(x1, x2)/s 
  plot(c(bg, ed), c(-0.5, max(y)+0.5), ty="n", xlab="position", ylab="index",
       cex.axis=1.3, cex.lab=1.3, main="", las=1)
  n <- length(x1)
  for (k in 1:n)  {
    lines(c(x1[k]/s, x2[k]/s), c(y[k],y[k]), lwd=0.3)
  }
}


plot.cluster.basic3 <- function(x1, x2, y, bounds, mycols)  {
  par(mai=c(0.5,0.5,0.3,0.2))
  s <- 1000000    ## scaling factor
  ## generate plot 
  bg <- min(x1, x2)/s 
  ed <- max(x1, x2)/s
  b <- bounds/s
  plot(c(bg, ed), c(-0.5, max(y)+0.5), ty="n", xlab="position", ylab="index",
       cex.axis=1.3, cex.lab=1.3, main="", las=1)
  n <- length(x1)

  j <- 1
  for (i in 1:nrow(bounds)) {
    ind <- which(x1/s >= b[i,1] & x2/s <= b[i,2])
    for (k in ind)  {
      lines(c(x1[k]/s, x2[k]/s), c(y[k],y[k]), lwd=3, col=mycols[j])
    }
    j <- j+1
  }
}


### Function that places empty spaces in front of each read to
### account for the start position
# pos = position
# seqs = read sequence

view.seq <- function(pos, seqs)  {
  n <- length(pos)
  pos <- c(0, cumsum(pos[2:n] - pos[1:(n-1)])) + 1
  rpos <- sapply(pos, rep, x=" ")
  rpos <- sapply(rpos, paste, collapse="")
  seq <- paste(rpos, seqs)
  seq
}

view.bpt.seq <- function(pos, seqs, bpt)  {
  n <- length(pos)
    if (n==1)  {
    rpos <- " "
  }  else  {
    pos <- c(0, cumsum(pos[2:n] - pos[1:(n-1)])) + 1
    rpos <- sapply(pos, rep, x=" ")
    rpos <- sapply(rpos, paste, collapse="")
  }
  seqs <- paste(rpos, seqs)
  ind <- which(bpt < nchar(seqs))
  bptseqs <- seqs
  if (length(ind)!=0)  
    bptseqs[ind] <- as.vector(unlist(sapply(seqs[ind], function(v) paste(substr(v,1,bpt), "|", substr(v, bpt+1, nchar(v))))))
  bptseqs
}

view.bpt.seq.sample <- function(pos, seqs, bpt, sample)  {
  n <- length(pos)
  if (n==1)  {
    rpos <- " "
  }  else  {
    pos <- c(0, cumsum(pos[2:n] - pos[1:(n-1)])) + 1
    rpos <- sapply(pos, rep, x=" ")
    rpos <- sapply(rpos, paste, collapse="")
  }
  seqs <- paste(rpos, seqs, sep="")
  ind <- which(bpt <= nchar(seqs)+1)
  bptseqs <- seqs
  if (length(ind)!=0)  
    bptseqs[ind] <- as.vector(unlist(sapply(seqs[ind], function(v) paste(substr(v,1,bpt), "|", substr(v, bpt+1, nchar(v))))))
  len <- nchar(bptseqs)
  m <- max(len) + 3
  blanks <- sapply(m-len, rep, x=" ")
  if (is.matrix(blanks))  {
    rblanks <- apply(blanks, 2, paste, collapse="")
  }  else  {
    rblanks <- sapply(blanks, paste, collapse="")
  }
  bptseqsample <- paste(bptseqs, rblanks, sample, sep=" ")
  as.matrix(bptseqsample)
}


### Get TP and FP rates for different mapping quality scores
# valid = validation matrix 0,1 or vector
# score = sum of mapping qualities
# s = steps for quantiles

get.tp.fp <- function(valid, score, s)  {
  qu <- round(quantile(score, probs = seq(0, 1, s)))
  tp <- tpp <- NULL
  k <- 1
  for (i in qu[-length(qu)])  {
    ind <- which(score>=i)
    if (is.matrix(valid) | is.data.frame(valid))  {
      tp[k] <- length(which(apply(valid[ind,], 1, function(v) any(v==1))))
      tpp[k] <- round(tp[k] / nrow(valid[ind,]),4)
    }  else  {
      tp[k] <- length(which(valid[ind]==1))
      tpp[k] <- round(tp[k] / length(valid[ind]),4)
    }
    k <- k+1
  }
  mq.cutoff <- data.frame(qu[-length(qu)], tp, 1-tpp)
  list(ms=qu[-length(qu)], tp=tp, fpp=1-tpp)
}


### Get TP and FP rates for different mapping quality scores using the union of both data sets
# valid = validation matrix 0,1 or vector
# score = sum of mapping qualities
# s = steps for quantiles

get.tp.fp.union <- function(valid, score, s, n)  {
  qu <- round(quantile(score, probs = seq(0, 1, s)))
  tp <- tpp <- NULL
  k <- 1
  for (i in qu[-length(qu)])  {
    ind <- which(score>=i)
    if (is.matrix(valid) | is.data.frame(valid))  {
      tp[k] <- length(which(apply(valid[ind,], 1, function(v) any(v==1))))
      tpp[k] <- round(tp[k] / (nrow(valid[ind,]) + n),4)
    }  else  {
      tp[k] <- length(which(valid[ind]==1))
      tpp[k] <- round(tp[k] / (length(valid[ind]) + n),4)
    }
    k <- k+1
  }
  mq.cutoff <- data.frame(qu[-length(qu)], tp, 1-tpp)
  list(ms=qu[-length(qu)], tp=tp, fpp=1-tpp)
}


### Convert matrix format to fasta

tab2fasta <- function (fn)
{
  ind <- c(1:nrow(fn))
  h.ind <- 2*ind-1
  s.ind <- 2*ind
  header <- as.character(fn[,1])
  sequence <- as.character(fn[,2])
  newfn <- NULL
  newfn[h.ind] <- header
  newfn[s.ind] <- sequence
  newfn
}

### Convert matrix format to fastq

tab2fastq <- function (fn)
{
  ind <- c(1:nrow(fn))
  h.ind <- 4*ind-3
  s.ind <- 4*ind-2
  p.ind <- 4*ind-1
  q.ind <- 4*ind
  header <- as.character(fn[,1])
  sequence <- as.character(fn[,2])
  quality <- as.character(fn[,3])
  newfn <- NULL
  newfn[h.ind] <- header
  newfn[s.ind] <- sequence
  newfn[p.ind] <- rep("+",length(ind))
  newfn[q.ind] <- quality
  newfn
}


##### Functions for genotyping #####

### Prepare and load the data
# fn = file name

load.data <- function(fn, quality=TRUE, cols="1-10,13")  {
  ## load reads and quality scores
  com <- paste("cut -f", cols, " ", fn, " > ~/bam_tmp.txt", sep="")
  ## com <- paste("cut -f1-10", fn, "> ~/bam_tmp.txt")
  system(com)
  bam <- "~/bam_tmp.txt"
  if (file.info(bam)$size != 0)  {
    f <- read.delim(bam, header=F, stringsAsFactor=FALSE)
    if (quality)  {
      com <- paste("cut -f11", fn, "> ~/qualities_tmp.txt")
      system(com)    
      qu <- scan("~/qualities_tmp.txt", what="character", na.strings="", sep="\n", quiet=TRUE)
      list(f=f, qu=qu)
    }  else {
      list(f=f)
    }
  }  else  {
    f <- NULL
    list(f=f)
  } 
}

load.data.washu <- function(flag, rcols, i, path)  {
  sam <- paste(path, flag, "_tmp", i,".sam", sep="")
  if (file.info(sam)$size != 0)  {
    com1 <- paste("cut -f", rcols, " ", sam," | sed 's/MF:i://' | sed 's/Aq:i://' | sed 's/AM:i://' > ", path, flag,"_reads", i,".sam", sep="") 
    com2 <- paste("cut -f11 ", sam, " > ", path, flag, "_quals", i,".sam", sep="")
    system(com1)
    system(com2)
    inreads <- paste(path, flag, "_reads", i,".sam", sep="")
    inquals <- paste(path, flag, "_quals", i,".sam", sep="")
    f <- read.delim(inreads, header=F, stringsAsFactor=FALSE)
    qu <- scan(inquals, what="character", na.strings="", sep="\n", quiet=TRUE)
    list(f=f, qu=qu)
  } 
}

### sort directory numerically by chr
# d = directory

sort.by.chr <- function(dd, col1=9, col2=10)  {
  chr <- vecstr2matstr(dd, "_")[,col1]
  chr <- gsub("chr","", chr)
  chr[chr=="X"] <- 23
  chr <- as.numeric(chr)
  start <- vecstr2matstr(dd, "_")[,col2]
  start <- as.numeric(vecstr2matstr(start, "-")[,1])
  ord <- order(chr, start)
  dd[ord]
}

## Select cutoff as the read length so that read goes across breakpoint
# f = matrix that contains reads, etc
# qu = vector of quality scores
# bpos = breakpoint position
# readlen = vector of read lengths
# microhom = sequence of microhomologies
# pflag = paired flag (col 11, if qualities are removed)

select.split.reads <- function(f, qu, bpos, readlen, microhom, left, pflag)  {
  ## select only mapped reads
  indp <- which(f[,pflag]!="MF:i:192")
  f <- f[indp,]
  qu <- qu[indp]
  readlen <- readlen[indp]
  ## length of microhomology at breakpoint
  microlen <- nchar(as.character(microhom))
  ## select reads within read length of breakpoint
  if (left)  {
    ind <- which(bpos-f[,4] + 1 < (readlen-microlen) & (bpos-f[,4]) > 0)      # left breakpoint
  }  else  {
    ind <- which(bpos-f[,4] + 1 < readlen & (bpos-f[,4]) > microlen)          # right breakpoint
  }
  list(f=f[ind,], qu=qu[ind])
}

## without base qualities
select.split.reads2 <- function(f, bpos, readlen, microhom, left, pflag)  {
  ## select only mapped reads
  indp <- which(f[,pflag]!="MF:i:192")
  f <- f[indp,]
  readlen <- readlen[indp]
  ## length of microhomology at breakpoint
  microlen <- nchar(as.character(microhom))
  ## select reads within read length of breakpoint
  if (left)  {
    ind <- which(bpos-f[,4] + 1 < (readlen-microlen) & (bpos-f[,4]) > 0)      # left breakpoint
  }  else  {
    ind <- which(bpos-f[,4] + 1 < readlen & (bpos-f[,4]) > microlen)          # right breakpoint
  }
  f[ind,]
}

### Select RPs within threshold around breakpoints
# f = matrix that contains reads, etc
# bpos.left = left breakpoint position
# bpos.right = right breakpoint position
# readlen = vector of read lengths
# pflag = paired flag (col 11, if qualities were removed)
# thresh = threshold of median +/- 5*sd

select.rps <- function(f, bpos.left, bpos.right, readlen, pflag, thresh)  {
  ## select only mapped reads
  indp <- which(f[,pflag]!="MF:i:192")
  f <- f[indp,]
  readlen <- readlen[indp]
  ## left breakpoint within RPs
  ind1 <- f[,4] > (bpos.left-thresh) & f[,9] > (bpos.left-f[,4]+readlen) & f[,4] < (bpos.left-readlen)
  ## right breakpoint within RPs
  ind2 <- f[,4] > (bpos.right-thresh) & f[,9] > (bpos.right-f[,4]+readlen) & f[,4] < (bpos.right-readlen)
  ## RPs within breakpoints - maybe better to use model for read depth
  # ind3 <- f[,4] > bpos.left & f[,4] < (bpos.right-readlen)
  ## RPs spanning both breakpoints
  ind4 <- f[,4] > (bpos.left-thresh) & f[,4] < bpos.left-readlen & f[,9] > (bpos.right-f[,4]+readlen)
  ## single reads around breakpoints
  #ind5 <- (f[,pflag]=="MF:i:32" | f[,pflag]=="MF:i:64") & f[,4] > (bpos.left-thresh) & f[,4] < (bpos.right+thresh)
  #ind <- ind1 | ind2 | ind3 | ind4 | ind5
  ind <- ind1 | ind2 | ind4
  f[ind,]
}

select.rps2 <- function(f, bpos, readlen, pflag, thresh, microhom)  {
  ## length of microhomology at breakpoint
  microlen <- nchar(as.character(microhom))
  ## select only mapped reads
  indp <- which(f[,pflag]!="MF:i:192")
  f <- f[indp,]
  readlen <- readlen[indp]
  ## breakpoint within RPs
  ind <- f[,4] > (bpos-thresh) & f[,9] > (bpos-f[,4]+readlen) & f[,4] < (bpos-readlen+microlen)
  f[ind,]
}

### Function that returns index of match and of mismatches for a given read and its mapping position
# r = read
# rpos = position of read in the genome
# bpos = break-point position of read in the genome
# ref = appropriate reference interval
# left = left breakpoint
# shift = distance between breakpoint and start of reference interval

which.match <- function(r, pos, ref, bmatch=TRUE)  {
  r <- as.character(r)
  n <- nchar(r)
  sref <- substr(ref, pos, pos+n-1)
  vref <- unlist(strsplit(sref, ""))
  ## F strand
  vr <- unlist(strsplit(r, ""))
  indm1 <- which(vr==vref)
  ## R strand
  vrr <- seq2revcompl(vr)
  indm2 <- which(vrr==vref)
  if (length(indm1) >= length(indm2))  {
    indm <- indm1
  }  else  {
    indm <- indm2
  }
  indmm <- setdiff(1:n, indm)
  if (bmatch) {
    ind <- indm
  }  else  {
    ind <- indmm
  }
  ind
}

  
### Collect all reads that support either reference or alternative reference
# l.bpos.ref = left breakpoint position on reference
# r.bpos.ref = right breakpoint position on reference
# l.bpos.alt = left breakpoint position on alternative reference
# r.bpos.alt = right breakpoint position on alternative reference
# reads.ref = reads that mapped against the reference genome
# reads.alt = reads that mapped against the alternative reference genome
# shift = length of sequence on left side of breakpoint
# microhom = microhomology at breakpoint

collect.reads <- function(l.bpos.ref, r.bpos.ref, l.bpos.alt, r.bpos.alt, reads.ref, reads.alt,
                          shift=100, microhom, pflag)  {
  ## get read lengths
  rreadlen <- nchar(as.character(reads.ref[,10]))
  areadlen <- nchar(as.character(reads.alt[,10]))
  ## select reads  
  l.rreads <- select.split.reads2(f=reads.ref, bpos=l.bpos.ref, readlen=rreadlen, microhom, left=TRUE, pflag)
  r.rreads <- select.split.reads2(f=reads.ref, bpos=r.bpos.ref, readlen=rreadlen, microhom, left=FALSE, pflag)
  l.areads <- select.split.reads2(f=reads.alt, bpos=l.bpos.alt, readlen=areadlen, microhom, left=TRUE, pflag)
  r.areads <- select.split.reads2(f=reads.alt, bpos=r.bpos.alt, readlen=areadlen, microhom, left=FALSE, pflag)
  ## adjust for start position on ~200bp selected reference
  l.rreads[,4] <- shift - (l.bpos.ref - as.numeric(l.rreads[,4]))
  r.rreads[,4] <- shift - (r.bpos.ref - as.numeric(r.rreads[,4])) + nchar(as.character(microhom)) + 1
  l.areads[,4] <- shift - (l.bpos.alt - as.numeric(l.areads[,4]))
  r.areads[,4] <- shift - (r.bpos.alt - as.numeric(r.areads[,4])) + nchar(as.character(microhom)) + 1
  if (any(l.rreads$f[,4]<=0) | any(r.rreads$f[,4]<=0) | any(l.areads$f[,4]<=0) | any(r.areads$f[,4]<=0))
    print("reference sequence not adequate")
  ## remove double entries
  allreads <- rbind(l.rreads, r.rreads, l.areads, r.areads)
  if (nrow(allreads)>1)  {
    id <- apply(allreads[,c(1,3,4)],1,paste,collapse="_")
    ord <- order(id)
    allreads.df <- data.frame(id[ord], allreads[ord,])
    freads <- get.unique.df(allreads.df, 1)[,-1]
  }  else  {
    freads <- allreads
  }
  freads
}
 
collect.reads.quality <- function(l.bpos.ref, r.bpos.ref, l.bpos.alt, r.bpos.alt, reads.ref, reads.alt,
                          q.ref, q.alt, shift=100, microhom, pflag)  {
  ## get read lengths
  rreadlen <- nchar(as.character(reads.ref[,10]))
  areadlen <- nchar(as.character(reads.alt[,10]))
  ## select reads
  l.rreads <- select.split.reads(f=reads.ref, qu=q.ref, l.bpos.ref, rreadlen, microhom, left=TRUE, pflag)  ## generate 2x reads for alternative ref
  r.rreads <- select.split.reads(f=reads.ref, qu=q.ref, r.bpos.ref, rreadlen, microhom, left=FALSE, pflag)
  l.areads <- select.split.reads(f=reads.alt, qu=q.alt, l.bpos.alt, areadlen, microhom, left=TRUE, pflag)
  r.areads <- select.split.reads(f=reads.alt, qu=q.alt, r.bpos.alt, areadlen, microhom, left=FALSE, pflag)
  ## adjust for start position on ~200bp selected reference
  l.rreads$f[,4] <- shift - (l.bpos.ref - as.numeric(l.rreads$f[,4]))
  r.rreads$f[,4] <- shift - (r.bpos.ref - as.numeric(r.rreads$f[,4])) + nchar(as.character(microhom)) + 1
  l.areads$f[,4] <- shift - (l.bpos.alt - as.numeric(l.areads$f[,4]))
  r.areads$f[,4] <- shift - (r.bpos.alt - as.numeric(r.areads$f[,4])) + nchar(as.character(microhom)) + 1
  ## check whether all positions are > 0
  if (any(l.rreads$f[,4]<=0) | any(r.rreads$f[,4]<=0) | any(l.areads$f[,4]<=0) | any(r.areads$f[,4]<=0))
    print("reference sequence not adequate")
  ## collate all reads and base qualities
  allreads <- rbind(l.rreads$f, r.rreads$f, l.areads$f, r.areads$f)
  allquals <- c(l.rreads$qu, r.rreads$qu, l.areads$qu, r.areads$qu)
  ## remove double entries
  if (!is.null(allreads))  {
    if (nrow(allreads)>1)  {
      id <- apply(allreads[,c(1,3,4)],1,paste,collapse="_")
      ord <- order(id)
      allreads.df <- data.frame(id[ord], allreads[ord,])
      freads <- get.unique.df(allreads.df, 1)[,-1]
      ## select associated base qualities
      ind.unique <- get.unique.df.id(allreads.df, 1)
      fquals <- allquals[ord][ind.unique]
    }  else  {
      freads <- allreads
      fquals <- allquals
    }
    list(f=freads, qu=fquals)
  }  else  {
    list(f=NULL, qu=NULL)
  }
}
 

### Collect aberrant read pairs
 
collect.rps <- function(l.bpos, r.bpos, reads, pflag=11, thresh)  {
  ## get read lengths
  rreadlen <- nchar(as.character(reads[,10]))
  ## select reads
  rreads <- select.rps(f=reads, bpos.left=l.bpos, bpos.right=r.bpos, readlen=rreadlen, pflag, thresh) 
  rreads
}
 
### Collect aberrant read pairs + singletons that span the breakpoint in alt-ref

collect.rps2 <- function(l.bpos.ref, r.bpos.ref, l.bpos.alt, reads.ref, reads.alt, microhom, pflag=11, thresh)  {
  ## get read lengths
  rreadlen <- nchar(as.character(reads.ref[,10]))
  ## select reads
  rreads <- select.rps(f=reads.ref, bpos.left=l.bpos.ref, bpos.right=r.bpos.ref, readlen=rreadlen, pflag, thresh) 
  areads <- select.rps2(f=reads.alt, bpos=l.bpos.alt, readlen=rreadlen, pflag, thresh, microhom)
  ## select reads that are single in ref and span breakpoint in alt-ref
  inds <- rreads[,11]=="MF:i:32" | rreads[,11]=="MF:i:64"
  if (length(which(inds))!=0) {
    indr <- as.vector(unlist(sapply(as.character(areads[,1]), function(v) which(v==as.character(rreads[inds,1])))))
    if (length(indr)!=0)  {
      frreads <- rbind(rreads[!inds,], rreads[inds,][indr,])
    }  else  {
      frreads <- rreads[!inds,]
    }
  }  else  {
    frreads <- rreads[!inds,]
  }
  frreads
}
 
### Select RPs within breakpoints
# f = matrix that contains reads, etc
# bpos.left = left breakpoint position
# bpos.right = right breakpoint position
# readlen = vector of read lengths
# pflag = paired flag (col 11, if qualities were removed)

collect.rd.ref <- function(f, bpos.left, bpos.right, pflag)  {
  ## select only mapped reads or mapped single reads (MF:i:0) and mapping quality > 0
  indp <- f[,pflag]!="MF:i:192" & f[,pflag]!="MF:i:64" & f[,pflag]!="MF:i:32" & f[,5]>0
  f <- f[indp,]
  readlen <- nchar(as.character(f[,10]))
  ## proper RPs or single reads within breakpoints, each read separately
  ind <- which(f[,4] > bpos.left & f[,4] < (bpos.right-readlen))
  rpnum <- length(ind)
  list(rpnum=rpnum, readlen=round(mean(readlen)))
}

### Collect all reads that support either reference or alternative reference
# l.bpos.ref = left breakpoint position on reference
# r.bpos.ref = right breakpoint position on reference
# l.bpos.alt = left breakpoint position on alternative reference
# r.bpos.alt = right breakpoint position on alternative reference
# reads.ref = reads that mapped against the reference genome
# reads.alt = reads that mapped against the alternative reference genome
# l.seqs.ref = sequence around left breakpoint (100bp left of breakpoint)
# r.seqs.ref = sequence around right breakpoint (100bp left of breakpoint)
# seqs.alt =  alternative sequence around deletion breakpoint (100bp left of breakpoint)
# shift = length of sequence on left side of breakpoint
# microhom = microhomology at breakpoint

count.reads.sample <- function(l.bpos.ref, r.bpos.ref, l.bpos.alt, r.bpos.alt, reads.ref, reads.alt, 
                               l.seqs.ref, r.seqs.ref, seqs.alt, shift=100, microhom, pflag)  {
  counts <- rep(0,4)
  ## select reads
  freads <- collect.reads(l.bpos.ref, r.bpos.ref, l.bpos.alt, r.bpos.alt, reads.ref, reads.alt,
                          shift, microhom, pflag)
  if (nrow(freads)>0)  {
    n <- nrow(freads)
    thresh <- 0.75 * nchar(as.character(freads[,10]))
    ## get counts
    m <- NULL
    if (n!=0)  {
      for (i in 1:n)  {
        m[1] <- length(which.match(r=freads[i,10], pos=freads[i,4], ref=l.seqs.ref, bmatch=TRUE))
        m[2] <- length(which.match(r=freads[i,10], pos=freads[i,4], ref=seqs.alt, bmatch=TRUE))
        m[3] <- length(which.match(r=freads[i,10], pos=freads[i,4], ref=seqs.alt, bmatch=TRUE))
        m[4] <- length(which.match(r=freads[i,10], pos=freads[i,4], ref=r.seqs.ref, bmatch=TRUE))
        indm <- which(m==max(m))
        if (m[indm][1] > thresh[i])  {
          counts[indm] <- counts[indm] + 1
        }
      }
    }
  }
  counts
}

### Calculate likelihood for each read
# f = read
# qu = quality of read
# bpos = breakpoint position
# bseqs = breakpoint sequence
# shift = length of sequence on each side of breakpoint
# left = left breakpoint
# bmatch = calculate likelihood for matching bases

get.read.lik <- function(f, qu, bseq, shift, uniform)  {
  ## mapping quality
  mq <- f[,5]
  mq[mq==0] <- 0.000001   ## maybe needs to be changed
  if (uniform==FALSE)  {
    mprob <- 10^(-mq/10)
    ## convert base qualities from ascii to decimal and to probs
    q.dec <- ascii2dec(qu)
    bprobs <- 10^(-q.dec/10)
  }  else  {
    mprob <- 10^(-60/10)                   
    bprobs <- rep(10^(-60/10), nchar(qu))   
  }
  ## get matches and mismatches
  seqmatch <- which.match(r=f[,10], pos=f[,4], ref=bseq, bmatch=TRUE)
  mismatch <- which.match(r=f[,10], pos=f[,4], ref=bseq, bmatch=FALSE)
  ## P(z | x, u) probability of a read z given the reference x and position u
  lik <- sum(log(1-bprobs[seqmatch])) + sum(log(bprobs[mismatch]/3))
  ## R(z, u) = P(z | x, u) * pQ / (1-pQ) likelihood that the read matches somewhere else
  rlik <- lik + (log(mprob) - log(1 - mprob)) 
  ## P(z | T) likelihood of a read 
  exp(lik) + exp(rlik)
}


### Get the probability of a read pair given the insert size distribution
# rpdist = outer distance between RPs
# m = median insert size
# sde = standard deviation

get.rp.prob <- function(rpdist, m, sde, ref)  {
  if (rpdist!=0)  {
    p <- dnorm(rpdist, mean=m, sd=sde)
  }  else  {
    if (ref)  {          ## maybe change to only use proper pairs
      p <- 0.00001
    }  else  {
      p <- 0.1
    }
  }
  p
}

### Calculate likelihood for a read at a SV breakpoint, using left and right breakpoint simultaneously (clean up on 15/6/2010)
# f = read
# qu = quality of read
# bpos = breakpoint position
# bseqs = breakpoint sequence
# shift = length of sequence on each side of breakpoint
# left = left breakpoint
# bmatch = calculate likelihood for matching bases
# uniform = if TRUE then uniform base and mapping qualities
# w = weight for hom (and ref)

get.sr.lik <- function(f, qu, seqs.ref.left, seqs.ref.right, seqs.alt, shift, uniform, w)  {
  liks <- NULL
  if (!is.null(f))  {
    ## calculate hom and ref likelihood for split reads
    lik.hom <- get.read.lik(f, qu, bseq=seqs.alt, shift, uniform)
    lik.ref.left <- get.read.lik(f, qu, bseq=seqs.ref.left, shift, uniform)
    lik.ref.right <- get.read.lik(f, qu, bseq=seqs.ref.right, shift, uniform)
    lik.ref <- max(c(lik.ref.left, lik.ref.right), na.rm=T)
    ## P(z | T) likelihood of a read given the genotype (hom, het, ref)
    liks[1] <- lik.hom                           
    liks[2] <- w*lik.hom + (1-w)*lik.ref
    liks[3] <- lik.ref
  }  else  {
    liks[1] <- 1/3     ## changed NA to log(0.333) on 9/2/2010
    liks[2] <- 1/3     ## changed to 0.333 on 15/6/2010
    liks[3] <- 1/3
  }
  n <- sum(liks, na.rm=TRUE)
  if (abs(n) > 1e-100)  {         ## changed on 17/6/2010, maybe not ideal choice yet
    probs <- liks/n
  }  else  {
    print("sum SR liks == 0")
    probs <- c(1/3,1/3,1/3)        ## added on 15/6/2010
  }
  probs  
}

### Calculate likelhood for RPs (clean up on 15/6/2010)

get.rp.lik <- function(f, svlen, m, sde, w)  {
  liks <- NULL
  rpdist <- f[,9]
  #if (!is.null(f))  {
  if (length(rpdist!=0))  {
    ## calculate hom and ref likelihood for read pairs
##     lik.hom <- get.rp.prob(rpdist, m+svlen, sde, ref=FALSE)
##     lik.ref <- get.rp.prob(rpdist, m, sde, ref=TRUE)
    dsde <- 2*sde
    lik.hom <- get.rp.prob(rpdist, m+svlen, dsde, ref=FALSE)     ## changed 19/6/2010, use greater standard deviation
    lik.ref <- get.rp.prob(rpdist, m, dsde, ref=TRUE)
    ## P(z | T) likelihood of a read pair given the genotype (hom, het, ref)
    liks[1] <- lik.hom                           
    liks[2] <- w*lik.hom + (1-w)*lik.ref
    liks[3] <- lik.ref
  }  else  {
    liks[1] <- 1/3     ## changed NA to log(0.333) on 9/2/2010
    liks[2] <- 1/3     ## changed to 0.333 on 15/6/2010
    liks[3] <- 1/3
  }
  n <- sum(liks, na.rm=TRUE)
  if (abs(n) > 1e-100)  {         ## changed on 17/6/2010
    probs <- liks/n
  }  else  {
    print("sum RP liks == 0")
    probs <- c(1/3,1/3,1/3)        ## added on 15/6/2010
  }
  probs  
}

### Get the probability of number of observed RPs within breakpoints given mapping coverage (clean up on 15/6/2010)

# mcov = mapping coverage
# l.bpos = left breakpoint
# r.bpos = right breakpoint
# k =  number of RPs within breakpoints
# readlen = read length

get.rd.lik <- function(mcov, l.bpos, r.bpos, k, readlen)  {
  ## compensate for not selecting RPs that cross the breakpoint
  svsize <- r.bpos - l.bpos + 1 - readlen               ## changed from 2*readlen to readlen - 15/6/2010
  ## lambda <- mcov/readlen
  lambda <- mcov/readlen   ## changed 19/6/2010 to account for overdispersion, changed 20/6/2010 back
  ## calculate hom and ref likelihood for read pairs
  if (svsize > 100) {                                     ## changed to svsize > 0 on 5/5/2010, changed to svsize > 100 on 15/6/2010
    lik.hom <- 0.000001^k
    lik.het <- dpois(k, (lambda/2)*svsize)    ## logs and normalise
    lik.ref <- dpois(k, lambda*svsize)
  }  else  {
    lik.hom <- 1/3
    lik.het <- 1/3
    lik.ref <- 1/3
  }
  ## P(z | T) likelihood of a read pair given the genotype (hom, het, ref)
  liks <- c(lik.hom, lik.het, lik.ref)
  n <- sum(liks, na.rm=TRUE)
  if (abs(n) > 1e-100)  {        ## changed on 17/6/2010
    probs <- liks/n
  }  else  {
    print("sum RD liks == 0")
    probs <- c(1/3,1/3,1/3)       ## added on 15/6/2010
  }
  probs  
}


### Combine SR, RP and RD liks (clean up on 15/6/2010) ###

get.lik.sample.sr.rp.rd <- function(l.bpos.ref, r.bpos.ref, l.bpos.alt, r.bpos.alt, reads.ref, reads.alt,
                           q.ref, q.alt, seqs.ref.left, seqs.ref.right, seqs.alt, shift, microhom,
                           uniform=FALSE, w=0.5, pflag=11, svlen, m, sde, thresh, mcov)  {
  probs <- NULL
  ## collect SR evidence across breakpoints
  reads <- collect.reads.quality(l.bpos.ref, r.bpos.ref, l.bpos.alt, r.bpos.alt, reads.ref, reads.alt,
                                  q.ref, q.alt, shift, microhom, pflag)
  ## collect RP evidence around breakpoints                                  ## needs to be tested !!!!!
  rps <- collect.rps(l.bpos.ref, r.bpos.ref, reads.ref, pflag, thresh)   ## change from rps to rps.ref (2010-05-27)
  if (is.null(rps))  {
    n2 <- 0
  }  else  {
    n2 <- nrow(rps)
  }
  ## collect RD evidence within breakpoints
  rd <- collect.rd.ref(f=reads.ref, l.bpos.ref, r.bpos.ref, pflag) 
  ## collate likelihoods
  liks1 <- liks2 <- NULL
  if (!is.null(reads))  {         ## moved down from above 20100615
    n1 <- length(reads$qu)
    ## get likelihood from SR
    if (n1!=0)  {
      for (i in 1:n1)  {
        liks1 <- rbind(liks1, get.sr.lik(f=reads$f[i,], qu=reads$qu[i], seqs.ref.left, seqs.ref.right, seqs.alt, shift, uniform, w))
      }
      lik1 <- apply(liks1, 2, prod, na.rm=TRUE)
      loglik1 <- log10(lik1/sum(lik1))              ## added on 17/6/2010
    }  else  {
      loglik1 <- c(-log10(3),-log10(3),-log10(3))    ## added on 4/6/2010
    }
  }  else  {
    logliks <- c(-log10(3), -log10(3), -log10(3))    ## moved up from below 20100615   ## added on 6/5/2010  
  }
  ## get likelihood from RP
  if (n2!=0)  {
    for (i in 1:n2)  {
      liks2 <- rbind(liks2, get.rp.lik(f=rps[i,], svlen, m, sde, w)) 
    }
    lik2 <- apply(liks2, 2, prod, na.rm=TRUE)
    loglik2 <- log10(lik2/sum(lik2))
  }  else  {
    loglik2 <- c(-log10(3),-log10(3),-log10(3))    ## added on 4/6/2010
  }
  ## get likelihood from read depth
  if (!is.na(rd$readlen))  {
    lik3 <- get.rd.lik(mcov,  l.bpos=l.bpos.ref, r.bpos=r.bpos.ref, k=rd$rpnum, readlen=rd$readlen)
    loglik3 <- log10(lik3)
  }  else  {
    loglik3 <- c(-log10(3),-log10(3),-log10(3))    ## added on 4/6/2010
  }
  ## combine likelihoods
  all.logliks <- rbind(loglik1, loglik2, loglik3)
  if (!is.null(all.logliks))  {
    logliks <- apply(all.logliks, 2, sum, na.rm=TRUE)   ## discard NA info from e.g. RD (5/5/2010)
    ## if all.logliks == 0
    inf <- logliks == -Inf                  ## changed is.infinite(logliks) to logliks == -Inf and moved it within IF clause (5/5/2010)
    if (any(inf))   print("warning: -Inf")
    logliks[inf] <- -800
  }  else  {
    logliks <- c(-log10(3), -log10(3), -log10(3))     ## changed NA to log10(0.333) on 9/2/2010  
  }
  list(logliks=logliks, loglik1=loglik1, loglik2=loglik2, loglik3=loglik3)
}


get.lik.sample.rd <- function(l.bpos.ref, r.bpos.ref, l.bpos.alt, r.bpos.alt, reads.ref, reads.alt,
                           q.ref, q.alt, seqs.ref.left, seqs.ref.right, seqs.alt, shift, microhom,
                           uniform=FALSE, w=0.5, pflag=11, svlen, m, sde, thresh, mcov)  {
  ## collect RD evidence within breakpoints
  rd <- collect.rd.ref(f=reads.ref, l.bpos.ref, r.bpos.ref, pflag) 
  ## get likelihood from read depth
  if (!is.na(rd$readlen))  {
    lik3 <- get.rd.lik(mcov, l.bpos=l.bpos.ref, r.bpos=r.bpos.ref, k=rd$rpnum, readlen=rd$readlen)
    print(lik3)
    loglik3 <- log10(lik3)
    print(loglik3)
  }  else  {
    loglik3 <- c(-log10(3),-log10(3),-log10(3))    ## added on 4/6/2010
  }
  loglik3
}




### Get genotype as max probs
# probs = genotypes likelihoods

get.genotype <- function(probs)  {
  n <- nrow(probs)
  len <- apply(probs, 1, function(v) length(unique(v)))
  probs[len==1] <- NA
  dat <- apply(probs, 1, which.max)
  len <- sapply(dat, length)
  dat[len==0] <- NA 
  gt <- unlist(dat)-1
  gt
}

get.genotype2 <- function(liks)  {
  #liks <- apply(liks, 1, function(v) v[is.infinite(v)] <- -400)
  n <- nrow(liks)
  gt <- rep(NA, nrow(liks))
  ord <- t(apply(liks, 1, order))
  ## quotient of optimal gt and sub-optimal gt
  k <- apply(cbind(liks,ord), 1, function(v) v[1:3][v[4:6][3]]/v[1:3][v[4:6][2]])
  ## if all three liks are the same, then NA, otherwise which.max-1 equals genotype
  indc <- which(k!=1)                          ## changed from k<1 on 8-6-2010
  ind <- apply(liks[indc,], 1, which.max)
  gt[indc] <- unlist(ind)-1
  gt
}

### Get LOD score as the quotient of most likely and second most likely

get.lod.score <- function(liks)  {
  n <- nrow(liks)
  gt <- rep(NA, nrow(liks))
  ord <- t(apply(liks, 1, order))
  ## quotient of optimal gt and sub-optimal gt
  k <- apply(cbind(liks,ord), 1, function(v) v[1:3][v[4:6][3]]/v[1:3][v[4:6][2]])
  abs(k)
}

### Get percentage of concordance of likelihoods with true genotypes
# probs = genotypes likelihoods
# gt = true genotypes

get.concordance <- function(probs, gt, na.rm=TRUE)  {
  n <- nrow(probs)
  #probs <- round(t(apply(probs, 1, function(v) { v / sum(v) } )),2)
  len <- apply(probs, 1, function(v) length(unique(v)))
  probs[len==1] <- NA
  dat <- apply(probs, 1, which.max)
  len <- sapply(dat, length)
  dat[len==0] <- NA 
  gtt <- unlist(dat)-1
  if (na.rm==TRUE)  {
    indna <- is.na(gt) | is.na(gtt)
    tp <- round(100*length(which(gt[!indna]==gtt[!indna]))/length(gt[!indna]), 1)
  }  else  {
    tp <- round(100*length(which(gt==gtt))/length(gt), 1)
  }
  tp
}

get.concordance2 <- function(probs, gt, thresh=1.3)  {
  n <- nrow(probs)
  sprobs <- apply(probs, 1, sort)
  if (!is.list(sprobs)) {
    sprobs <- t(sprobs)
  }  else  {
    len <- sapply(sprobs, length)
    sprobs[len==0] <- list(c(NA, NA, NA))
    sprobs <- convert.list2mat(sprobs, cnum=3)
  }
  lodscore <- apply(sprobs, 1, function(v) v[2]/v[3])
  ## which index is max
  indm <- apply(probs, 1, which.max)
  len <- sapply(indm, length)
  indm[len==0] <- NA
  indna <- is.na(gtt)
  ## adjust index to fit genotype call
  gtt <- unlist(indm)-1
  indu <- !is.na(gtt) | lodscore > thresh
  ind <- indu & !is.na(indu)
  tp <- round(100*length(which(gt[ind]==gtt[ind]))/length(gtt[ind]), 1)
  tp
}

### Get sample index for discordant likelihoods
# probs = genotypes likelihoods
# gt = true genotypes

get.discordance <- function(probs, gt)  {
  n <- nrow(probs)
  len <- apply(probs, 1, function(v) length(unique(v)))
  probs[len==1] <- NA
  dat <- apply(probs, 1, which.max)
  len <- sapply(dat, length)
  dat[len==0] <- NA
  gtt <- unlist(dat)-1
  fp <- gt!=gtt 
  fp
}

### Normalise likelihoods
# probs = genotypes likelihoods

normalise.rows <- function(probs)  {
  round(t(apply(probs, 1, function(v) { v / sum(v) } )),2)
}

#### Compare genotype probs before and after imputation
## dat1 = genotype probs before imputation
## dat2 = genotype probs after imputation
## delid = deletion id (= rs number)

compare.ps <- function(dat1,dat2,delid) {
  pin <- as.matrix(dat1[delid == dat1[,1],-(1:3)])     # changed dat2 to dat1 ??
  pout <- as.matrix(dat2[delid == dat2[,1],-(1:3)])

  pin <- matrix(pin,ncol=3,byrow=TRUE)
  pinsc <- t(apply(pin,1,function(v) v/sum(v)))
  pout <- matrix(pout,ncol=3,byrow=TRUE)
  poutsc <- t(apply(pout,1,function(v) v/sum(v)))

  inout <- cbind(pinsc,poutsc)
  inout
}

#### Chi-square test for SNP genotypes
## df = data frame
## i = first row number
## j = secund row number

chi.compare <- function(df,i,j) {
  x <- as.numeric(df[i,-(1:3)])
  y <- as.numeric(df[j,-(1:3)])
  ind <- is.na(x) | is.na(y)
  kx <- length(unique(x[!ind]))
  ky <- length(unique(y[!ind]))
  if (kx > 1 & ky > 1)  {
    r <- chisq.test(x,y)
    s <- r$statistic
  }  else  {
    s <- NA
  }
  s
}

#### Table to store all chi-square tests for SNP genotypes
## df = data frame

chi.table <- function(df) {
  n <- nrow(df)
  m <- matrix(0,n,n)
  for (i in 1:n) {
    #print(paste("i =", i))
    for (j in (i):n) {
      #print(paste("j =", j))
      chival <- chi.compare(df,i,j)
      #cat(i,j,chival$statistic,"\n")
      if (!is.na(chival))  {
##         m[i,j] <- round(chival$statistic, 1)
##         m[j,i] <- round(chival$statistic, 1)
        m[i,j] <- round(chival, 1)
        m[j,i] <- round(chival, 1)
      }  else  {
        m[i,j] <- NA
        m[j,i] <- NA
      }
    }
  }
  m
}


####### Analyse genotyping ######

### Function that calculates overall concordance
## gts = genotype set that contains golden set and set to be evaluated
## indx = golden set
## indy = set to be evaluated
calc.concord <- function(gts)  {
  indx <- grep("NA[0-9]*.x", names(gts), perl=TRUE)
  indy <- grep("NA[0-9]*.y", names(gts), perl=TRUE)
  con <- NULL
  for (i in 1:nrow(gts))  {
    gtp <- factor(paste(gts[i,indx], gts[i,indy], sep=""), levels=c(00,01,02,10,11,12,20,21,22))
    tab <- table(gtp)
    mat <- matrix(tab, nrow=3, ncol=3, byrow=TRUE)
    con[i] <- round(100*sum(c(mat[1,1], mat[2,2], mat[3,3])) / sum(mat),1)
  }
  con
}

### Function that calculates the non-reference discordance
## gts = genotype set that contains golden set and set to be evaluated
## indx = golden set
## indy = set to be evaluated
non.ref.discord <- function(gts)  {
  indx <- grep("NA[0-9]*.x", names(gts), perl=TRUE)
  indy <- grep("NA[0-9]*.y", names(gts), perl=TRUE)
  nrd <- NULL
  for (i in 1:nrow(gts))  {
    gtp <- factor(paste(gts[i,indx], gts[i,indy], sep=""), levels=c(00,01,02,10,11,12,20,21,22))
    tab <- table(gtp)
    mat <- matrix(tab, nrow=3, ncol=3, byrow=TRUE)
    if (sum(mat[-3,-3])!=0)  {
      nrd[i] <- round(100*sum(c(mat[1,2:3], mat[1,c(1,3)], mat[3,1:2])) / (sum(mat) - mat[3,3]),1)
    }  else  {
      nrd[i] <- NA
    }
  }
  nrd
}

### Function that calculates the non-reference sensitivity
## gts = genotype set that contains golden set and set to be evaluated
## indx = golden set
## indy = set to be evaluated
non.ref.sens <- function(gts)  {
  indx <- grep("NA[0-9]*.x", names(gts), perl=TRUE)
  indy <- grep("NA[0-9]*.y", names(gts), perl=TRUE)
  nrd <- NULL
  for (i in 1:nrow(gts))  {
    gtp <- factor(paste(gts[i,indx], gts[i,indy], sep=""), levels=c(00,01,02,10,11,12,20,21,22))
    tab <- table(gtp)
    mat <- matrix(tab, nrow=3, ncol=3, byrow=TRUE)
    n <- sum(mat[1:3,1:2]) + sum(is.na(gts[i,indx]))
    if (n!=0)  {
      nrd[i] <- round(100*sum(c(mat[1,1:2], mat[2,1:2])) / n, 1)
    }  else  {
      nrd[i] <- NA
    }
  }
  nrd
}

### Function that calculates overall concordance
## gts = genotype set that contains golden set and set to be evaluated
## indx = golden set
## indy = set to be evaluated
calc.concord.sample <- function(gts)  {
  indx <- grep("NA[0-9]*.x", names(gts), perl=TRUE)
  indy <- grep("NA[0-9]*.y", names(gts), perl=TRUE)
  con <- NULL
  for (j in 1:length(indx))  {
    gtp <- factor(paste(gts[,indx][,j], gts[,indy][,j], sep=""), levels=c(00,01,02,10,11,12,20,21,22))
    tab <- table(gtp)
    mat <- matrix(tab, nrow=3, ncol=3, byrow=TRUE)
    con[j] <- round(100*sum(c(mat[1,1], mat[2,2], mat[3,3])) / sum(mat),1)
  }
  id <- gsub(".x", "", names(gts)[indx])
  dat <- data.frame(id, con)
  names(dat) <- c("sample","pcon")
  dat
}



### Function that selects entry regarding keyword from INFO in VCF format
## info = VCF INFO field
## keyword = keyword
get.vcf.info.all <- function(info, keyword, num=TRUE)  {
  keyword <- paste(keyword, "=", sep="")
  info.list <- vecstr2matstr(info, ";")
  ind <- sapply(info.list, grep, pattern=keyword)
  len <- sapply(ind, length)
  ind[len==0] <- NA
  ind <- unlist(ind)
  n <- length(info.list)
  v <- NULL
  for (i in 1:n)  {
    v[i] <- info.list[[i]][ind[i]]
  }
  v <- sapply(v, edit.info, keyword)
  if (num)  {
    v <- sapply(v, as.numeric)
  }
  as.vector(unlist(v))
}

### Function that selects entry regarding keyword from INFO in VCF format
## info = VCF INFO field
## keyword = keyword
get.vcf.info.old <- function(vcfinfo, keyword, num=TRUE)  {
  info <- vecstr2matstr(vcfinfo, ";")
  if (is.matrix(info) & nrow(info)!=1)  {
    info <- list(info)
  }
  info <- sapply(info, function(v) v[grep(keyword, v)])
  keyword <- paste(keyword, "=", sep="")
  info <- gsub(keyword, "", info)
  if (num)  {
    info <- as.numeric(info)
  }
  info
}

### Function that selects entry regarding keyword from INFO in VCF format
## info = VCF INFO field
## keyword = keyword
get.vcf.info.old <- function(vcfinfo, keyword, num=TRUE, islist=FALSE)  {
  info <- vecstr2matstr(vcfinfo, ";")
  if (is.matrix(info) & nrow(info)!=1)  {
    info <- list(info)
  }
  if (keyword=="END")  {
    keyword <- paste("^", keyword, "=", sep="")
  }  else  {
    keyword <- paste("^", keyword, "=", sep="")
  }
  info <- sapply(info, function(v) v[grep(keyword, v, perl=TRUE)])
  info <- gsub(keyword, "", info)
  if (islist)  {
    info <- vecstr2liststr(info, ",")
  }
  if (num)  {
    info <- as.numeric(info)
  }
  info
}

### Function that selects entry regarding keyword from INFO in VCF format
## info = VCF INFO field
## keyword = keyword
get.vcf.info <- function(vcfinfo, keyword, num=TRUE)  {
  info <- vecstr2matstr(vcfinfo, ";")
  if (is.matrix(info) & nrow(info)!=1)  {
    info <- list(info)
  }
  keyword <- paste("^", keyword, "=", sep="")
  info <- sapply(info, function(v) v[grep(keyword, v, perl=TRUE)])
  info <- gsub(keyword, "", info)
  if (num)  {
    info <- as.numeric(info)
  }
  info
}


### Function that edits the entry regarding the keyword, used in get.info
edit.info <- function(v, keyword)  {
  v <- gsub(keyword, "", v)
  ind <- grep(",", v)
  if (length(ind)!=0)  {
    v <- as.numeric(unlist(strsplit(v, ",")))
    v[is.na(v)] <- c(NA, NA)
  }
  v
}

### Get read depth for given window
# pos = positions
# win = window size
get.rd.win <- function(pos, win)  {
  n <- length(pos)
  a <- pos[1]
  b <- pos[n]
  s1 <- seq(a, b-win, win)
  s2 <- seq(a+win, b, win)
  rd <- apply(cbind(s1,s2), 1, function(v) sum(pos>=v[1] & pos<v[2]))
  list(rd=rd, pos=(s1+s2)/2)
}


### Get k-mers

get.kmer.probs <- function(s, m, j)  {
  k <- vecstr2matstr(s, " ")
  n <- lapply(k, function(v) as.numeric(vecstr2matstr(v,":")[,2]))
  s <- sapply(n, sum)
  p <- lapply(n, function(v) round(v/sum(v),4))
  t <- (1/m) + j * (1/m) * (1 - 1/m)
  r <- sapply(p, function(v) all(v<=t))
  list(n=n, p=p, r=r)
}


### Get FDR (Benjamini-Hochberg 1995)
## p = p-values
## q = FDR at level q
get.fdr <- function(p, q)  {
  m <- length(p)
  sp <- sort(p)
  k <- which(sp <= (q * 1:m / m))
  if (length(k)!=0)  {
    n <- max(k)
    p <- max(sp[k])
  }  else  {
    n <- 0
    p <- 0
  }
  list(n=n, p=p)
}
