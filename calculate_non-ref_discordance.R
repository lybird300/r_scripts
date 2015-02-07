
## This R script calculates the non-reference discordance of two sets of genotypes.
## One input file needs to be in VCF format --vcf=<input.vcf>, but just the row with the variant that is compared
## The other input file needs to be in PED format --ped=<input.ped>
## The script assumes that the samples are in the same order.
## Author: Klaudia Walter, 2012-08-01

## Usage
## R CMD BATCH --vanilla --slave --args --vcf=<input.vcf> --ped=<input.ped> calculate_non-ref_discordance.R



## Function that splits string by given pattern, returns vector
# string = string that needs splitting
# pat = pattern by which to split

str2vec <- function(string, pat) 
{
  charvec <- unlist(strsplit(string, pat))
  charvec
}

## Function that splits vector of strings by given pattern, returns matrix
# vecstring = vector of string 
# pat = pattern by which to split
# returns matrix/list that contains split up strings

vecstr2matstr <- function(vecstring, pat)
{
  matstring <- t(apply(as.matrix(vecstring), 1, str2vec, pat))
  matstring
}

## Reads input.vcf and input.ped
get.Args <- function() {
   ## Extract argument names and corresponding values
   my_args_list   <- strsplit(grep("=", gsub("--", "", commandArgs()), value = T), "=")
   my_args        <- lapply(my_args_list, function(x) x[2])
   names(my_args) <- lapply(my_args_list, function(x) x[1])
   ## Assign argument values to global variables in R
   for (arg_name in names(my_args)) {
     eval(parse(text=paste(arg_name, " <<- my_args[[arg_name]]", sep="")))
   }
}
get.Args()

## Read data
vcf.dat <- read.delim(vcf, stringsAsFactors=F, header=F)
ped.dat <- read.delim(ped, stringsAsFactors=F, header=F, sep=" ")

## Check number of samples
nvcf <- ncol(vcf.dat) - 9
nped <- nrow(ped.dat)
  
if (nvcf!=nped)
  print("Number of samples is not the same!")


### Get VCF genotypes
ref <- vcf.dat[,4]
alt <- vcf.dat[,5]
gt0 <- paste(ref, ref, sep="")
gt1 <- paste(ref, alt, sep="")
gt2 <- paste(alt, alt, sep="")
vgts <- vecstr2matstr(vcf.dat[,10:ncol(vcf.dat)], ":")[,1]
vgtsm <- matrix("N", nrow=length(vgts), ncol=2)
vgtsm[vgts=="0/0" | vgts=="0|0",1] <- ref
vgtsm[vgts=="0/0" | vgts=="0|0",2] <- ref
vgtsm[vgts=="0/1" | vgts=="0|1" | vgts=="1/0" | vgts=="1|0",1] <- ref
vgtsm[vgts=="0/1" | vgts=="0|1" | vgts=="1/0" | vgts=="1|0",2] <- alt
vgtsm[vgts=="1/1" | vgts=="1|1",1] <- alt
vgtsm[vgts=="1/1" | vgts=="1|1",2] <- alt

### Get Plink genotypes
pgtsm <- ped.dat[,7:8]

### Remove Ref/Ref calls that are shared and remove unknown genotypes
indrm <- (vgtsm[,1]==ref & vgtsm[,2]==ref & pgtsm[,1]==ref & pgtsm[,2]==ref) | (vgtsm[,1]=="N" | vgtsm[,2]=="N" | pgtsm[,1]=="N" | pgtsm[,2]=="N")
vgtsr <- vgtsm[!indrm,]
pgtsr <- pgtsm[!indrm,]

### Calculate non-reference discordance of remaining genotypes
indd <- (vgtsr[,1]==pgtsr[,1] & vgtsr[,2]==pgtsr[,2]) | (vgtsr[,1]==pgtsr[,2] & vgtsr[,2]==pgtsr[,1])
nrd <- round(100 * sum(!indd) / nrow(vgtsr),1)
print(paste("Non-reference discordance = ", nrd, "%", sep=""))


