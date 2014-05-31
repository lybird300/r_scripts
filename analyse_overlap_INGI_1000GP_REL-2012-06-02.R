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
rm(list=ls())
#isolates
population <- "FVG"
population <- "VBI"
i_pops <- c("FVG","VBI")
#outbred
population <- "UK10K"
population <- "TGP"
o_pops <- c("UK10K","TGP")

path <- "~/fvg_seq/variant_refinemet/annotations/UK10K/ANNOTATED_OVERLAP/MAF/"
path <- "~/fvg_seq/variant_refinemet/annotations/1TGP/ANNOTATED_OVERLAP/MAF/"
table_name <- paste(path,population,"_all_chr.csv",sep="")
pop_table <- read.table(table_name,header=T,stringsAsFactors=F, comment.char="")

#remove monomorphic sites
dim(pop_table)
pop_table_no_mono <- pop_table[-which(pop_table$AF ==1 | pop_table$AF == 0),]
dim(pop_table_no_mono)
pop_table_mono <- pop_table[which(pop_table$AF ==1 | pop_table$AF == 0),]
dim(pop_table_mono)

source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
maf_classes <- c(0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30)
all_maf_classes <- split_bins(maf_classes,pop_table_no_mono,population)

mode <- "inbred"

maf_shared <- NULL
for(i in 2:(length(maf_classes) + 1)){
  class_maf_name <- paste(population,'_maf_lte_',maf_classes[i],sep='')
  filename <- paste(class_maf_name,'table.txt',sep='_')
    
  if (maf_classes[i - 1] == maf_classes[length(maf_classes)]){
    if (maf_classes[length(maf_classes)] != 0.5) {        
      class_maf_name <- paste(population,'_maf_gte_',maf_classes[i-1],sep='')
      filename <- paste(class_maf_name,'table.txt',sep='_')
      maf_table <- read.table(filename,sep="\t",header=T,stringsAsFactors=F, comment.char="")
      i=i-1
    }else{break}
  }else{
    maf_table <- read.table(filename,sep="\t",header=T,stringsAsFactors=F, comment.char="")
  }
  #extract overlapping sites
  #modify according to which population we are using as reference to calculate the frequencies
  
  total <- length(maf_table$POS)
  total_pop_over <- NULL
  pop_head <- NULL

  if (mode == 'inbred'){
    sel_pops <- i_pops
  }else{
    sel_pops <- o_pops
  }

  for(pop in sel_pops){
    current_pop_name <- paste(pop,"maf",sep="_")
    current_over_name <- paste("over",pop,sep="_")
    current_p_name <- paste("p",pop,sep="_")

    current_pop_index <- which(colnames(maf_table) %in%  pop)
    current_pop_maf <- maf_table[which(maf_table[,current_pop_index]==1),]
    current_over <- length(current_pop_maf$POS)
    current_p <- round(100*current_over/total,2)
    # assign(current_pop_name,current_pop_maf)
    # assign(current_over_name,current_over)
    # assign(current_p_name,current_p)
    total_pop_over <- cbind(total_pop_over,current_over,current_p)
    pop_head <- c(pop_head,paste("overlap",pop,sep="_"),paste("overlap_",pop,"(%)",sep=""))
  }

  # over_tgp <- length(tgp_maf$POS)
  # over_uk10k <- length(uk10k_maf$POS)

  # p_tgp <- round(100*over_tgp/total,2)
  # p_uk10k <- round(100*over_uk10k/total,2)

  # maf_shared <- rbind(maf_shared,c(maf_classes[i],total,over_tgp,over_uk10k,p_tgp,p_uk10k))
  maf_shared <- rbind(maf_shared,c(maf_classes[i],total,total_pop_over))
}
maf_shared <- as.data.frame(maf_shared)
# colnames(maf_shared) <- c("MAF","total","overlap_tgp","overlap_uk10k","overlap_tgp(%)","overlap_uk10k(%)")

colnames(maf_shared) <- c("MAF","total",pop_head)

write.table(maf_shared,file=paste(population,"maf_shared.txt",sep="_"), row.names=F,col.names=T, quote=F, sep="\t")
##### Plots ######

unuk10k <- maf_shared$total - maf_shared$overlap_uk10k
ung1k <- maf_shared$total - maf_shared$overlap_tgp
xmax <- 3000000
whitebar <- xmax - (maf_shared$overlap_uk10k + unuk10k)
fscale <- 100000
sxmax <- xmax/fscale

pg1k <- paste(round(100 * maf_shared$overlap_tgp / maf_shared$total,1),"%", sep="")
puk10k <- paste(round(100 * maf_shared$overlap_uk10k / maf_shared$total,1),"%", sep="")

pdf(paste("shared_variants_uk10k_1000GP_",population,".pdf",sep=""), width=8, height=4, pointsize=12)
par(mfrow=c(1,2), mai=c(1,1.2,0.1,0.1))
slab <- c(0,5,10,15,20)
# labnames=c("singletons    ","doubletons   ","AF<=5%     ","AF<=10%    ","AF>10%     ")
# maf_classes <- c(0,0.005,0.01,0.02,0.05,0.10,0.20,0.30)
# maf_classes <- c(0,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30)
# labnames=c("AF<=0.005","1%  ","2% ","5% ","10%  ","20%  ","30%  ","AF>30% ")
labnames=c("MAF<=0.5%","1%  ","2% ","5% ","10%  ","15%  ","20%  ","25%  ","30%  ","MAF>30% ")
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

#############################################
#for local use
#isolates
population <- "VBI"
population <- "FVG"
i_pops <- c("FVG","VBI")
#outbred
population <- "UK10K"
population <- "TGP"
o_pops <- c("UK10K","TGP")

shared_table <- paste(population,"maf_shared.txt",sep="_")
maf_shared <- read.table(shared_table,header=T)

colnames(maf_shared) <- c("MAF","total","overlap_tgp","overlap_uk10k","overlap_tgp_perc","overlap_uk10k_perc")

############################################################################################################
#plot for each population with overlap with outbred population (there should be 2 plots: FVG vs 1000G and UK10K, VBI vs 1000G and UK10K)
#calculate the number of UNSHARED sites 
unuk10k <- maf_shared$total - maf_shared$overlap_uk10k
ung1k <- maf_shared$total - maf_shared$overlap_tgp
#set xaxis limit
xmax <- max(maf_shared$total)*2
#calculate the white gap to leave aside the bar
whitebar <- xmax - (maf_shared$overlap_uk10k + unuk10k)
fscale <- 100000
sxmax <- xmax/fscale

pg1k <- paste(round(100 * maf_shared$overlap_tgp / maf_shared$total,1),"%", sep="")
puk10k <- paste(round(100 * maf_shared$overlap_uk10k / maf_shared$total,1),"%", sep="")
# mycols <- c("lightskyblue","steelblue3","steelblue4")
mycols <- c("red","black")
pdf(paste("shared_variants_uk10k_1000GP_",population,".pdf",sep=""), width=8, height=4, pointsize=12)
  par(mfrow=c(1,2), mai=c(1,1.2,0.1,0.1))
  slab <- seq(0,round(sxmax/1.5),5)
  # slab <- c(0,5,10,15,20,25)
  labnames=c("MAF<=0.05%","1%  ","2% ","5% ","10%  ","15%  ","20%  ","25%  ","30%  ","MAF>30% ")
  #left-hand plot
  m <- barplot(-rbind(maf_shared$overlap_uk10k,unuk10k)/fscale,horiz=T,
              col=rev(mycols),
              xlab=paste("number of variants (x","100000",")"),
              axes=F,
              xlim=c(-sxmax+1,0),
              border=rev(mycols),
              )

  axis(1, labels=slab, at=-slab)
  if (population == "FVG"){
    text((-maf_shared$total/70000)-3, m, puk10k) #-> FVG plot
  }else if(population == "VBI"){
    text((-maf_shared$total/100000)-15, m, puk10k) #-> VBI plot
  }
  legend(x="left",y="center", legend=c(paste("unique",population),"shared"), fill=mycols, cex=0.8, border=mycols)
  mtext("D", side=3, line=0, at=xmax, cex=1.5)
  #right-hand plot
  barplot(rbind(maf_shared$overlap_tgp,ung1k)/fscale, horiz=T,
          names.arg=labnames,
          col=rev(mycols),
          las=1,
          xlab=paste("number of variants (x","100000",")"),
          axes=F, xlim=c(0,sxmax+1), border=rev(mycols))
  axis(1, labels=slab, at=slab)
  if (population == "FVG"){
    text((maf_shared$total/70000)+3, m, pg1k) #-> FVG plot
  }else if(population == "VBI"){
    text((maf_shared$total/100000)+15, m, pg1k) #-> VBI plot
  }
dev.off()

#######################################################################################################
#plot for each population with reciprocal overlap with outbred population (there should be 4 plots: FVG vs 1000G, FVG vs UK10K, VBI vs 1000G, VBI vs UK10K)
#isolates
# population <- "VBI"
# population <- "FVG"
i_pops <- c("FVG","VBI")
#outbred
o_pops <- c("UK10K","TGP")

inbred_shared_table <- paste(population,"maf_shared.txt",sep="_")
inbred_maf_shared <- read.table(inbred_shared_table,header=T)
inbred_maf_shared$MAF <- as.character(inbred_maf_shared$MAF)
inbred_maf_shared$MAF[10] <- ">0.3"
assign(paste(population,"maf_shared",sep="_"),inbred_maf_shared)

for(pop in o_pops){
  outbred_shared_table <- paste("/home/max/Work/OUTBRED/",pop,"_maf_shared.txt",sep="")
  outbred_maf_shared <- read.table(outbred_shared_table,header=T)
  outbred_maf_shared$MAF <- as.character(outbred_maf_shared$MAF)
  outbred_maf_shared$MAF[10] <- ">0.3"
  assign(paste(pop,"maf_shared",sep="_"),outbred_maf_shared)  
}

o_population <- "UK10K"
o_population <- "TGP"

if(population == "VBI"){
  VBI_TGP <- merge(VBI_maf_shared,TGP_maf_shared,by.x="MAF",by.y="MAF",sort=FALSE)
  colnames(VBI_TGP) <- c("MAF","tot_VBI","overlap_tgp","overlap_uk10k","overlap_tgp_perc","overlap_uk10k_perc","total_TGP","overlap_TGP_FVG","overlap_TGP_FVG_perc","overlap_TGP_VBI","overlap_TGP_VBI_perc")
  VBI_TGP_UK10K <- merge(VBI_TGP,UK10K_maf_shared,by.x="MAF",by.y="MAF",sort=FALSE)
  colnames(VBI_TGP_UK10K) <- c("MAF","tot_VBI","overlap_tgp","overlap_uk10k","overlap_tgp_perc","overlap_uk10k_perc","total_TGP","overlap_TGP_FVG","overlap_TGP_FVG_perc","overlap_TGP_VBI","overlap_TGP_VBI_perc","total_UK10K","overlap_UK10K_FVG","overlap_UK10K_FVG_perc","overlap_UK10K_VBI","overlap_UK10K_VBI_perc")
  
  if (o_population == "UK10K"){
    current_iso_un <- VBI_TGP_UK10K$tot_VBI - VBI_TGP_UK10K$overlap_uk10k
    current_out_un <- VBI_TGP_UK10K$total_UK10K - VBI_TGP_UK10K$overlap_UK10K_VBI
    #set xaxis limit
    xmax <- max(VBI_TGP_UK10K$total_UK10K)*1.5
    p_uk10k_iso <- paste(VBI_TGP_UK10K$overlap_uk10k_perc,"%", sep="")
    p_iso_uk10k <- paste(VBI_TGP_UK10K$overlap_UK10K_VBI_perc,"%", sep="")
  }else if(o_population=="TGP"){
    current_iso_un <- VBI_TGP_UK10K$tot_VBI - VBI_TGP_UK10K$overlap_tgp
    current_out_un <- VBI_TGP_UK10K$total_TGP - VBI_TGP_UK10K$overlap_TGP_VBI
    #set xaxis limit
    xmax <- max(VBI_TGP_UK10K$total_TGP)*1.5
    p_tgp_iso <- paste(VBI_TGP_UK10K$overlap_tgp_perc,"%", sep="")
    p_iso_tgp <- paste(VBI_TGP_UK10K$overlap_TGP_VBI_perc,"%", sep="")
    current_overlap <- VBI_TGP_UK10K
  }
  current_overlap <- VBI_TGP_UK10K


}else{
  FVG_TGP <- merge(FVG_maf_shared,TGP_maf_shared,by.x="MAF",by.y="MAF",sort=FALSE)
  colnames(FVG_TGP) <- c("MAF","tot_FVG","overlap_tgp","overlap_uk10k","overlap_tgp_perc","overlap_uk10k_perc","total_TGP","overlap_TGP_FVG","overlap_TGP_FVG_perc","overlap_TGP_VBI","overlap_TGP_VBI_perc")
  FVG_TGP_UK10K <- merge(FVG_TGP,UK10K_maf_shared,by.x="MAF",by.y="MAF",sort=FALSE)
  colnames(FVG_TGP_UK10K) <- c("MAF","tot_FVG","overlap_tgp","overlap_uk10k","overlap_tgp_perc","overlap_uk10k_perc","total_TGP","overlap_TGP_FVG","overlap_TGP_FVG_perc","overlap_TGP_VBI","overlap_TGP_VBI_perc","total_UK10K","overlap_UK10K_FVG","overlap_UK10K_FVG_perc","overlap_UK10K_VBI","overlap_UK10K_VBI_perc")
  if (o_population == "UK10K"){
    current_iso_un <- FVG_TGP_UK10K$tot_FVG - FVG_TGP_UK10K$overlap_uk10k
    current_out_un <- FVG_TGP_UK10K$total_UK10K - FVG_TGP_UK10K$overlap_UK10K_FVG
    #set xaxis limit
    xmax <- max(FVG_TGP_UK10K$total_UK10K)*1.5
    p_uk10k_iso <- paste(FVG_TGP_UK10K$overlap_uk10k_perc,"%", sep="")
    p_iso_uk10k <- paste(FVG_TGP_UK10K$overlap_UK10K_FVG_perc,"%", sep="")
  }else if(o_population=="TGP"){
    current_iso_un <- FVG_TGP_UK10K$tot_FVG - FVG_TGP_UK10K$overlap_tgp
    current_out_un <- FVG_TGP_UK10K$total_TGP - FVG_TGP_UK10K$overlap_TGP_FVG
    #set xaxis limit
    xmax <- max(FVG_TGP_UK10K$total_TGP)*1.5
    p_tgp_iso <- paste(FVG_TGP_UK10K$overlap_tgp_perc,"%", sep="")
    p_iso_tgp <- paste(FVG_TGP_UK10K$overlap_TGP_FVG_perc,"%", sep="")
  }
  current_overlap <- FVG_TGP_UK10K
}


#calculate the white gap to leave aside the bar
# whitebar <- xmax - (maf_shared$overlap_uk10k + unuk10k)
fscale <- 100000
sxmax <- xmax/fscale
mycols <- c("lightskyblue","steelblue3","steelblue4")
# mycols <- c("red","black")

####################################################
#Isolates Vs UK10K
if(o_population == "UK10K"){
  pdf(paste("shared_variants_",o_population,"_",population,".pdf",sep=""), width=8, height=4, pointsize=12)
    par(mfrow=c(1,2), mai=c(1,1,0.1,0.1),cex=0.8)
    slab <- seq(0,round(sxmax/1.5),40)
    # slab <- c(0,5,10,15,20,25)
    labnames=c("MAF<=0.05%","1%  ","2% ","5% ","10%  ","15%  ","20%  ","25%  ","30%  ","MAF>30% ")
    #left-hand plot

    m <- barplot(-rbind(current_overlap$overlap_uk10k, current_iso_un)/fscale,horiz=T,
                col=mycols[2:1],
                xlab=paste("number of variants (x",fscale,")"),
                axes=F,
                xlim=c(-sxmax+1,0),
                border=mycols[2:1],
                )

    axis(1, labels=slab, at=-slab)
    if (population == "FVG"){
      text((-current_overlap$tot_FVG/fscale)-60, m, p_uk10k_iso) #-> FVG plot
    }else if(population == "VBI"){
      text((-current_overlap$tot_VBI/fscale)-60, m, p_uk10k_iso) #-> VBI plot
    }
    legend(x="left",y="center", legend=c(paste("unique ",population,"-cohorts",sep=""), "shared",paste("unique",o_population)), fill=mycols, cex=0.8, border=mycols, box.col="white")
    # mtext("D", side=3, line=0, at=xmax, cex=1.5)
    #right-hand plot
    barplot(rbind(current_overlap$overlap_UK10K_VBI,current_out_un)/(fscale), horiz=T,
            names.arg=labnames, col=mycols[2:3], las=1,
            xlab=paste("number of variants (x",fscale,")"),
            axes=F, xlim=c(0,sxmax+1), border=mycols[2:3])
    axis(1, labels=slab, at=slab)
    if (population == "FVG"){
      text((current_overlap$total_UK10K/fscale)+60, m, p_iso_uk10k) #-> FVG plot
    }else if(population == "VBI"){
      text((current_overlap$total_UK10K/fscale)+60, m, p_iso_uk10k) #-> VBI plot
    }
    # legend(14,2, legend=c("shared","unique"), fill=c("black","red"), cex=0.8)

  dev.off()
}
######################################################################
#Isolates vs 1000G
if(o_population == "TGP"){
  pdf(paste("shared_variants_",o_population,"_",population,".pdf",sep=""), width=8, height=4, pointsize=12)
    par(mfrow=c(1,2), mai=c(1,1,0.1,0.1),cex=0.8)
    slab <- seq(0,round(sxmax/1.5),40)
    # slab <- c(0,5,10,15,20,25)
    labnames=c("MAF<=0.05%","1%  ","2% ","5% ","10%  ","15%  ","20%  ","25%  ","30%  ","MAF>30% ")
    #left-hand plot

    m <- barplot(-rbind(current_overlap$overlap_tgp, current_iso_un)/fscale,horiz=T,
                col=mycols[2:1],
                xlab=paste("number of variants (x",fscale,")"),
                axes=F,
                xlim=c(-sxmax+1,0),
                border=mycols[2:1],
                )

    axis(1, labels=slab, at=-slab)
    if (population == "FVG"){
      text((-current_overlap$tot_FVG/fscale)-60, m, p_tgp_iso) #-> FVG plot
    }else if(population == "VBI"){
      text((-current_overlap$tot_VBI/fscale)-60, m, p_tgp_iso) #-> VBI plot
    }
    legend(x="left",y="center", legend=c(paste("unique ",population,"-cohorts",sep=""), "shared",paste("unique",o_population)), fill=mycols, cex=0.8, border=mycols, box.col="white")
    # mtext("D", side=3, line=0, at=xmax, cex=1.5)
    #right-hand plot
    barplot(rbind(current_overlap$overlap_TGP_VBI,current_out_un)/(fscale), horiz=T,
            names.arg=labnames, col=mycols[2:3], las=1,
            xlab=paste("number of variants (x",fscale,")"),
            axes=F, xlim=c(0,sxmax+1), border=mycols[2:3])
    axis(1, labels=slab, at=slab)
    if (population == "FVG"){
      text((current_overlap$total_TGP/fscale)+60, m, p_iso_tgp) #-> FVG plot
    }else if(population == "VBI"){
      text((current_overlap$total_TGP/fscale)+60, m, p_iso_tgp) #-> VBI plot
    }
    # legend(14,2, legend=c("shared","unique"), fill=c("black","red"), cex=0.8)

  dev.off()
}




barplot(-rbind(uk10k4$overlap, unl4)/wfscale, horiz=T)

# in the plot below, then this makes the left-hand plot easier. 
# When I talked to V one day, it made me realise that I had made it too complicated.


mycols <- c("lightskyblue","steelblue3","steelblue4")
xlim <- max(c(uk10k4$total, g1k4$total, uk10k5$total, g1k5$total))/wfscale+3.5
slab <- c(0,5,10,15,20)
labnames=c("AC=1    ","AC=2    ",expression(paste("AF", ""<=5, "%  ")),expression(paste("AF", ""<=10, "% ")),"AF>10%  ")
m <- barplot(-rbind(uk10k4$overlap, unl4)/wfscale, horiz=T, col=mycols[2:1], las=1, 
             axes=F, border=mycols[2:1], xlim=c(-xlim,0), cex.names=1.5,
             xlab="number of variants (x1M)", cex.lab=1.8)
text(-uk10k4$total/wfscale-1.45, m, pl4, cex=1.5)
axis(1, labels=slab, at=-slab, cex.axis=1.3)
legend(x="topleft", legend=c("unique UK10K-cohorts", "shared","unique 1000GP"), fill=mycols, cex=1.3, border=mycols, box.col="white")
mtext("D", side=3, line=0, at=-xlim, cex=1.5)
barplot(rbind(g1k4$overlap, unr4)/wfscale, horiz=T, names.arg=labnames, col=mycols[2:3], las=1,
        axes=F,  border=mycols[2:3], xlim=c(0,xlim), cex.names=1.5, xlab="number of variants (x1M)", cex.lab=1.8)
text(g1k4$total/wfscale+1.45, m, pr4, cex=1.5)
axis(1, labels=slab, at=slab, cex.axis=1.3)