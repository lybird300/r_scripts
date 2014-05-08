#r script to plot for enza
rm(list=ls())
########################################################################
#plot 1 DAF spectrun for each population
#upload data for each population:

chr <- "22"
pops <- c("CEU","TSI","VBI","FVG")

all_pop_DAF <- NULL

for (pop in pops) {

pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")

pop_table <- read.table(pop_table_file,header=F)
colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
pop_table$DAF <- as.numeric(as.character(pop_table$DAF))

#remove nas
pop_table <- pop_table[-which(is.na(pop_table$DAF)),]
current_hist <- hist(pop_table$DAF,plot=F)

# current_pop_DAF <- cbind(as.data.frame(current_hist$breaks),as.data.frame(current_hist$counts),as.data.frame(rep(pop,length(current_hist$counts))))
# current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=(current_hist$counts/length(pop_table$POS)),pop=rep(pop,length(current_hist$counts))))
#use relative frequency sites per bin count/total variants number
current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pop,length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS))))

all_pop_DAF <- rbind(all_pop_DAF,current_pop_DAF)
}

all_pop_DAF$breaks <- as.numeric(as.character(all_pop_DAF$breaks))
all_pop_DAF$counts <- as.numeric(as.character(all_pop_DAF$counts))
# all_pop_DAF$counts <- all_pop_DAF$counts*length(pop_table$POS)
all_pop_DAF$rel_count <- as.numeric(as.character(all_pop_DAF$rel_count))
all_pop_DAF$pop <- as.character(all_pop_DAF$pop)

all_barplot_DAF <- NULL

for(i in 1:length(all_pop_DAF$pop)){
  current_row <- cbind(rep(all_pop_DAF[i,1],all_pop_DAF[i,2]),rep(as.character(all_pop_DAF[i,3]),all_pop_DAF[i,2]))
  
  all_barplot_DAF <- rbind(all_barplot_DAF,current_row)  
}

all_cols <- NULL

for(i in 1:length(pops)){
pop_col <- cbind(colors()[62+(i*10)],pops[i])
all_cols <- rbind(all_cols,pop_col)
}

jpeg("all_pop_DAF.jpg",width=1600, height=800,pointsize = 20)
  barplot(table(all_barplot_DAF[,2],all_barplot_DAF[,1]),beside=T,col=all_cols[,1], legend=all_cols[,2], main="DAF in all populations", xlab="DAF",ylab="Frequency")
dev.off()

#now upload the private and shared plots and add them to the previous plot
merged_daf <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/INGI_chr22.merged_daf.tab",header=F)
colnames(merged_daf) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","FVG")

merged_daf$CEU <- as.numeric(as.character(merged_daf$CEU))
merged_daf$TSI <- as.numeric(as.character(merged_daf$TSI))
merged_daf$VBI <- as.numeric(as.character(merged_daf$VBI))
merged_daf$FVG <- as.numeric(as.character(merged_daf$FVG))


pop_VBI_private <- read.table(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/VBI_private_chr",chr,".merged_daf.tab",sep=""),header=F)
pop_FVG_private <- read.table(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/FVG_private_chr",chr,".merged_daf.tab",sep=""),header=F)
pop_VBI_shared <- read.table(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/VBI_shared_chr",chr,".merged_daf.tab",sep=""),header = F)
pop_FVG_shared <- read.table(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/FVG_shared_chr",chr,".merged_daf.tab",sep=""),header = F)

colnames(pop_VBI_private) <- c("CHR","POZ","POS","VT","CEU","TSI","DAF","FVG")
pop_VBI_private$CEU <- as.numeric(as.character(pop_VBI_private$CEU))
pop_VBI_private$TSI <- as.numeric(as.character(pop_VBI_private$TSI))
pop_VBI_private$DAF <- as.numeric(as.character(pop_VBI_private$DAF))
pop_VBI_private$FVG <- as.numeric(as.character(pop_VBI_private$FVG))


colnames(pop_FVG_private) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","DAF")
pop_FVG_private$CEU <- as.numeric(as.character(pop_FVG_private$CEU))
pop_FVG_private$TSI <- as.numeric(as.character(pop_FVG_private$TSI))
pop_FVG_private$VBI <- as.numeric(as.character(pop_FVG_private$VBI))
pop_FVG_private$DAF <- as.numeric(as.character(pop_FVG_private$DAF))


colnames(pop_VBI_shared) <- c("CHR","POZ","POS","VT","CEU","TSI","DAF","FVG")
pop_VBI_shared$CEU <- as.numeric(as.character(pop_VBI_shared$CEU))
pop_VBI_shared$TSI <- as.numeric(as.character(pop_VBI_shared$TSI))
pop_VBI_shared$DAF <- as.numeric(as.character(pop_VBI_shared$DAF))
pop_VBI_shared$FVG <- as.numeric(as.character(pop_VBI_shared$FVG))


colnames(pop_FVG_shared) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","DAF")
pop_FVG_shared$CEU <- as.numeric(as.character(pop_FVG_shared$CEU))
pop_FVG_shared$TSI <- as.numeric(as.character(pop_FVG_shared$TSI))
pop_FVG_shared$VBI <- as.numeric(as.character(pop_FVG_shared$VBI))
pop_FVG_shared$DAF <- as.numeric(as.character(pop_FVG_shared$DAF))

pops_ingi <- c("VBI_p","FVG_p","VBI_s","FVG_s")
pops_ingi_files <- c("pop_VBI_private","pop_FVG_private","pop_VBI_shared","pop_FVG_shared")

for (k in 1:length(pops_ingi_files)){
current_pop <- get(pops_ingi_files[k])

current_hist <- hist(current_pop$DAF,plot=F)

#added relative frequency sites per bin count/total variants number
current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pops_ingi[k],length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS))))
current_pop_DAF$breaks <- as.numeric(as.character(current_pop_DAF$breaks))
current_pop_DAF$counts <- as.numeric(as.character(current_pop_DAF$counts))
current_pop_DAF$rel_count <- as.numeric(as.character(current_pop_DAF$rel_count))
current_pop_DAF$pop <- as.character(current_pop_DAF$pop)

#add again to all_pop_DAF to have a complete plot
all_pop_DAF <- rbind(all_pop_DAF,current_pop_DAF)
  
}

#create the object to plot
all_barplot_DAF_sp <- NULL

for(i in 1:length(all_pop_DAF$pop)){
  current_row <- cbind(rep(all_pop_DAF[i,1],all_pop_DAF[i,2]),rep(as.character(all_pop_DAF[i,3]),all_pop_DAF[i,2]))
  
  all_barplot_DAF_sp <- rbind(all_barplot_DAF_sp,current_row)  
}

all_cols <- NULL
pops2 <- c(pops,pops_ingi)

for(i in 1:length(pops2)){
pop_col <- cbind(colors()[62+(i*10)],pops2[i])
all_cols <- rbind(all_cols,pop_col)
}

jpeg("all_pop_DAF_sp.jpg",width=1600, height=800,pointsize = 20)
  barplot(table(all_barplot_DAF_sp[,2],all_barplot_DAF_sp[,1]),beside=T,col=all_cols[,1], legend=all_cols[,2], main="DAF in all populations", xlab="DAF",ylab="Frequency")
dev.off()

###########################################################################################
#Plot 2: venn diagram with overlap between all populations and categories

# require(venneuler)

# #we need to do all the cases for intersections
# for(l in 1:length(pops)){
#   pop_index <- grep(pops[l],colnames(merged_daf))

#   current_gt0 <- merged_daf[which(merged_daf[,pop_index] > 0),]
#   current_shared_all <- merged_daf[which(merged_daf[,pop_index] > 0 & merged_daf[,pop_index+1] > 0 & merged_daf[,pop_index+2] > 0 & merged_daf[,pop_index+3] > 0),]
#   current_private <- merged_daf[which(merged_daf[,pop_index] > 0 & merged_daf[,pop_index+1] == 0 & merged_daf[,pop_index+2] == 0 & merged_daf[,pop_index+3] == 0),]
# }
# # v <- venneuler(c(FVG=length(merged_daf[which(merged_daf$FVG > 0),]$CHR),
# #   VBI=length(merged_daf[which(merged_daf$VBI > 0),]$CHR),
# #   "VBI&FVG"=length(merged_daf[which(merged_daf$FVG > 0 & merged_daf$VBI > 0),]$CHR)))

# # v <- draw.quad.venn(length(merged_daf[which(merged_daf$CEU > 0),]$CHR),
# v <- draw.pairwise.venn(length(merged_daf[which(merged_daf$CEU > 0),]$CHR),
#   length(merged_daf[which(merged_daf$TSI > 0),]$CHR),
#   length(merged_daf[which(merged_daf$TSI > 0 & merged_daf$CEU > 0),]$CHR))
#   # length(merged_daf[which(merged_daf$VBI > 0),]$CHR),
#   # length(merged_daf[which(merged_daf$FVG > 0),]$CHR),  
#   # n12=length(merged_daf[which(merged_daf$FVG > 0 & merged_daf$VBI > 0 & merged_daf$TSI > 0 & merged_daf$CEU > 0),]$CHR)),
#   # n13=length(merged_daf[which(merged_daf$FVG > 0 & merged_daf$VBI > 0 & merged_daf$TSI > 0 & merged_daf$CEU > 0),]$CHR)),
#   # n14=length(merged_daf[which(merged_daf$FVG > 0 & merged_daf$VBI > 0 & merged_daf$TSI > 0 & merged_daf$CEU > 0),]$CHR)),
#   # n1234=length(merged_daf[which(merged_daf$FVG > 0 & merged_daf$VBI > 0 & merged_daf$TSI > 0 & merged_daf$CEU > 0),]$CHR)),

# # jpeg(paste("DAF.jpg",sep="_"),width=800, height=800,pointsize = 20)
# jpeg("VENN_DAF.jpg",width=800, height=800,pointsize = 20)
#   grid.draw(v)
# dev.off()


###########################################################################################
#Plot 3 barplot for privte and shared variants for each INGI population by category
#use also here the relative frequency

funk_cat <-read.csv("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/consequences.list",header=F)
funk_cat$V1 <- as.character(funk_cat$V1)

chr <- "22"
pop_is <- c("VBI","FVG")


for (pop in pop_is) {

  conseq_count_tot <- NULL

  pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")

  pop_table <- read.table(pop_table_file,header=T)
  # colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
  current_private <- get(paste("pop",pop,"private",sep="_"))
  current_shared <- get(paste("pop",pop,"shared",sep="_"))

  for(cat in funk_cat$V1){

    current_grep <- pop_table[(grep(cat,pop_table$INFO)),]

    #merge the grepped data with the private
    current_grep_private <- merge(current_private,current_grep,by.x="POS",by.y="POS")

    current_grep_shared <- merge(current_shared,current_grep,by.x="POS",by.y="POS")

    current_conseq <- cbind(category=as.character(cat),private=length(current_grep_private$POS),shared=length(current_grep_shared$POS))
    conseq_count_tot <- rbind(conseq_count_tot,current_conseq)

  }

  conseq_count_tot <- as.data.frame(conseq_count_tot)
  conseq_count_tot$category <- as.character(conseq_count_tot$category)
  conseq_count_tot$private <- as.numeric(as.character(conseq_count_tot$private))
  conseq_count_tot$shared <- as.numeric(as.character(conseq_count_tot$shared))
  conseq_count_tot$N_shared <- rep(length(current_shared$POS),length(conseq_count_tot$category))
  conseq_count_tot$N_private <- rep(length(current_private$POS),length(conseq_count_tot$category))
  # conseq_count_tot$N_shared <- as.numeric(as.character(conseq_count_tot$N_shared))
  # conseq_count_tot$N_private <- as.numeric(as.character(conseq_count_tot$N_private))

  #lets check if the differences are significative with a chisquare test
  for(i in 1:length(conseq_count_tot$category)){
    print(conseq_count_tot$category[i])
    A=matrix(c(conseq_count_tot$shared[i], conseq_count_tot$N_shared[i], conseq_count_tot$private[i], conseq_count_tot$N_private[i]), nrow=2, ncol=2, byrow=T)
    pippo=chisq.test(A)
    conseq_count_tot$p[i] <- pippo$p.value
  }

  write.table(conseq_count_tot,file=paste(pop,"_consequences_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

  all_barplot_conseq <- NULL

  for(i in 1:length(conseq_count_tot$category)){
    for(j in 2:length(colnames(conseq_count_tot))){
      current_row <- cbind(rep(as.character(conseq_count_tot[i,1]),conseq_count_tot[i,j]),rep(colnames(conseq_count_tot)[j],conseq_count_tot[i,j]))
      all_barplot_conseq <- rbind(all_barplot_conseq,current_row)  
    }
  }

  #we need to use the same expedient used for the maf count
  jpeg(paste(pop,"all_conseq_DAF.jpg",sep="_"),width=1800, height=800,pointsize = 20)
    barplot(table(all_barplot_conseq[,2],all_barplot_conseq[,1]),beside=T,col=c("red","blue"), legend=c("Private","Shared"), main="DAF in private/shared sites for different consequences annotations", xlab="Annotations",ylab="Relative frequency", las=2 )
  dev.off()
}

###########################################################################################
#Plot 4 For each population (VBI FVG)  density distribution + wilcoxon test  (shared e private)  GERP SCORE  [stratify for  genic/intergenic? functional categories?]

###########################################################################################
#Plot 5 SIFT Polyphen  https://faculty.washington.edu/wjs18/GS561/cSNPs_lab.html     [stratify for  genic/intergenic? functional categories?]

###########################################################################################


#  <- merged_daf[which(merged_daf$VBI > 0 & merged_daf$TSI == 0 & merged_daf$CEU ==0),]
#  <- merged_daf[which(merged_daf$FVG > 0 & merged_daf$TSI == 0 & merged_daf$CEU ==0),]


# merged_daf[which(merged_daf$VBI > 0 & (merged_daf$TSI > 0 | merged_daf$CEU > 0)),]
# merged_daf[which(merged_daf$FVG > 0 & (merged_daf$TSI > 0 | merged_daf$CEU > 0)),]

#we need to extract all cases of sharing in order to do the venn diagram
#CEU
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'
# #TSI
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'
# #VBI
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'
# #FVG
# awk '$6 != "na" && $5 != "na" && $7 != "na" && $7 > 0 && $6 >0 && $5 > 0'