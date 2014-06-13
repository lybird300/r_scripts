#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
########################################################################
#plot 1 DAF spectrun for each population
#upload data for each population:

chr <- "22"
pops <- c("CEU","TSI","VBI","FVG")

all_pop_DAF <- NULL

for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  pop_table_file <- paste(pop,".chr",chr,".tab.gz",sep="")

  pop_table <- read.table(pop_table_file,header=F,skip=1,stringsAsFactors=F, comment.char="")
  colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  pop_table$DAF <- as.numeric(as.character(pop_table$DAF))

  #remove nas
  pop_table <- pop_table[-which(is.na(pop_table$DAF)),]
  #remove fixed variants (the ones with DAF == 0 or == 1)
  pop_table <- pop_table[which(pop_table$DAF != 0),]
  pop_table <- pop_table[which(pop_table$DAF != 1),]

  current_hist <- hist(pop_table$DAF,plot=F)

  # current_pop_DAF <- cbind(as.data.frame(current_hist$breaks),as.data.frame(current_hist$counts),as.data.frame(rep(pop,length(current_hist$counts))))
  # current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=(current_hist$counts/length(pop_table$POS)),pop=rep(pop,length(current_hist$counts))))
  #use relative frequency sites per bin count/total variants number
  current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pop,length(current_hist$counts)),rel_count=(current_hist$counts/length(pop_table$POS))))
  
  all_pop_DAF <- rbind(all_pop_DAF,current_pop_DAF)

}

all_pop_DAF$breaks <- as.numeric(as.character(all_pop_DAF$breaks))
all_pop_DAF$counts <- as.numeric(as.character(all_pop_DAF$counts))
all_pop_DAF$rel_count <- as.numeric(as.character(all_pop_DAF$rel_count))
all_pop_DAF$pop <- as.character(all_pop_DAF$pop)

all_pop_DAF_table <- all_pop_DAF[which(all_pop_DAF$pop == "CEU"),1]
all_pop_DAF_table_names <- c("breaks")

for(pop_i in pops){
  actual_pop <- all_pop_DAF[which(all_pop_DAF$pop == pop_i),]
  all_pop_DAF_table_names <- c(all_pop_DAF_table_names,paste(pop_i,"count",sep="_"),paste(pop_i,"rel_count",sep="_"))

  all_pop_DAF_table <- cbind(all_pop_DAF_table,cbind(actual_pop$counts,actual_pop$rel_count*100))
}
all_pop_DAF_table <- as.data.frame(all_pop_DAF_table)
colnames(all_pop_DAF_table) <- all_pop_DAF_table_names

write.table(all_pop_DAF_table,file=paste("All_pop_DAF_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#absolute count
barplot3 <- as.matrix(t(all_pop_DAF_table[,c(2,4,6,8)]))

#relative site count
barplot4 <- as.matrix(t(all_pop_DAF_table[,c(3,5,7,9)]))

all_cols <- NULL

for(i in 1:length(pops)){
  if(pops[i] == "CEU"){
    cur_col <- colors()[41]
  }
  if(pops[i] == "FVG"){
    cur_col <- colors()[52]
  }
  if(pops[i] == "FVG_p"){
    cur_col <- colors()[56]
  }
  if(pops[i] == "FVG_s"){
    cur_col <- colors()[61]
  }
  if(pops[i] == "TSI"){
    cur_col <- colors()[72]
  }
  if(pops[i] == "VBI"){
    cur_col <- colors()[81]
  }
  if(pops[i] == "VBI_p"){
    cur_col <- colors()[85]
  }
  if(pops[i] == "VBI_s"){
    cur_col <- colors()[89]
  }
  
  pop_col <- cbind(cur_col,pops[i])
  all_cols <- rbind(all_cols,pop_col)
}

jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot3,
    beside=T,
    main="DAF in all populations",
    names.arg=all_pop_DAF_table$breaks,
    xlab="",
    ylab="",col=all_cols[,1])
    legend("topright",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=length(all_cols[,1]))
    mtext(1, text = "Annotations", line = 4,cex=1.4)
    mtext(2, text = "Site count", line = 4,cex=1.4)
dev.off()

jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_rf.jpg",sep=""),width=1800, height=800)
  par(oma=c(3,3,3,3),cex=1.4)
  barplot(barplot4,
    beside=T,
    main="DAF in all populations",
    names.arg=all_pop_DAF_table$breaks,
    xlab="",
    ylab="",col=all_cols[,1])
    legend("topright",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=length(all_cols[,1]))
    mtext(1, text = "Annotations", line = 4,cex=1.4)
    mtext(2, text = "Relative Frequency (N sites/Tot sites with DAF information )(%)", line = 4,cex=1.4)
dev.off()


#now upload the private and shared plots and add them to the previous plot
merged_daf <- read.table(paste0(base_folder,"INPUT_FILES/INGI_chr",chr,".merged_daf.tab.gz"),skip=1,header=F)
colnames(merged_daf) <- c("CHR","POZ","POS","VT","CEU","TSI","VBI","FVG")

merged_daf$CEU <- as.numeric(as.character(merged_daf$CEU))
merged_daf$TSI <- as.numeric(as.character(merged_daf$TSI))
merged_daf$VBI <- as.numeric(as.character(merged_daf$VBI))
merged_daf$FVG <- as.numeric(as.character(merged_daf$FVG))

pop_VBI_private <- read.table(paste(base_folder,"INPUT_FILES/VBI_private_chr",chr,".merged_daf.tab.gz",sep=""),header=F)
pop_FVG_private <- read.table(paste(base_folder,"INPUT_FILES/FVG_private_chr",chr,".merged_daf.tab.gz",sep=""),header=F)
pop_VBI_shared <- read.table(paste(base_folder,"INPUT_FILES/VBI_shared_chr",chr,".merged_daf.tab.gz",sep=""),header = F)
pop_FVG_shared <- read.table(paste(base_folder,"INPUT_FILES/FVG_shared_chr",chr,".merged_daf.tab.gz",sep=""),header = F)

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

summary(pop_VBI_private)
summary(pop_FVG_private)
summary(pop_VBI_shared)
summary(pop_FVG_shared)

dim(pop_VBI_private)
dim(pop_FVG_private)
dim(pop_VBI_shared)
dim(pop_FVG_shared)

#this is useless now...
# pop_VBI_private <- pop_VBI_private[which(pop_VBI_private$DAF != 1),]
# pop_FVG_private <- pop_FVG_private[which(pop_FVG_private$DAF != 1),]
# pop_VBI_shared <- pop_VBI_shared[which(pop_VBI_shared$DAF != 1),]
# pop_FVG_shared <- pop_FVG_shared[which(pop_FVG_shared$DAF != 1),]

pops_ingi <- c("VBI_p","FVG_p","VBI_s","FVG_s")
pops_ingi_files <- c("pop_VBI_private","pop_FVG_private","pop_VBI_shared","pop_FVG_shared")

for (k in 1:length(pops_ingi_files)){
  current_pop <- get(pops_ingi_files[k])

  current_hist <- hist(current_pop$DAF,plot=F)

  #added relative frequency sites per bin count/total variants number for the current category! shared or private!!
  current_pop_DAF <- as.data.frame(cbind(breaks=current_hist$breaks[2:length(current_hist$breaks)],counts=current_hist$counts,pop=rep(pops_ingi[k],length(current_hist$counts)),rel_count=(current_hist$counts/length(current_pop$POS))))
  current_pop_DAF$breaks <- as.numeric(as.character(current_pop_DAF$breaks))
  current_pop_DAF$counts <- as.numeric(as.character(current_pop_DAF$counts))
  current_pop_DAF$rel_count <- as.numeric(as.character(current_pop_DAF$rel_count))
  current_pop_DAF$pop <- as.character(current_pop_DAF$pop)

  #add again to all_pop_DAF to have a complete plot
  all_pop_DAF <- rbind(all_pop_DAF,current_pop_DAF)
  
}

all_pop_DAF_table_sp <- all_pop_DAF[which(all_pop_DAF$pop == "CEU"),1]
all_pop_DAF_table_sp_names <- c("breaks")

pops2 <- sort(c(pops,pops_ingi))

for(pop_i2 in pops2){
  actual_pop2 <- all_pop_DAF[which(all_pop_DAF$pop == pop_i2),]
  all_pop_DAF_table_sp_names <- c(all_pop_DAF_table_sp_names,paste(pop_i2,"count",sep="_"),paste(pop_i2,"rel_count",sep="_"))

  all_pop_DAF_table_sp <- cbind(all_pop_DAF_table_sp,cbind(actual_pop2$counts,actual_pop2$rel_count*100))
}

all_pop_DAF_table_sp <- as.data.frame(all_pop_DAF_table_sp)
colnames(all_pop_DAF_table_sp) <- all_pop_DAF_table_sp_names

write.table(all_pop_DAF_table_sp,file=paste("All_pop_DAF_sp_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#absolute count
barplot5 <- as.matrix(t(all_pop_DAF_table_sp[,c(2,4,6,8,10,12,14,16)]))

#relative site count
barplot6 <- as.matrix(t(all_pop_DAF_table_sp[,c(3,5,7,9,11,13,15,17)]))

#relative count but only for private sites
barplot7 <- as.matrix(t(all_pop_DAF_table_sp[,c(3,5,9,11,13,15)]))

all_cols <- NULL
# pops2 <- sort(c(pops,pops_ingi))

for(i in 1:length(pops2)){
  if(pops2[i] == "CEU"){
    cur_col <- colors()[130]
  }
  if(pops2[i] == "FVG"){
    cur_col <- colors()[517]
  }
  if(pops2[i] == "FVG_p"){
    cur_col <- colors()[50]
  }
  if(pops2[i] == "FVG_s"){
    cur_col <- colors()[81]
  }
  if(pops2[i] == "TSI"){
    cur_col <- colors()[30]
  }
  if(pops2[i] == "VBI"){
    cur_col <- colors()[421]
  }
  if(pops2[i] == "VBI_p"){
    cur_col <- colors()[33]
  }
  if(pops2[i] == "VBI_s"){
    cur_col <- colors()[36]
  }
  
  pop_col <- cbind(cur_col,pops2[i])
  all_cols <- rbind(all_cols,pop_col)
}

jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_sp.jpg",sep=""),width=1800, height=800)
    par(oma=c(3,3,3,3),cex=1.4)
    barplot(barplot5,
      beside=T,
      main="DAF in all populations",
      names.arg=all_pop_DAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[,1]
      )
      mtext(1, text = "DAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Site count", line = 4,cex=1.4)
      legend(x="top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=8)
  dev.off()

  jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_sp_rf.jpg",sep=""),width=1800, height=800)
      par(oma=c(3,3,3,3),cex=1.4)
    barplot(barplot6,
      beside=T,
      main="DAF in all populations",
      names.arg=all_pop_DAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[,1],ylim=c(0,60))
      mtext(1, text = "DAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Relative Frequency (N sites/Tot sites per category)(%)", line = 4,cex=1.4)
      legend("top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[,1],legend=all_cols[,2], ncol=8)
  dev.off()

  #plot only private, not shared
  jpeg(paste(base_folder,"PLOTS/1_all_pop_DAF_p_rf.jpg",sep=""),width=1800, height=800)
      par(oma=c(3,3,3,3),cex=1.4)
    barplot(barplot7,
      beside=T,
      main="DAF in all populations",
      names.arg=all_pop_DAF_table_sp$breaks*100,
      xlab="",
      ylab="",col=all_cols[c(1,2,3,5,6,7),1],ylim=c(0,60))
      mtext(1, text = "DAF (%)", line = 4,cex=1.4)
      mtext(2, text = "Relative Frequency (N sites/Tot sites per category)(%)", line = 4,cex=1.4)
      legend("top",pch =c(rep(19,length(all_cols[,1]))),col=all_cols[c(1,2,3,5,6,7),1],legend=all_cols[c(1,2,3,5,6,7),2], ncol=6)
  dev.off()

###########################################################################################
#Plot 2: venn diagram with overlap between all populations and categories

require(gplots)
require(VennDiagram)

#We need to get everything that has DAF > 0 for each population, this means that in each row we want to keep stuff if there is at least a DAF > 0
#But since we're trying to show how ALL variants split in the 4 populations, we NEED TO DO NOTHING HERE!!!WE NEED ALL VARIANTS!!!

# #we are going to remove first the fixed sites from isolate populations
# merged_daf_diag_not_fixed_iso <- merged_daf[which(merged_daf$FVG != 1 & merged_daf$VBI != 1 ),]
# #now we'll remove from this dataset all the fixed sites with DAF == 0
# merged_daf_diag_not_fixed_iso <- merged_daf_diag_not_fixed_iso[which(merged_daf_diag_not_fixed_iso$FVG != 0 | merged_daf_diag_not_fixed_iso$VBI != 0 ),]

# merged_daf_diag_not_fixed <- merged_daf[which(merged_daf$CEU != 1 & merged_daf$FVG != 1 & merged_daf$TSI != 1 & merged_daf$VBI != 1 ),]
# merged_daf_diag_not_fixed <- merged_daf_diag_not_fixed[which(merged_daf_diag_not_fixed$CEU != 0 | merged_daf_diag_not_fixed$FVG != 0 | merged_daf_diag_not_fixed$TSI != 0 | merged_daf_diag_not_fixed$VBI != 0 ),]

# merged_daf_diag <- merged_daf_diag_not_fixed[,-c(1:4)]

# merged_daf_diag <- merged_daf_diag_not_fixed_iso[,-c(1:4)]
merged_daf_diag <- merged_daf[,-c(1:4)]

merged_daf_diag[merged_daf_diag>0]=1
aree=colSums(merged_daf_diag,na.rm=T)

# #easy venn diagramm EVER
# res.all <- list(which(merged_daf_diag[,1]>0),
#   which(merged_daf_diag[,2]>0),
#   which(merged_daf_diag[,3]>0),
#   which(merged_daf_diag[,4]>0))
# names(res.all) <- colnames(merged_daf_diag)

#we need to do all the cases for intersections
res.all2 <- NULL
for(i in c(0,1)){
  for(j in c(0,1)){
    for(k in c(0,1)){
      for(h in c(0,1)){
        if(any(c(i,j,k,h)==1)){
          res=length(which(merged_daf_diag$CEU==i & merged_daf_diag$TSI==j & merged_daf_diag$VBI==k & merged_daf_diag$FVG==h))
          res.all2=rbind(res.all2,c(res,paste(which(c(i,j,k,h)==1),collapse="")))
        }
      }
    }
  }
}

res.all3 <- as.data.frame(res.all2[order(res.all2[,2]),])
res.all3$V1 <- as.numeric(as.character(res.all3$V1))
res.all3$V2 <- as.character(res.all3$V2)
res.all3$V3 <- as.character(res.all3$V2)
colnames(res.all3) <- c("n_sites","cat_code","pop")

for(i in 1:length(res.all3$cat_code)){
  for(j in names(aree)){
    cat <- res.all3$pop[i]
    newcat <- sub(as.character(which(names(aree)==j)),paste(j,";",sep=""),cat)
    res.all3$pop[i] <- newcat
  }
}

all_cols <- NULL
cats <- colnames(merged_daf_diag)

for(i in 1:length(cats)){
  if(cats[i] == "CEU"){
    cur_col <- colors()[130]
  }
  if(cats[i] == "FVG"){
    cur_col <- colors()[517]
  }
  if(cats[i] == "FVG_p"){
    cur_col <- colors()[50]
  }
  if(cats[i] == "FVG_s"){
    cur_col <- colors()[81]
  }
  if(cats[i] == "TSI"){
    cur_col <- colors()[30]
  }
  if(cats[i] == "VBI"){
    cur_col <- colors()[421]
  }
  if(cats[i] == "VBI_p"){
    cur_col <- colors()[33]
  }
  if(cats[i] == "VBI_s"){
    cur_col <- colors()[36]
  }
  
  pop_col <- cbind(cur_col,cats[i])
  all_cols <- rbind(all_cols,pop_col)
}

#write the table:
write.table(res.all3,file=paste(base_folder,"RESULTS/VENN/All_pop_intersect_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#find all indexes
idx1=grep("1",res.all2[,2])
idx2=grep("2",res.all2[,2])
idx3=grep("3",res.all2[,2])
idx4=grep("4",res.all2[,2])

v <- draw.quad.venn(area1=as.numeric(sum(as.numeric(res.all2[idx1,1]))),
  area2=as.numeric(sum(as.numeric(res.all2[idx2,1]))),
  area3=as.numeric(sum(as.numeric(res.all2[idx3,1]))),
  area4=as.numeric(sum(as.numeric(res.all2[idx4,1]))),
  n12=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx2),1]))),
  n13=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx3),1]))),
  n14=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx4),1]))),
  n23=as.numeric(sum(as.numeric(res.all2[intersect(idx2,idx3),1]))),
  n24=as.numeric(sum(as.numeric(res.all2[intersect(idx2,idx4),1]))),
  n34=as.numeric(sum(as.numeric(res.all2[intersect(idx3,idx4),1]))),
  n123=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx2),idx3),1]))),
  n124=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx2),idx4),1]))),
  n134=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx3),idx4),1]))),
  n234=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx2,idx3),idx4),1]))),
  n1234=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(idx1,idx2),idx3),idx4),1]))),
  ind=FALSE, category=colnames(merged_daf_diag),fill=all_cols[,1])

jpeg(paste(base_folder,"PLOTS/2_all_pop_DAF_VENN.jpg",sep=""),width=800, height=800,pointsize = 20)
  grid.draw(v)
dev.off()


###########################################################################################
#Plot 3 barplot for privte and shared variants for each INGI population by category
#use also here the relative frequency

funk_cat <-read.csv(paste(base_folder,"INPUT_FILES/consequences.list",sep=""),header=F)
funk_cat$V1 <- as.character(funk_cat$V1)

chr <- "22"
pop_is <- c("VBI","FVG")

for (pop in pop_is) {

  conseq_count_tot <- NULL

  pop_table_file <- paste(base_folder,"INPUT_FILES/",pop,".chr",chr,".tab.gz",sep="")

  pop_table <- read.table(pop_table_file,header=F,skip=1)
  colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
  current_private <- get(paste("pop",pop,"private",sep="_"))
  current_shared <- get(paste("pop",pop,"shared",sep="_"))

  for(cat in funk_cat$V1){

    current_grep <- pop_table[(grep(cat,pop_table$INFO)),]

    #merge the grepped data with the private
    current_grep_private <- merge(current_private,current_grep,by.x="POS",by.y="POS")

    #merge the grepped data with the shared
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
  conseq_count_tot$rel_private <- (conseq_count_tot$private/conseq_count_tot$N_private)*100
  conseq_count_tot$rel_shared <- (conseq_count_tot$shared/conseq_count_tot$N_shared)*100

  #lets check if the differences are significative with a chisquare test
  for(i in 1:length(conseq_count_tot$category)){
    print(conseq_count_tot$category[i])
    A=matrix(c(conseq_count_tot$shared[i], conseq_count_tot$N_shared[i], conseq_count_tot$private[i], conseq_count_tot$N_private[i]), nrow=2, ncol=2, byrow=T)
    pippo=chisq.test(A)
    conseq_count_tot$p[i] <- pippo$p.value
  }


  write.table(conseq_count_tot,file=paste(base_folder,"RESULTS/CONSEQUENCES/",pop,"_consequences_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)
  
  #site count
  barplot1 <- as.matrix(t(conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),c(2:3)]))

  #relative site count
  barplot2 <- as.matrix(t(conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),c(6:7)]))

  jpeg(paste(base_folder,"PLOTS/3_",pop,"_all_conseq_DAF.jpg",sep=""),width=1800, height=800)
    par(oma=c(9,3,3,3),cex=1.4)
    barplot(barplot1,
      beside=T,
      col=c("red","blue"),
      legend=c(paste("Private (tot:",conseq_count_tot$N_private[1],")",sep=""),paste("Shared (tot:",conseq_count_tot$N_shared[1],")",sep="")),
      main=paste(pop," DAF in private/shared sites categories for different consequences annotation classes",sep=""),
      names.arg=conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),1],
      xlab="",las=2,
      ylab="")
      mtext(1, text = "Annotations", line = 11)
      mtext(2, text = "Frequency", line = 4)
  dev.off()

  jpeg(paste(base_folder,"PLOTS/3_",pop,"_all_conseq_DAF_rf.jpg",sep=""),width=1800, height=800)
    par(oma=c(9,3,3,3),cex=1.4)
    barplot(barplot2,
      beside=T,
      col=c("red","blue"),
      legend=c("Private","Shared"),
      main=paste(pop," DAF in private/shared sites categories for different consequences annotation classes",sep=""),
      names.arg=conseq_count_tot[-which(conseq_count_tot$private == 0 & conseq_count_tot$shared == 0),1],
      xlab="",las=2,
      ylab="")
      mtext(1, text = "Annotations", line = 11)
      mtext(2, text = "Relative Frequency (N sites/Tot sites in category)(%)", line = 4)
  dev.off()
}

###########################################################################################
#Plot 4 For each population (VBI FVG)  density distribution + wilcoxon test  (shared e private)  GERP SCORE  [stratify for  genic/intergenic? functional categories?]
# We'll plot the gerp score distribution for shared and private
pop_is <- c("VBI","FVG")

for (pop in pop_is) {

  conseq_count_tot <- NULL

  pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/INPUT_FILES/",pop,".chr",chr,".tab",sep="")

  pop_table <- read.table(pop_table_file,header=T)
  # colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
  pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
  current_private <- get(paste("pop",pop,"private",sep="_"))
  current_shared <- get(paste("pop",pop,"shared",sep="_"))

  #we need to get the information of the gerp score
  current_grep_GERP <- pop_table[(grep("GERP",pop_table$INFO)),]
  
  #merge the grepped data with the private
  current_grep_GERP_private <- merge(current_private,current_grep_GERP,by.x="POS",by.y="POS")
  
  #merge the grepped data with the shared
  current_grep_GERP_shared <- merge(current_shared,current_grep_GERP,by.x="POS",by.y="POS")
  
  current_conseq <- cbind(category=as.character(cat),private=length(current_grep_private$POS),shared=length(current_grep_shared$POS))
  conseq_count_tot <- rbind(conseq_count_tot,current_conseq)

}

###########################################################################################
#Plot 5 SIFT Polyphen  https://faculty.washington.edu/wjs18/GS561/cSNPs_lab.html     [stratify for  genic/intergenic? functional categories?]

###########################################################################################


###########################################################################################
#Plot 6 ROH
chr <- "22"
pops <- c("CEU","TSI","VBI","FVG")

# input_format <- "PLINK"
input_format <- "BEAGLE"
xmax <- NULL
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    pop_roh_file <- paste(pop,"roh.length.hbd",sep=".")
    pop_roh_table <- read.table(pop_roh_file,header=FALSE,sep=" ")
    colnames(pop_roh_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_roh_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_roh_table <- read.table(pop_roh_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_roh_table <- pop_roh_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_roh_table$IID2 <- as.character(pop_roh_table$IID2)
  pop_roh_table$IID1 <- as.character(pop_roh_table$IID1)
  current_pop_table_name <- paste(pop,"roh",sep="_")
  assign(current_pop_table_name,pop_roh_table)

  tot_roh <- tapply(pop_roh_table$LENGTH,pop_roh_table$IID1,sum)
  tot_roh <- as.data.frame(tot_roh, row.names=NULL)
  tot_roh$ID <- rownames(tot_roh)
  colnames(tot_roh) <- c("ROH_tot","ID")
  tot_roh$ID <- as.character(tot_roh$ID)
  tot_roh$ROH_tot <- as.numeric(as.character(tot_roh$ROH_tot))

  assign(paste(pop,"tot_roh",sep="_"),tot_roh)
  
  assign(paste("M",pop,sep="_"),ecdf(tot_roh$ROH_tot))
  xmax <- c(xmax,summary(ecdf(tot_roh$ROH_tot))[6])
  # jpeg(paste("6_roh_",pop".jpg",sep=""),width=800, height=800)
  # plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity")
  # lines(M_FVG,FVG_tot_roh$ROH_tot,col="red")
  # lines(M_TSI,TSI_tot_roh$ROH_tot,col="green")
  # lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue")
  # legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
  # dev.off()

}

# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"PLOTS/6_roh_all.jpg",sep=""),width=800, height=800)
  par(lwd=2)
  plot(M_CEU,CEU_tot_roh$ROH_tot,main="",xlab="Total ROH homozigosity", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  lines(M_FVG,FVG_tot_roh$ROH_tot,col="red", verticals=TRUE, pch=46)
  lines(M_TSI,TSI_tot_roh$ROH_tot,col="green", verticals=TRUE, pch=46)
  lines(M_VBI,VBI_tot_roh$ROH_tot,col="blue", verticals=TRUE, pch=46)
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()

jpeg(paste(base_folder,"PLOTS/6_roh_iso.jpg",sep=""),width=800, height=800)
  plot(M_FVG,FVG_tot_roh$ROH_tot,verticals=TRUE, pch=46,main="",xlab="Total ROH homozigosity",col="red")
  lines(M_VBI,VBI_tot_roh$ROH_tot,verticals=TRUE, pch=46,col="blue")
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("FVG","VBI"),col=c("red","blue"),ncol=4)
dev.off()



# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/7_roh.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"PLOTS/7_roh.jpg",sep=""),width=800, height=800)
  plot(density(VBI_tot_roh$ROH_tot),main="",xlab="Total ROH homozigosity", col="blue")
  lines(density(FVG_tot_roh$ROH_tot),col="red")
  lines(density(TSI_tot_roh$ROH_tot),col="green")
  lines(density(CEU_tot_roh$ROH_tot))
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()

######################################################################################################
# Plot 7: IBD
###############################################################################
#Same as roh, but we consider samples's pairs:
# input_format <- "PLINK"
chr <- "22"
pops <- c("CEU","TSI","VBI","FVG")

input_format <- "BEAGLE"
xmax <- NULL
for (pop in pops) {

  # pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")
  if (input_format == "BEAGLE"){
    #this bit read BEAGLE output
    pop_ibd_file <- paste(pop,"roh.length.ibd",sep=".")
    pop_ibd_table <- read.table(pop_ibd_file,header=FALSE,sep=" ")
    colnames(pop_ibd_table) <- c("IID1","I1","IID2","I2","CHROM","START","END","LOD","PAIR","LENGTH")
    
  }else if(input_format == "PLINK"){
    #to read plink input we need this:
    pop_ibd_file <- paste(pop,"roh.hom.indiv",sep=".")
    pop_ibd_table <- read.table(pop_ibd_file,header=TRUE)
    colnames(pop_roh_table) <- c("IID1","IID2","PHE","NSEG","LENGTH","KBAVG")
    #remove samples with NSEG == 0
    pop_ibd_table <- pop_ibd_table[which(pop_roh_table$NSEG != 0),]
    
  }
  pop_ibd_table$IID2 <- as.character(pop_ibd_table$IID2)
  pop_ibd_table$IID1 <- as.character(pop_ibd_table$IID1)
  pop_ibd_table$PAIR <- as.character(pop_ibd_table$PAIR)
  current_pop_table_name <- paste(pop,"ibd",sep="_")
  assign(current_pop_table_name,pop_ibd_table)

  tot_ibd <- tapply(pop_ibd_table$LENGTH,pop_ibd_table$PAIR,sum)
  tot_ibd <- as.data.frame(tot_ibd, row.names=NULL)
  tot_ibd$ID <- rownames(tot_ibd)
  colnames(tot_ibd) <- c("IBD_tot","ID")
  tot_ibd$ID <- as.character(tot_ibd$ID)
  tot_ibd$IBD_tot <- as.numeric(as.character(tot_ibd$IBD_tot))

  assign(paste(pop,"tot_ibd",sep="_"),tot_ibd)
  
  assign(paste("M_ibd",pop,sep="_"),ecdf(tot_ibd$IBD_tot))
  xmax <- c(xmax,summary(ecdf(tot_ibd$IBD_tot))[6])

}

# jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PLOTS/6_roh.jpg",sep=""),width=800, height=800)
jpeg(paste(base_folder,"PLOTS/8_ibd.jpg",sep=""),width=800, height=800)
  par(lwd=2)
  plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
  lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col="red", verticals=TRUE, pch=46)
  lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col="green", verticals=TRUE, pch=46)
  lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col="blue", verticals=TRUE, pch=46)
  legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
dev.off()



#plot the ROH cumulative dstribution (for length ?)

##########################################################################################################################
# PLOT 8: DAF spectrum for different functional annotations

rm(list=ls())
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
#outbred
o_pops <- c("TSI","CEU")
#isolates
i_pops <- c("FVG","VBI")
all_pops <- c(i_pops,o_pops)

chr <- "22"
source('/nfs/users/nfs_m/mc14/Work/r_scripts/maf_bins_splitter.r')
# maf_classes <- c(0,0.01,0.02,0.05,0.10,0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1)
maf_classes <- c(0,0.05,0.10,0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1)
conseq_classes <- paste(base_folder,"INPUT_FILES/consequences.list",sep="")
conseq <- read.table(conseq_classes,header=F)
conseq$V1 <- as.character(conseq$V1)
# conseq <- c("missense_variant","synonymous_variant")
# cons <- conseq$V1[1]
# pop <- "FVG"

# for (cons in conseq) {
for (cons in conseq$V1) {
  for (pop in o_pops){
    file_path <- paste(base_folder,"INPUT_FILES/CHR",chr,"_no_fixed/",pop,"/",sep="")
    file_name <- paste(pop,".",cons,".",chr,".tab.gz",sep="")
    pop_table <- read.table(paste(file_path,file_name,sep=""),header=T,stringsAsFactors=F, comment.char="")
    #replace header because we know already how it is and we need to have a column called "CHROM", it's MANDATORY!
    colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","REC","ALC","DAC","MAC","DAF","MAF")

    if(length(pop_table$DAF) > 0) {

      pop_table$DAC <- as.numeric(as.character(pop_table$DAC))
      pop_table$MAC <- as.numeric(as.character(pop_table$MAC))
      pop_table$DAF <- as.numeric(as.character(pop_table$DAF))
      pop_table$MAF <- as.numeric(as.character(pop_table$MAF))
      
      #write a cute output
      outdir <- paste(pop,"_",chr,sep="")
      system(paste("mkdir -p ",pop,"_",chr,"/summaries",sep=""))

      #remove sites without the DAF info
      dim(pop_table)
      pop_table_no_na <- pop_table[-which(is.na(pop_table$DAF)),]
      dim(pop_table_no_na)
      summary(pop_table_no_na)
      pop_table_no_mono <- pop_table_no_na[-which(pop_table_no_na$DAF ==1 | pop_table_no_na$DAF == 0),]
      all_maf_classes <- split_bins(maf_classes,pop_table_no_mono,file_name,"DAF",outdir)
      gc()
      
      sink(paste(outdir,"/summaries/",file_name,'_maf_bin_resume.txt',sep=""))
      print(all_maf_classes)
      sink()
    }

  }
  
}


########################################
#Plot 9: some check stats on kinship data
i_pops <- c("FVG","VBI")

for (i_pop in i_pops){
  if(i_pop == "FVG"){
    new_kinship <- read.table("/lustre/scratch113/projects/fvg_seq/20140319/20140401_VQSR2.5_reapply_v138_excl/20140517_RELEASE/PLINK/KINSHIP/FVG_20140517.all.no_1st_deg.vcf.ibs0",header=T)
  }else if (i_pop == "VBI"){
    new_kinship <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/PLINK/KINSHIP/VBI_20140518.all.no_1st_deg.vcf.ibs0",header=T)
  }

  #check how many 1st degree samples we have still
  head(new_kinship[order(new_kinship$Kinship,decreasing=T),])

  # summaries
  summary(new_kinship)

  #boxplot for kinship
  jpeg(paste(base_folder,"PLOTS/9_kinship_",i_pop,".jpg",sep=""),width=800, height=800)
    par(lwd=2)
    boxplot(new_kinship$Kinship)
    # plot(M_ibd_CEU,CEU_tot_ibd$IBD_tot,main="",xlab="IBD segments length per couple", xlim=c(0,max(xmax)), verticals=TRUE, pch=46)
    # lines(M_ibd_FVG,FVG_tot_ibd$IBD_tot,col="red", verticals=TRUE, pch=46)
    # lines(M_ibd_TSI,TSI_tot_ibd$IBD_tot,col="green", verticals=TRUE, pch=46)
    # lines(M_ibd_VBI,VBI_tot_ibd$IBD_tot,col="blue", verticals=TRUE, pch=46)
    # legend("bottomright",pch =c(rep(19,length(pops))),legend=c("CEU","FVG","TSI","VBI"),col=c("black","red","green","blue"),ncol=4)
  dev.off()

  #plot also a heatmap

  
}



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