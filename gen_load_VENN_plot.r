#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
# base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/INPUT_FILES/FIVE_POPS/WG/"
chr <- "20"

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

##################update 14/05/2015#############################
# We want to plot a VENN diagramm matching the definition we gave of shred/novel/private,
# based on frequency, so we need to filter out,for each population all sites with MAF==0 (what about NA??)
# we have to read the tablewith our population data
merged_daf <- read.table(paste(base_folder,"CHR",chr,"/INGI_chr",chr,".merged_maf.tab.gz",sep=""),header=T,comment.char = "")

merged_daf_diag <- merged_daf[,-c(1:4)]

merged_daf_diag[merged_daf_diag>0]=1 #this is the fastest way to remove all sites with maf==0
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
        for(l in c(0,1)){
          if(any(c(i,j,k,h,l)==1)){
            res=length(which(merged_daf_diag$CEU==i & merged_daf_diag$TSI==j & merged_daf_diag$VBI==k & merged_daf_diag$FVG==h & merged_daf_diag$CARL==h))
            res.all2=rbind(res.all2,c(res,paste(which(c(i,j,k,h,l)==1),collapse="")))
          }
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

pops <- colnames(merged_daf_diag)
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
all_cols <-col_pop(pops)

#write the table:
write.table(res.all3,file=paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/VENN/All_pop_intersect_count_chr",chr,".txt",sep=""),sep="\t",col.names=T,quote=F,row.names=F)

#find all indexes
idx1=grep("1",res.all2[,2])
idx2=grep("2",res.all2[,2])
idx3=grep("3",res.all2[,2])
idx4=grep("4",res.all2[,2])
idx5=grep("5",res.all2[,2])

v <- draw.quintuple.venn(area1=as.numeric(sum(as.numeric(res.all2[idx1,1]))),
  area2=as.numeric(sum(as.numeric(res.all2[idx2,1]))),
  area3=as.numeric(sum(as.numeric(res.all2[idx3,1]))),
  area4=as.numeric(sum(as.numeric(res.all2[idx4,1]))),
  area5=as.numeric(sum(as.numeric(res.all2[idx5,1]))),
  n12=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx2),1]))),
  n13=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx3),1]))),
  n14=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx4),1]))),
  n15=as.numeric(sum(as.numeric(res.all2[intersect(idx1,idx5),1]))),
  n23=as.numeric(sum(as.numeric(res.all2[intersect(idx2,idx3),1]))),
  n24=as.numeric(sum(as.numeric(res.all2[intersect(idx2,idx4),1]))),
  n25=as.numeric(sum(as.numeric(res.all2[intersect(idx2,idx5),1]))),
  n34=as.numeric(sum(as.numeric(res.all2[intersect(idx3,idx4),1]))),
  n35=as.numeric(sum(as.numeric(res.all2[intersect(idx3,idx5),1]))),
  n45=as.numeric(sum(as.numeric(res.all2[intersect(idx4,idx5),1]))),
  n123=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx2),idx3),1]))),
  n124=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx2),idx4),1]))),
  n125=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx2),idx5),1]))),
  n134=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx3),idx4),1]))),
  n135=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx3),idx5),1]))),
  n145=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx1,idx4),idx5),1]))),
  n234=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx2,idx3),idx4),1]))),
  n235=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx2,idx3),idx5),1]))),
  n245=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx2,idx4),idx5),1]))),
  n345=as.numeric(sum(as.numeric(res.all2[intersect(intersect(idx3,idx4),idx5),1]))),
  n1234=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(idx1,idx2),idx3),idx4),1]))),
  n1235=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(idx1,idx2),idx3),idx5),1]))),
  n1245=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(idx1,idx2),idx4),idx5),1]))),
  n1345=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(idx1,idx3),idx4),idx5),1]))),
  n2345=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(idx2,idx3),idx4),idx5),1]))),
  n12345=as.numeric(sum(as.numeric(res.all2[intersect(intersect(intersect(intersect(idx1,idx2),idx3),idx4),idx5),1]))),
  ind=FALSE, category=pops,fill=all_cols)

jpeg(paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/RESULTS/VENN/2_5pop_DAF_VENN_chr",chr,".jpg",sep=""),width=800, height=800,pointsize = 20)
  grid.draw(v)
dev.off()



n1235, n1245, n1345, n2345, n12345

