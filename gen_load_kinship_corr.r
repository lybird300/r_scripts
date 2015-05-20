#Calculate correlation kinship gen-ped
rm(list=ls())
library(kinship2)

fvg_base_ped <- "/nfs/users/nfs_m/mc14/Work/SANGER/FVG/pedigree_FVG/VILLAGES"
# fvr_ped <- read.table('/home/max/Work/FVG/pedigree_FVG/SEQ/resia_ped.csv.seq.list')
fvr_ped <- read.table(paste(fvg_base_ped,'/resia_ped_sorted_manmod.csv',sep=""))
#had to manually fix the pedigree data
fvr_pedigree <- pedigree(id=fvr_ped$V1,dadid=fvr_ped$V2,momid=fvr_ped$V3,sex=fvr_ped$V4,miss=0)
fvr_kinship <- kinship(fvr_pedigree)

fve_ped <- read.table(paste(fvg_base_ped,'/erto_ped.csv',sep=""),skip=1)
fve_pedigree <- pedigree(id=fve_ped$V1,dadid=fve_ped$V2,momid=fve_ped$V3,sex=fve_ped$V4,miss=0)
#I had to manually remove mothers from DUOS relationships (id with no father but mother info)
fve_kinship <- kinship(fve_pedigree)

fvi_ped <- read.table(paste(fvg_base_ped,'/illegio_ped_uniq.csv',sep=""))
fvi_pedigree <- pedigree(id=fvi_ped$V1,dadid=fvi_ped$V2,momid=fvi_ped$V3,sex=fvi_ped$V4,miss=0)
#added manually founders
fvi_kinship <- kinship(fvi_pedigree)

fvs_ped <- read.table(paste(fvg_base_ped,'/sauris_ped_uniq.csv',sep=""))
fvs_pedigree <- pedigree(id=fvs_ped$V1,dadid=fvs_ped$V2,momid=fvs_ped$V3,sex=fvs_ped$V4,miss=0)
#removed duplicates founders
fvs_kinship <- kinship(fvs_pedigree)

################################################
#now load genomic kinship done with king
# fvg_seq_list <- read.table('/home/max/Work/FVG/SEQ_CALLING/FVG_seq_sorted.list',h=F)
fvg_seq_list <- read.table('/lustre/scratch113/teams/soranzo/users/mc14/INGI_FVG/SEQ_CALLING/FVG_seq_sorted.singlecol.keeplist',h=F)
colnames(fvg_seq_list) <- c("ID")
fvg_seq_list$ID <- as.character(fvg_seq_list$ID)


#load genomic kinship previously calculated
# base_folder <- "/home/max/Work/Analyses/PURGE_INBREEDING/KINSHIP"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP"

fve_kinship_king_seq <- read.table(paste(base_folder,"/FVG_Erto.keeplist.ibs0",sep=""),header=T)
fvi_kinship_king_seq <- read.table(paste(base_folder,"/FVG_Illegio.keeplist.ibs0",sep=""),header=T)
fvr_kinship_king_seq <- read.table(paste(base_folder,"/FVG_Resia.keeplist.ibs0",sep=""),header=T)
fvs_kinship_king_seq <- read.table(paste(base_folder,"/FVG_Sauris.keeplist.ibs0",sep=""),header=T)

head(fve_kinship_king_seq)
head(fvi_kinship_king_seq)
head(fvr_kinship_king_seq)
head(fvs_kinship_king_seq)

#load conversion table data
# conv_fvg <- read.table("/home/max/Work/Analyses/PURGE_INBREEDING/WGS_FVG_id_conversion_complete_sorted.list")
conv_fvg <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/WGS_FVG_id_conversion_complete_sorted.list")
colnames(conv_fvg) <- c("SEQ","CLIN")
conv_fvg$SEQ <- as.character(conv_fvg$SEQ)
conv_fvg$CLIN <- as.character(conv_fvg$CLIN)

#clean tables from useless columns
kin_set <- c("fve_kinship_king_seq","fvi_kinship_king_seq","fvr_kinship_king_seq","fvs_kinship_king_seq")
for (set in kin_set){
	current_set <- get(set)
	current_set$FID1 <- NULL
	current_set$FID2 <- NULL
	current_set$N_IBS0 <- NULL
	current_set$N_IBS1 <- NULL
	current_set$N_IBS2 <- NULL
	current_set$IBS <- NULL
	current_set$SE_IBS <- NULL
	current_set$N_HetHet <- NULL
	current_set$N_Het1 <- NULL
	current_set$N_Het2 <- NULL
	current_set$Distance <- NULL
	current_set$SE_Dist <- NULL
	# colnames(current_set) <- c("ID1","ID2","KIN","KEY")
	#convert ids
	current_merged1 <- merge(current_set,conv_fvg,by.x="ID1",by.y="SEQ",all.x)
	current_merged <- merge(current_merged1,conv_fvg,by.x="ID2",by.y="SEQ",all.x)
	assign(set,current_merged)
}


library(reshape2)
# kin_ped_set <- c("fve_kinship_seq","fvi_kinship_seq","fvr_kinship_seq","fvs_kinship_seq","fve_kinship","fvi_kinship","fvr_kinship","fvs_kinship")
kin_ped_set <- c("fve_kinship","fvi_kinship","fvr_kinship","fvs_kinship")
for (ped_set in kin_ped_set){
	current_ped_set <- get(ped_set)
	#convert ids
	current_ped_set_melted <- melt(current_ped_set, varnames = c('ID1', 'ID2'), na.rm = TRUE)
	# current_ped_set_melted <- current_ped_set_melted[which(current_ped_set_melted$ID1 != current_ped_set_melted$ID2),]
	# current_ped_set_melted <- current_ped_set_melted[which(paste(current_ped_set_melted$ID1,current_ped_set_melted$ID2,sep="_") != paste(current_ped_set_melted$ID2,current_ped_set_melted$ID1,sep="_")),]
	assign(paste(ped_set,"_melted",sep=""),current_ped_set_melted)
}

#split kinships in villages extracting only sequenced samples
# fvr_kinship_ped_seq_melted2 <- fvr_kinship_melted[fvr_kinship_melted$ID1 %in% fvg_seq_list$ID,]
# fvr_kinship_ped_seq_melted <- fvr_kinship_ped_seq_melted2[fvr_kinship_ped_seq_melted2$ID2 %in% fvg_seq_list$ID,]
to_r_seq_fve <- which(!(fve_kinship_king_seq$CLIN.x %in% fvg_seq_list$ID) | !(fve_kinship_king_seq$CLIN.y %in% fvg_seq_list$ID))
to_r_seq_fvi <- which(!(fvi_kinship_king_seq$CLIN.x %in% fvg_seq_list$ID) | !(fvi_kinship_king_seq$CLIN.y %in% fvg_seq_list$ID))
to_r_seq_fvr <- which(!(fvr_kinship_king_seq$CLIN.x %in% fvg_seq_list$ID) | !(fvr_kinship_king_seq$CLIN.y %in% fvg_seq_list$ID))
to_r_seq_fvs <- which(!(fvs_kinship_king_seq$CLIN.x %in% fvg_seq_list$ID) | !(fvs_kinship_king_seq$CLIN.y %in% fvg_seq_list$ID))
#just check that the previous index arrays are empty

to_r_ped_fve <- which(!(fve_kinship_melted$ID1 %in% fvg_seq_list$ID) | !(fve_kinship_melted$ID2 %in% fvg_seq_list$ID))
to_r_ped_fvi <- which(!(fvi_kinship_melted$ID1 %in% fvg_seq_list$ID) | !(fvi_kinship_melted$ID2 %in% fvg_seq_list$ID))
to_r_ped_fvr <- which(!(fvr_kinship_melted$ID1 %in% fvg_seq_list$ID) | !(fvr_kinship_melted$ID2 %in% fvg_seq_list$ID))
to_r_ped_fvs <- which(!(fvs_kinship_melted$ID1 %in% fvg_seq_list$ID) | !(fvs_kinship_melted$ID2 %in% fvg_seq_list$ID))
 
fve_kinship_ped_seq_melted <- fve_kinship_melted[-to_r_ped_fve,]
fvi_kinship_ped_seq_melted <- fvi_kinship_melted[-to_r_ped_fvi,]
fvr_kinship_ped_seq_melted <- fvr_kinship_melted[-to_r_ped_fvr,]
fvs_kinship_ped_seq_melted <- fvs_kinship_melted[-to_r_ped_fvs,]

# fvr_kinship_ped_seq_melted <- fvr_kinship_melted[fvr_kinship_melted$ID1 %in% fvg_seq_list$ID | fvr_kinship_melted$ID2 %in% fvg_seq_list$ID,]
# fve_kinship_ped_seq_melted <- fve_kinship_melted[fve_kinship_melted$ID1 %in% fvg_seq_list$ID | fve_kinship_melted$ID2 %in% fvg_seq_list$ID,]
# fvi_kinship_ped_seq_melted <- fvi_kinship_melted[fvi_kinship_melted$ID1 %in% fvg_seq_list$ID | fvi_kinship_melted$ID2 %in% fvg_seq_list$ID,]
# fvs_kinship_ped_seq_melted <- fvs_kinship_melted[fvs_kinship_melted$ID1 %in% fvg_seq_list$ID | fvs_kinship_melted$ID2 %in% fvg_seq_list$ID,]
# dim(fvr_kinship_ped_seq_melted)
# dim(fve_kinship_ped_seq_melted)
# dim(fvi_kinship_ped_seq_melted)
# dim(fvs_kinship_ped_seq_melted)

# dim(fvr_kinship)
# dim(fve_kinship)
# dim(fvi_kinship)
# dim(fvs_kinship)

dim(fvr_kinship_ped_seq_melted)
dim(fve_kinship_ped_seq_melted)
dim(fvi_kinship_ped_seq_melted)
dim(fvs_kinship_ped_seq_melted)
dim(fve_kinship_melted)
dim(fvi_kinship_melted)
dim(fvr_kinship_melted)
dim(fvs_kinship_melted)

# kin_ped_set_melted <- c("fve_kinship_seq_melted","fvi_kinship_seq_melted","fvr_kinship_seq_melted","fvs_kinship_seq_melted","fve_kinship_melted","fvi_kinship_melted","fvr_kinship_melted","fvs_kinship_melted")
kin_ped_set_melted <- c("fvr_kinship_ped_seq_melted","fve_kinship_ped_seq_melted","fvi_kinship_ped_seq_melted","fvs_kinship_ped_seq_melted","fve_kinship_melted","fvi_kinship_melted","fvr_kinship_melted","fvs_kinship_melted")

for (ped_set in kin_ped_set_melted){
	current_ped_set <- get(ped_set)
	#write tables
	write.table(current_ped_set,file=paste(getwd(),"KIN_PED",ped_set,sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
}


dim(fve_kinship_king_seq)
dim(fve_kinship_ped_seq_melted)
fve_kinship_ped_seq_melted$key <- paste(apply(fve_kinship_ped_seq_melted[,1:2],1,min),apply(fve_kinship_ped_seq_melted[,1:2],1,max),sep="_")
fve_kinship_ped_seq_melted <- fve_kinship_ped_seq_melted[!duplicated(fve_kinship_ped_seq_melted$key),]
fve_kinship_ped_seq_melted <- fve_kinship_ped_seq_melted[!fve_kinship_ped_seq_melted$ID1 == fve_kinship_ped_seq_melted$ID2,]
dim(fve_kinship_ped_seq_melted)

dim(fvi_kinship_king_seq)
dim(fvi_kinship_ped_seq_melted)
fvi_kinship_ped_seq_melted$key <- paste(apply(fvi_kinship_ped_seq_melted[,1:2],1,min),apply(fvi_kinship_ped_seq_melted[,1:2],1,max),sep="_")
fvi_kinship_ped_seq_melted <- fvi_kinship_ped_seq_melted[!duplicated(fvi_kinship_ped_seq_melted$key),]
fvi_kinship_ped_seq_melted <- fvi_kinship_ped_seq_melted[!fvi_kinship_ped_seq_melted$ID1 == fvi_kinship_ped_seq_melted$ID2,]
dim(fvi_kinship_ped_seq_melted)

dim(fvr_kinship_king_seq)
dim(fvr_kinship_ped_seq_melted)
fvr_kinship_ped_seq_melted$key <- paste(apply(fvr_kinship_ped_seq_melted[,1:2],1,min),apply(fvr_kinship_ped_seq_melted[,1:2],1,max),sep="_")
fvr_kinship_ped_seq_melted <- fvr_kinship_ped_seq_melted[!duplicated(fvr_kinship_ped_seq_melted$key),]
fvr_kinship_ped_seq_melted <- fvr_kinship_ped_seq_melted[!fvr_kinship_ped_seq_melted$ID1 == fvr_kinship_ped_seq_melted$ID2,]
dim(fvr_kinship_ped_seq_melted)

dim(fvs_kinship_king_seq)
dim(fvs_kinship_ped_seq_melted)
fvs_kinship_ped_seq_melted$key <- paste(apply(fvs_kinship_ped_seq_melted[,1:2],1,min),apply(fvs_kinship_ped_seq_melted[,1:2],1,max),sep="_")
fvs_kinship_ped_seq_melted <- fvs_kinship_ped_seq_melted[!duplicated(fvs_kinship_ped_seq_melted$key),]
fvs_kinship_ped_seq_melted <- fvs_kinship_ped_seq_melted[!fvs_kinship_ped_seq_melted$ID1 == fvs_kinship_ped_seq_melted$ID2,]
dim(fvs_kinship_ped_seq_melted)

########trying to extract the same samples from ped kinship and seq kinship....

length(unique(c(fve_kinship_ped_seq_melted$ID1,fve_kinship_ped_seq_melted$ID2)))
length(unique(c(fve_kinship_king_seq$CLIN.x,fve_kinship_king_seq$CLIN.y)))


length(unique(c(fvi_kinship_ped_seq_melted$ID1,fvi_kinship_ped_seq_melted$ID2)))
length(unique(c(fvi_kinship_king_seq$CLIN.x,fvi_kinship_king_seq$CLIN.y)))

length(unique(c(fvr_kinship_ped_seq_melted$ID1,fvr_kinship_ped_seq_melted$ID2)))
length(unique(c(fvr_kinship_king_seq$CLIN.x,fvr_kinship_king_seq$CLIN.y)))

length(unique(c(fvs_kinship_ped_seq_melted$ID1,fvs_kinship_ped_seq_melted$ID2)))
length(unique(c(fvs_kinship_king_seq$CLIN.x,fvs_kinship_king_seq$CLIN.y)))

############################# we don't have the same individuals, it seeems..not all sequenced samples are in the pedigree kinship ###################################################
fve_ped_list <- c(fve_kinship_ped_seq_melted$ID1, fve_kinship_ped_seq_melted$ID2)
to_r_seq_fve <- which(!(fve_kinship_king_seq$CLIN.x %in% fve_ped_list) | !(fve_kinship_king_seq$CLIN.y %in% fve_ped_list))
fve_kinship_king_seq2 <- fve_kinship_king_seq[-to_r_seq_fve,]

fvi_ped_list <- c(fvi_kinship_ped_seq_melted$ID1, fvi_kinship_ped_seq_melted$ID2)
to_r_seq_fvi <- which(!(fvi_kinship_king_seq$CLIN.x %in% fvi_ped_list) | !(fvi_kinship_king_seq$CLIN.y %in% fvi_ped_list))
fvi_kinship_king_seq2 <- fvi_kinship_king_seq[-to_r_seq_fvi,]

fvr_ped_list <- c(fvr_kinship_ped_seq_melted$ID1, fvr_kinship_ped_seq_melted$ID2)
to_r_seq_fvr <- which(!(fvr_kinship_king_seq$CLIN.x %in% fvr_ped_list) | !(fvr_kinship_king_seq$CLIN.y %in% fvr_ped_list))

if(length(to_r_seq_fvr) > 0) {
	fvr_kinship_king_seq2 <- fvr_kinship_king_seq[-to_r_seq_fvr,]
	
}else{
	fvr_kinship_king_seq2 <- fvr_kinship_king_seq
}

fvs_ped_list <- c(fvs_kinship_ped_seq_melted$ID1, fvs_kinship_ped_seq_melted$ID2)
to_r_seq_fvs <- which(!(fvs_kinship_king_seq$CLIN.x %in% fvs_ped_list) | !(fvs_kinship_king_seq$CLIN.y %in% fvs_ped_list))
fvs_kinship_king_seq2 <- fvs_kinship_king_seq[-to_r_seq_fvs,]


# fve_kinship_king_seq2 <- (fve_kinship_king_seq[which(fve_kinship_king_seq$CLIN.x %in% fve_kinship_ped_seq_melted$ID1),][which(fve_kinship_king_seq$CLIN.y %in% fve_kinship_ped_seq_melted$ID2),])
# fve_kinship_king_seq2 <- (fve_kinship_king_seq[which(fve_kinship_king_seq$CLIN.x %in% fve_kinship_ped_seq_melted$ID1),][which(fve_kinship_king_seq$CLIN.y %in% fve_kinship_ped_seq_melted$ID2),])
# fve_kinship_king_seq2 <- fve_kinship_king_seq2[!is.na(fve_kinship_king_seq2$ID2),]

# fvi_kinship_king_seq2 <- (fvi_kinship_king_seq[which(fvi_kinship_king_seq$CLIN.x %in% fvi_kinship_ped_seq_melted$ID1),][which(fvi_kinship_king_seq$CLIN.y %in% fvi_kinship_ped_seq_melted$ID2),])
# fvi_kinship_king_seq2 <- fvi_kinship_king_seq2[!is.na(fvi_kinship_king_seq2$ID2),]

# fvr_kinship_king_seq2 <- (fvr_kinship_king_seq[which(fvr_kinship_king_seq$CLIN.x %in% fvr_kinship_ped_seq_melted$ID1),][which(fvr_kinship_king_seq$CLIN.y %in% fvr_kinship_ped_seq_melted$ID2),])
# fvr_kinship_king_seq2 <- fvr_kinship_king_seq2[!is.na(fvr_kinship_king_seq2$ID2),]

# fvs_kinship_king_seq2 <- (fvs_kinship_king_seq[which(fvs_kinship_king_seq$CLIN.x %in% fvs_kinship_ped_seq_melted$ID1),][which(fvs_kinship_king_seq$CLIN.y %in% fvs_kinship_ped_seq_melted$ID2),])
# fvs_kinship_king_seq2 <- fvs_kinship_king_seq2[!is.na(fvs_kinship_king_seq2$ID2),]


length(unique(c(fve_kinship_king_seq2$CLIN.x,fve_kinship_king_seq2$CLIN.y)))
length(unique(c(fve_kinship_ped_seq_melted$ID1,fve_kinship_ped_seq_melted$ID2)))

length(unique(c(fvi_kinship_king_seq2$CLIN.x,fvi_kinship_king_seq2$CLIN.y)))
length(unique(c(fvi_kinship_ped_seq_melted$ID1,fvi_kinship_ped_seq_melted$ID2)))
# fvi_kinship_ped_seq_melted[fvi_kinship_ped_seq_melted$ID1 == fvi_kinship_ped_seq_melted$ID2,]
# fvi_kinship_king_seq2[fvi_kinship_king_seq2$CLIN.x == fvi_kinship_king_seq2$CLIN.y,]

length(unique(c(fvr_kinship_king_seq2$CLIN.x,fvr_kinship_king_seq2$CLIN.y)))
length(unique(c(fvr_kinship_ped_seq_melted$ID1,fvr_kinship_ped_seq_melted$ID2)))

length(unique(c(fvs_kinship_king_seq2$CLIN.x,fvs_kinship_king_seq2$CLIN.y)))
length(unique(c(fvs_kinship_ped_seq_melted$ID1,fvs_kinship_ped_seq_melted$ID2)))

############now we have everything with the same samples .... ###############


fve_kinship_king_seq2$key <- paste(apply(fve_kinship_king_seq2[,4:5],1,min),apply(fve_kinship_king_seq2[,4:5],1,max),sep="_")
fvi_kinship_king_seq2$key <- paste(apply(fvi_kinship_king_seq2[,4:5],1,min),apply(fvi_kinship_king_seq2[,4:5],1,max),sep="_")
fvr_kinship_king_seq2$key <- paste(apply(fvr_kinship_king_seq2[,4:5],1,min),apply(fvr_kinship_king_seq2[,4:5],1,max),sep="_")
fvs_kinship_king_seq2$key <- paste(apply(fvs_kinship_king_seq2[,4:5],1,min),apply(fvs_kinship_king_seq2[,4:5],1,max),sep="_")
###############################################################

# fve_kinship_king_seq3 <- fve_kinship_king_seq2[!duplicated(fve_kinship_king_seq2$key),]
# fvi_kinship_king_seq3 <- fvi_kinship_king_seq2[!duplicated(fvi_kinship_king_seq2$key),]
# fvr_kinship_king_seq3 <- fvr_kinship_king_seq2[!duplicated(fvr_kinship_king_seq2$key),]
# fvs_kinship_king_seq3 <- fvs_kinship_king_seq2[!duplicated(fvs_kinship_king_seq2$key),]


# length(unique(fve_kinship_king_seq3$CLIN.y))
# length(unique(fvi_kinship_king_seq3$CLIN.y))
# length(unique(fvr_kinship_king_seq3$CLIN.y))
# length(unique(fvs_kinship_king_seq3$CLIN.y))

# length(unique(fve_kinship_ped_seq_melted$ID2))
# length(unique(fvi_kinship_ped_seq_melted$ID2))
# length(unique(fvr_kinship_ped_seq_melted$ID2))
# length(unique(fvs_kinship_ped_seq_melted$ID2))

#sort dataframes
fve_kinship_king_seq2 <- fve_kinship_king_seq2[with(fve_kinship_king_seq2,order(key)),]
fve_kinship_ped_seq_melted <- fve_kinship_ped_seq_melted[with(fve_kinship_ped_seq_melted,order(key)),]

fvi_kinship_king_seq2 <- fvi_kinship_king_seq2[with(fvi_kinship_king_seq2,order(key)),]
fvi_kinship_ped_seq_melted <- fvi_kinship_ped_seq_melted[with(fvi_kinship_ped_seq_melted,order(key)),]

fvr_kinship_king_seq2 <- fvr_kinship_king_seq2[with(fvr_kinship_king_seq2,order(key)),]
fvr_kinship_ped_seq_melted <- fvr_kinship_ped_seq_melted[with(fvr_kinship_ped_seq_melted,order(key)),]

fvs_kinship_king_seq2 <- fvs_kinship_king_seq2[with(fvs_kinship_king_seq2,order(key)),]
fvs_kinship_ped_seq_melted <- fvs_kinship_ped_seq_melted[with(fvs_kinship_ped_seq_melted,order(key)),]

dim(fvs_kinship_king_seq2)
dim(fvs_kinship_ped_seq_melted)

#####now write this final kinbship tables to use for correlation checks
king_kinship_set <- c("fve_kinship_king_seq2","fvi_kinship_king_seq2","fvr_kinship_king_seq2","fvs_kinship_king_seq2")
ped_kinship_set <- c("fve_kinship_ped_seq_melted","fvi_kinship_ped_seq_melted","fvr_kinship_ped_seq_melted","fvs_kinship_ped_seq_melted")

for (ped_set in king_kinship_set){
	current_ped_set <- get(ped_set)
	current_ped_set$ID1 <- NULL
	current_ped_set$ID2 <- NULL
	#write tables
	write.table(current_ped_set,file=paste(getwd(),"/KIN_GEN/",ped_set,".final_overlap",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
}
for (ped_set in ped_kinship_set){
	current_ped_set <- get(ped_set)
	#write tables
	write.table(current_ped_set,file=paste(getwd(),"/KIN_PED/",ped_set,".final_overlap",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
}
#####################################################################################
#Correlation between kinship

fve_kinship_gen <- read.table(paste(getwd(),"/KIN_GEN/fve_kinship_king_seq2.final_overlap",sep=""),header=F)
fvi_kinship_gen <- read.table(paste(getwd(),"/KIN_GEN/fvi_kinship_king_seq2.final_overlap",sep=""),header=F)
fvr_kinship_gen <- read.table(paste(getwd(),"/KIN_GEN/fvr_kinship_king_seq2.final_overlap",sep=""),header=F)
fvs_kinship_gen <- read.table(paste(getwd(),"/KIN_GEN/fvs_kinship_king_seq2.final_overlap",sep=""),header=F)

fve_kinship_ped <- read.table(paste(getwd(),"/KIN_PED/fve_kinship_ped_seq_melted.final_overlap",sep=""),header=F)
fvi_kinship_ped <- read.table(paste(getwd(),"/KIN_PED/fvi_kinship_ped_seq_melted.final_overlap",sep=""),header=F)
fvr_kinship_ped <- read.table(paste(getwd(),"/KIN_PED/fvr_kinship_ped_seq_melted.final_overlap",sep=""),header=F)
fvs_kinship_ped <- read.table(paste(getwd(),"/KIN_PED/fvs_kinship_ped_seq_melted.final_overlap",sep=""),header=F)


colnames(fve_kinship_gen) <- c("KIN","ID1","ID2","KEY")
colnames(fve_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fve_gen_ped_merged <- merge(fve_kinship_gen,fve_kinship_ped,by.x="KEY",by.y="KEY",all.x)


colnames(fvi_kinship_gen) <- c("KIN","ID1","ID2","KEY")
colnames(fvi_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fvi_gen_ped_merged <- merge(fvi_kinship_gen,fvi_kinship_ped,by.x="KEY",by.y="KEY",all.x)

colnames(fvr_kinship_gen) <- c("KIN","ID1","ID2","KEY")
colnames(fvr_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fvr_gen_ped_merged <- merge(fvr_kinship_gen,fvr_kinship_ped,by.x="KEY",by.y="KEY",all.x)

colnames(fvs_kinship_gen) <- c("KIN","ID1","ID2","KEY")
colnames(fvs_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fvs_gen_ped_merged <- merge(fvs_kinship_gen,fvs_kinship_ped,by.x="KEY",by.y="KEY",all.x)

#conversion table kinship2 / KING
# MZT: 0.5 > 0.3535
# 1st degree: 0.25 (0.1767,0.3535)
# 2nd degree: 0.125 (0.0883,0.1767)
# 3rd degree: 0.0625 (0.0441,0.0883)
# Unrelated: 0 < 0.0441


villages <- c("fve","fvi","fvr","fvs")
# villages <- c("fve","fvi","fvr","fvs","vbi","carl","fvg")

for (village in villages){
	current_village_set <- paste(village,"_gen_ped_merged",sep="")
	current_village <- get(current_village_set)
	current_village$KIN.cast <- 0
	current_village[current_village$KIN.gen >= 0.3535,"KIN.cast"] <- 0.5 
	current_village[current_village$KIN.gen < 0.3535 & current_village$KIN.gen >= 0.1767,"KIN.cast"] <- 0.25 
	current_village[current_village$KIN.gen < 0.1767 & current_village$KIN.gen >= 0.0883,"KIN.cast"] <- 0.125 
	current_village[current_village$KIN.gen < 0.0883 & current_village$KIN.gen >= 0.0441,"KIN.cast"] <- 0.0625 
	current_village[current_village$KIN.gen <= 0.0441,"KIN.cast"] <- 0
	colnames(current_village) <- c("KEY","KIN.gen","ID1.x","ID2.x","ID1.y","ID2.y","KIN.ped","KIN.cast")
	assign(current_village_set,current_village)
}

summary(fve_gen_ped_merged)
summary(fvi_gen_ped_merged)
summary(fvr_gen_ped_merged)
summary(fvs_gen_ped_merged)

library(Hmisc)
library(corrgram)
# corrgram(fve_gen_ped_merged[,c(4,7,8)],lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# corrgram(fve_gen_ped_merged[,c(4,7,8)],lower.panel=panel.ellipse,upper.panel=panel.pts, text.panel=panel.txt,diag.panel=panel.minmax)

rcorr(fve_gen_ped_merged$KIN.gen,fve_gen_ped_merged$KIN.ped)
rcorr(fve_gen_ped_merged$KIN.cast,fve_gen_ped_merged$KIN.ped)
corrgram(fve_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)

rcorr(fvi_gen_ped_merged$KIN.gen,fvi_gen_ped_merged$KIN.ped)
rcorr(fvi_gen_ped_merged$KIN.cast,fvi_gen_ped_merged$KIN.ped)
corrgram(fvi_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)

rcorr(fvr_gen_ped_merged$KIN.gen,fvr_gen_ped_merged$KIN.ped)
rcorr(fvr_gen_ped_merged$KIN.cast,fvr_gen_ped_merged$KIN.ped)
corrgram(fvr_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)

rcorr(fvs_gen_ped_merged$KIN.gen,fvs_gen_ped_merged$KIN.ped)
rcorr(fvs_gen_ped_merged$KIN.cast,fvs_gen_ped_merged$KIN.ped)
corrgram(fvs_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)

cor_fve <- cor.test(fve_gen_ped_merged$KIN.gen,fve_gen_ped_merged$KIN.ped,alternative='g')
cor_fvi <- cor.test(fvi_gen_ped_merged$KIN.gen,fvi_gen_ped_merged$KIN.ped,alternative='g')
cor_fvr <- cor.test(fvr_gen_ped_merged$KIN.gen,fvr_gen_ped_merged$KIN.ped,alternative='g')
cor_fvs <- cor.test(fvs_gen_ped_merged$KIN.gen,fvs_gen_ped_merged$KIN.ped,alternative='g')

cor_fve_c <- cor.test(fve_gen_ped_merged$KIN.cast,fve_gen_ped_merged$KIN.ped,alternative='g')
cor_fvi_c <- cor.test(fvi_gen_ped_merged$KIN.cast,fvi_gen_ped_merged$KIN.ped,alternative='g')
cor_fvr_c <- cor.test(fvr_gen_ped_merged$KIN.cast,fvr_gen_ped_merged$KIN.ped,alternative='g')
cor_fvs_c <- cor.test(fvs_gen_ped_merged$KIN.cast,fvs_gen_ped_merged$KIN.ped,alternative='g')

# > cor_fve_c

#         Pearson's product-moment correlation

# data:  fve_gen_ped_merged$KIN.cast and fve_gen_ped_merged$KIN.ped
# t = 48.5384, df = 1429, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.7719512 1.0000000
# sample estimates:
#       cor 
# 0.7889581 

# > cor_fvi_c

#         Pearson's product-moment correlation

# data:  fvi_gen_ped_merged$KIN.cast and fvi_gen_ped_merged$KIN.ped
# t = 95.1738, df = 1483, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.9207357 1.0000000
# sample estimates:
#       cor 
# 0.9269908 

# > cor_fvr_c

#         Pearson's product-moment correlation

# data:  fvr_gen_ped_merged$KIN.cast and fvr_gen_ped_merged$KIN.ped
# t = 56.9334, df = 2483, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.7378112 1.0000000
# sample estimates:
#       cor 
# 0.7524915 

# > cor_fvs_c

#         Pearson's product-moment correlation

# data:  fvs_gen_ped_merged$KIN.cast and fvs_gen_ped_merged$KIN.ped
# t = 38.4834, df = 1376, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.697926 1.000000
# sample estimates:
#       cor 
# 0.7199802 

# Kinship correlation values from cor.test(),Pearson's coefficient
# Fvg-E: r= 0.6517779, pval = <2.2e-16
# Fvg-I: r= 0.775517, pval = <2.2e-16
# Fvg-R: r= 0.5116486 , pval = <2.2e-16
# Fvg-S: r= 0.5232983, pval = <2.2e-16

# > cor_fve

#         Pearson's product-moment correlation

# data:  fve_gen_ped_merged$KIN.gen and fve_gen_ped_merged$KIN.ped
# t = 32.4872, df = 1429, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.6260273 1.0000000
# sample estimates:
#       cor 
# 0.6517779 

# > cor_fvi

#         Pearson's product-moment correlation

# data:  fvi_gen_ped_merged$KIN.gen and fvi_gen_ped_merged$KIN.ped
# t = 47.3051, df = 1483, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.7579145 1.0000000
# sample estimates:
#      cor 
# 0.775517 

# > cor_fvr

#         Pearson's product-moment correlation

# data:  fvr_gen_ped_merged$KIN.gen and fvr_gen_ped_merged$KIN.ped
# t = 29.6735, df = 2483, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.486866 1.000000
# sample estimates:
#       cor 
# 0.5116486 

# > cor_fvs

#         Pearson's product-moment correlation

# data:  fvs_gen_ped_merged$KIN.gen and fvs_gen_ped_merged$KIN.ped
# t = 22.7794, df = 1376, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.4903437 1.0000000
# sample estimates:
#       cor 
# 0.5232983 

###################Merge all FVG together
fvg_gen_ped_merged <- rbind(fve_gen_ped_merged,fvi_gen_ped_merged,fvr_gen_ped_merged,fvs_gen_ped_merged)
rcorr(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped)
cor(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped)
cor.test(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped,method="spearman")
cor_fvg <- cor.test(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped,alternative="g")

corrgram(fvg_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)

###########
# Kinship VB
# vbi_kinship_ped <- read.table("/home/max/Work/Analyses/PURGE_INBREEDING/KINSHIP/KIN_GEN/VBI.keeplist.ibs0.kinship.conv.kin.over_ped",header=F)
# vbi_pedigree <- pedigree(id=vbi_kinship_ped$V1,dadid=vbi_kinship_ped$V2,momid=vbi_kinship_ped$V3,sex=vbi_kinship_ped$V4,miss=0)
# vbi_kinship <- kinship(vbi_pedigree)

# library(reshape2)
# vbi_kinship_melted <- melt(vbi_kinship, varnames = c('ID1', 'ID2'), na.rm = TRUE)

# write.table(vbi_kinship_melted,file=paste(getwd(),"vbi_kinship_melted",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
vbi_kinship_ped <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/KIN_PED/vbi_kinship_melted.sorted.pkin.over_gen",header=F)
vbi_kinship_gen <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/KIN_GEN/VBI.keeplist.ibs0.kinship.conv.kin.over_ped",header=F)

colnames(vbi_kinship_gen) <- c("ID1","ID2","KIN","KEY")
colnames(vbi_kinship_ped) <- c("ID1","ID2","KIN","KEY")
vbi_gen_ped_merged <- merge(vbi_kinship_gen,vbi_kinship_ped,by.x="KEY",by.y="KEY",all.x)
colnames(vbi_gen_ped_merged) <- c("KEY","ID1.x","ID2.x","KIN.gen","ID1.y","ID2.y","KIN.ped")

# villages <- c("fve","fvi","fvr","fvs")
villages <- c("vbi")

for (village in villages){
	current_village_set <- paste(village,"_gen_ped_merged",sep="")
	current_village <- get(current_village_set)
	current_village$KIN.cast <- 0
	current_village[current_village$KIN.gen >= 0.3535,"KIN.cast"] <- 0.5 
	current_village[current_village$KIN.gen < 0.3535 & current_village$KIN.gen >= 0.1767,"KIN.cast"] <- 0.25 
	current_village[current_village$KIN.gen < 0.1767 & current_village$KIN.gen >= 0.0883,"KIN.cast"] <- 0.125 
	current_village[current_village$KIN.gen < 0.0883 & current_village$KIN.gen >= 0.0441,"KIN.cast"] <- 0.0625 
	current_village[current_village$KIN.gen <= 0.0441,"KIN.cast"] <- 0
	colnames(current_village) <- c("KEY","ID1.x","ID2.x","KIN.gen","ID1.y","ID2.y","KIN.ped","KIN.cast")
	assign(current_village_set,current_village)
}


cor_vbi <- cor.test(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped,alternative="g")
cor.test(vbi_gen_ped_merged$KIN.cast,vbi_gen_ped_merged$KIN.ped,alternative="g")
# > cor.test(vbi_gen_ped_merged$KIN.cast,vbi_gen_ped_merged$KIN.ped,alternative="g")

#         Pearson's product-moment correlation

# data:  vbi_gen_ped_merged$KIN.cast and vbi_gen_ped_merged$KIN.ped
# t = 119.7181, df = 12003, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.7307975 1.0000000
# sample estimates:
#       cor 
# 0.7377167 


# Kinship CARL
carl_kinship_ped <- read.table("/nfs/users/nfs_m/mc14/Work/SANGER/CARLANTINO/pedigree_CARL/carlantino_sorted_nofam.ped",header=F)
carl_pedigree <- pedigree(id=carl_kinship_ped$V1,dadid=carl_kinship_ped$V2,momid=carl_kinship_ped$V3,sex=carl_kinship_ped$V4,miss=0)
carl_kinship <- kinship(carl_pedigree)

library(reshape2)
carl_kinship_melted <- melt(carl_kinship, varnames = c('ID1', 'ID2'), na.rm = TRUE)

write.table(carl_kinship_melted,file=paste(getwd(),"carl_kinship_melted",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)

carl_kinship_ped <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/KIN_PED/carl_kinship_melted.sorted.pkin.over_gen",header=F)
carl_kinship_gen <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/KIN_GEN/CARL.keeplist.ibs0.kinship.conv.kin.over_ped",header=F)

colnames(carl_kinship_gen) <- c("ID1","ID2","KIN.gen","KEY")
colnames(carl_kinship_ped) <- c("ID1","ID2","KIN.ped","KEY")
carl_gen_ped_merged <- merge(carl_kinship_gen,carl_kinship_ped,by.x="KEY",by.y="KEY",all.x)

villages <- c("carl")

for (village in villages){
	current_village_set <- paste(village,"_gen_ped_merged",sep="")
	current_village <- get(current_village_set)
	current_village$KIN.cast <- 0
	current_village[current_village$KIN.gen >= 0.3535,"KIN.cast"] <- 0.5 
	current_village[current_village$KIN.gen < 0.3535 & current_village$KIN.gen >= 0.1767,"KIN.cast"] <- 0.25 
	current_village[current_village$KIN.gen < 0.1767 & current_village$KIN.gen >= 0.0883,"KIN.cast"] <- 0.125 
	current_village[current_village$KIN.gen < 0.0883 & current_village$KIN.gen >= 0.0441,"KIN.cast"] <- 0.0625 
	current_village[current_village$KIN.gen <= 0.0441,"KIN.cast"] <- 0
	colnames(current_village) <- c("KEY","ID1.x","ID2.x","KIN.gen","ID1.y","ID2.y","KIN.ped","KIN.cast")
	assign(current_village_set,current_village)
}


cor_carl <- cor.test(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped,alternative="g")
cor.test(carl_gen_ped_merged$KIN.cast,carl_gen_ped_merged$KIN.ped,alternative="g")
# corrgram(carl_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# > cor.test(carl_gen_ped_merged$KIN.cast,carl_gen_ped_merged$KIN.ped,alternative="g")

#         Pearson's product-moment correlation

# data:  carl_gen_ped_merged$KIN.cast and carl_gen_ped_merged$KIN.ped
# t = 198.8691, df = 4093, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.9494812 1.0000000
# sample estimates:
#       cor 
# 0.9519527 

# rcor_fvg <- rcorr(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped)
rcor_vbi <- rcorr(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped)
rcor_carl <- rcorr(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped)

# cov_fvg <- cov(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped)
cov_vbi <- cov(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped)
cov_carl <- cov(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped)

# Kinship correlation values from cor.test(),Pearson's coefficient
# Fvg: r= 0.611911,pval = <2.2e-16
# Vbi: r= 0.5021468,pval = <2.2e-16
# Carl: r= 0.8717032,pval = <2.2e-16

# > cor_vbi

# 	Pearson's product-moment correlation

# data:  vbi_gen_ped_merged$KIN.x and vbi_gen_ped_merged$KIN.y
# t = 63.6164, df = 12003, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.4908341 1.0000000
# sample estimates:
#       cor 
# 0.5021468 

# > cor_carl

# 	Pearson's product-moment correlation

# data:  carl_gen_ped_merged$KIN.x and carl_gen_ped_merged$KIN.y
# t = 113.8055, df = 4093, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.8653884 1.0000000
# sample estimates:
#       cor 
# 0.8717032 

# > 
# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvg_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
# 	corrgram(fvg_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# dev.off()

# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/vbi_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
# 	corrgram(vbi_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# dev.off()

# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/carl_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
# 	corrgram(carl_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# dev.off()

# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fve_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
# 	corrgram(fve_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# dev.off()
# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvi_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
# 	corrgram(fvi_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# dev.off()
# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvr_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
# 	corrgram(fvr_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# dev.off()
# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvs_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
# 	corrgram(fvs_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
# dev.off()

#read unrelated lists
vbi_unrel <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/VBI_unrelated.list")
vbi_conv <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140319/20140402_VQSR2.5_reapply_138_excl/20140518_RELEASE/VBI_SEQ2HSR_table.txt")
vbi_unrel$V1 <- as.character(vbi_unrel$V1)
colnames(vbi_conv) <- c("SEQ","CLIN")
vbi_unrel_conv <-vbi_conv[which(vbi_conv$SEQ %in% vbi_unrel$V1),]

carl_unrel <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/CARL_unrelated.list")
carl_conv <- read.table("/lustre/scratch113/projects/carl_seq/20140410/20140410_VQSR2.5_reapply_138_excl/20140710_RELEASE/conversion_table_EGa_clinic_id.txt")
carl_unrel$V1 <- as.character(carl_unrel$V1)
colnames(carl_conv) <- c("SEQ","CLIN")
carl_unrel_conv <- carl_conv[which(carl_conv$SEQ %in% carl_unrel$V1),]


fvg_unrel <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/FVG_unrelated.list")
fve_unrel <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Erto_unrelated.list")
fvi_unrel <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Illegio_unrelated.list")
fvr_unrel <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Resia_unrelated.list")
fvs_unrel <- read.table("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/listpop/Sauris_unrelated.list")

fvg_unrel$V1 <- as.character(fvg_unrel$V1)
fve_unrel$V1 <- as.character(fve_unrel$V1)
fvi_unrel$V1 <- as.character(fvi_unrel$V1)
fvr_unrel$V1 <- as.character(fvr_unrel$V1)
fvs_unrel$V1 <- as.character(fvs_unrel$V1)


# fvg_unrel_conv <- conv_fvg[which(fvg_unrel$V1 %in% conv_fvg$SEQ),]

fve_unrel_conv <- conv_fvg[which(conv_fvg$SEQ %in% fve_unrel$V1),]
fvi_unrel_conv <- conv_fvg[which(conv_fvg$SEQ %in% fvi_unrel$V1 ),]
fvr_unrel_conv <- conv_fvg[which(conv_fvg$SEQ %in% fvr_unrel$V1 ),]
fvs_unrel_conv <- conv_fvg[which(conv_fvg$SEQ %in% fvs_unrel$V1 ),]

dim(fvg_unrel_conv)

dim(fve_unrel_conv)
dim(fvi_unrel_conv)
dim(fvr_unrel_conv)
dim(fvs_unrel_conv)

# fvg_unrel_conv$key <- paste(apply(fvg_unrel_conv[,1:2],1,min),apply(fvg_unrel_conv[,1:2],1,max),sep="_")
# fve_unrel_conv$key <- paste(apply(fve_unrel_conv[,1:2],1,min),apply(fve_unrel_conv[,1:2],1,max),sep="_")
# fvi_unrel_conv$key <- paste(apply(fvi_unrel_conv[,1:2],1,min),apply(fvi_unrel_conv[,1:2],1,max),sep="_")
# fvr_unrel_conv$key <- paste(apply(fvr_unrel_conv[,1:2],1,min),apply(fvr_unrel_conv[,1:2],1,max),sep="_")
# fvs_unrel_conv$key <- paste(apply(fvs_unrel_conv[,1:2],1,min),apply(fvs_unrel_conv[,1:2],1,max),sep="_")

# fvg_unrel_conv <- fvg_unrel_conv[!duplicated(fvg_unrel_conv$key),]
# fve_unrel_conv <- fve_unrel_conv[!duplicated(fve_unrel_conv$key),]
# fvi_unrel_conv <- fvi_unrel_conv[!duplicated(fvi_unrel_conv$key),]
# fvr_unrel_conv <- fvr_unrel_conv[!duplicated(fvr_unrel_conv$key),]
# fvs_unrel_conv <- fvs_unrel_conv[!duplicated(fvs_unrel_conv$key),]

to_r_merged_vbi <- which(!(vbi_gen_ped_merged$ID1.x %in% vbi_unrel_conv$CLIN) | !(vbi_gen_ped_merged$ID2.x %in% vbi_unrel_conv$CLIN))
to_r_merged_carl <- which(!(carl_gen_ped_merged$ID1.x %in% carl_unrel_conv$CLIN) | !(carl_gen_ped_merged$ID2.x %in% carl_unrel_conv$CLIN))

vbi_gen_ped_merged_unrel <- vbi_gen_ped_merged[-to_r_merged_vbi,]
carl_gen_ped_merged_unrel <- carl_gen_ped_merged[-to_r_merged_carl,]

# fvg_gen_ped_merged_unrel <- subset(fvg_gen_ped_merged, fvg_gen_ped_merged$ID1.x %in% fvg_unrel_conv$CLIN | fvg_gen_ped_merged$ID2.x %in% fvg_unrel_conv$CLIN )

to_r_merged_fve <- which(!(fve_gen_ped_merged$ID1.x %in% fve_unrel_conv$CLIN) | !(fve_gen_ped_merged$ID2.x %in% fve_unrel_conv$CLIN))
to_r_merged_fvi <- which(!(fvi_gen_ped_merged$ID1.x %in% fvi_unrel_conv$CLIN) | !(fvi_gen_ped_merged$ID2.x %in% fvi_unrel_conv$CLIN))
to_r_merged_fvr <- which(!(fvr_gen_ped_merged$ID1.x %in% fvr_unrel_conv$CLIN) | !(fvr_gen_ped_merged$ID2.x %in% fvr_unrel_conv$CLIN))
to_r_merged_fvs <- which(!(fvs_gen_ped_merged$ID1.x %in% fvs_unrel_conv$CLIN) | !(fvs_gen_ped_merged$ID2.x %in% fvs_unrel_conv$CLIN))

fve_gen_ped_merged_unrel <- fve_gen_ped_merged[-to_r_merged_fve,]
fvi_gen_ped_merged_unrel <- fvi_gen_ped_merged[-to_r_merged_fvi,]
fvr_gen_ped_merged_unrel <- fvr_gen_ped_merged[-to_r_merged_fvr,]
fvs_gen_ped_merged_unrel <- fvs_gen_ped_merged[-to_r_merged_fvs,]

# if(length(to_r_seq_fvr) > 0) {
# 	fve_gen_ped_merged_unrel <- fve_gen_ped_merged[-to_r_merged_fve,]
# 	fve_gen_ped_merged_unrel <- fve_gen_ped_merged[-to_r_merged_fve,]
# 	fve_gen_ped_merged_unrel <- fve_gen_ped_merged[-to_r_merged_fve,]
# 	fve_gen_ped_merged_unrel <- fve_gen_ped_merged[-to_r_merged_fve,]
# }else{
# 	fve_gen_ped_merged_unrel <- fve_gen_ped_merged
# }

dim(fve_gen_ped_merged_unrel)
dim(fvi_gen_ped_merged_unrel)
dim(fvr_gen_ped_merged_unrel)
dim(fvs_gen_ped_merged_unrel)
dim(vbi_gen_ped_merged_unrel)
dim(carl_gen_ped_merged_unrel)
# fve_gen_ped_merged_unrel <- subset(fve_gen_ped_merged, fvg_gen_ped_merged$ID1.x %in% fve_unrel_conv$CLIN | fvg_gen_ped_merged$ID2.x %in% fve_unrel_conv$CLIN)
# fvi_gen_ped_merged_unrel <- subset(fvi_gen_ped_merged, fvg_gen_ped_merged$ID1.x %in% fvi_unrel_conv$CLIN | fvg_gen_ped_merged$ID2.x %in% fvi_unrel_conv$CLIN)
# fvr_gen_ped_merged_unrel <- subset(fvr_gen_ped_merged, fvg_gen_ped_merged$ID1.x %in% fvr_unrel_conv$CLIN | fvg_gen_ped_merged$ID2.x %in% fvr_unrel_conv$CLIN)
# fvs_gen_ped_merged_unrel <- subset(fvs_gen_ped_merged, fvg_gen_ped_merged$ID1.x %in% fvs_unrel_conv$CLIN | fvg_gen_ped_merged$ID2.x %in% fvs_unrel_conv$CLIN)


# jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvg_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
# 	plot(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped,cex=0.5)
# 	points(fvg_gen_ped_merged_unrel$KIN.gen,fvg_gen_ped_merged_unrel$KIN.ped,col="red")
# dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/05162015/PLOTS/vbi_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped,cex=0.5)
	points(vbi_gen_ped_merged_unrel$KIN.gen,vbi_gen_ped_merged_unrel$KIN.ped,col="red")
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/05162015/PLOTS/carl_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped,cex=0.5)
	points(carl_gen_ped_merged_unrel$KIN.gen,carl_gen_ped_merged_unrel$KIN.ped,col="red")
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/05162015/PLOTS/fve_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fve_gen_ped_merged$KIN.gen,fve_gen_ped_merged$KIN.ped,cex=0.5)
	points(fve_gen_ped_merged_unrel$KIN.gen,fve_gen_ped_merged_unrel$KIN.ped,col="red")
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/05162015/PLOTS/fvi_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fvi_gen_ped_merged$KIN.gen,fvi_gen_ped_merged$KIN.ped,cex=0.5)
	points(fvi_gen_ped_merged_unrel$KIN.gen,fvi_gen_ped_merged_unrel$KIN.ped,col="red")
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/05162015/PLOTS/fvr_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fvr_gen_ped_merged$KIN.gen,fvr_gen_ped_merged$KIN.ped,cex=0.5)
	points(fvr_gen_ped_merged_unrel$KIN.gen,fvr_gen_ped_merged_unrel$KIN.ped,col="red")
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/05162015/PLOTS/fvs_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fvs_gen_ped_merged$KIN.gen,fvs_gen_ped_merged$KIN.ped,cex=0.5)
	points(fvs_gen_ped_merged_unrel$KIN.gen,fvs_gen_ped_merged_unrel$KIN.ped,col="red")
dev.off()

vbi_gen_ped_merged$cohort <- "VBI"
carl_gen_ped_merged$cohort <- "CARL"
fve_gen_ped_merged$cohort <- "FVG-E"
fvi_gen_ped_merged$cohort <- "FVG-I"
fvr_gen_ped_merged$cohort <- "FVG-R"
fvs_gen_ped_merged$cohort <- "FVG-S"


vbi_gen_ped_merged_unrel$cohort <- "VBI"
carl_gen_ped_merged_unrel$cohort <- "CARL"
fve_gen_ped_merged_unrel$cohort <- "FVG-E"
fvi_gen_ped_merged_unrel$cohort <- "FVG-I"
fvr_gen_ped_merged_unrel$cohort <- "FVG-R"
fvs_gen_ped_merged_unrel$cohort <- "FVG-S"

villages_n <- c("carl","vbi","fve","fvi","fvr","fvs")

for (village in villages_n){
	current_village_set <- paste(village,"_gen_ped_merged_unrel",sep="")
	# current_village_set <- paste(village,"_gen_ped_merged",sep="")
	current_village <- get(current_village_set)
	current_village <- current_village[colnames(vbi_gen_ped_merged_unrel)]
	current_village$KEY <- as.character(current_village$KEY)
	current_village$ID1.x <- as.character(current_village$ID1.x)
	current_village$ID1.y <- as.character(current_village$ID1.y)
	current_village$ID2.x <- as.character(current_village$ID2.x)
	current_village$ID2.y <- as.character(current_village$ID2.y)
	assign(current_village_set,current_village)
}


all_gen_ped_merged <- rbind(vbi_gen_ped_merged,carl_gen_ped_merged,fve_gen_ped_merged,fvi_gen_ped_merged,fvr_gen_ped_merged,fvs_gen_ped_merged)
all_gen_ped_merged_unrel <- rbind(vbi_gen_ped_merged_unrel,carl_gen_ped_merged_unrel,fve_gen_ped_merged_unrel,fvi_gen_ped_merged_unrel,fvr_gen_ped_merged_unrel,fvs_gen_ped_merged_unrel)

all_gen_ped_merged$rel <- as.factor("All")
all_gen_ped_merged_unrel$rel <- as.factor("Unrelated")

shapes <- as.data.frame(t(c(ALL="1",Unrel="4")))

names(shapes) <- c("ALL","Unrel")

all_gen_ped <- rbind(all_gen_ped_merged,all_gen_ped_merged_unrel)
all_gen_ped_nodups <- all_gen_ped[!duplicated(all_gen_ped$KEY,fromLast=T),]
all_gen_ped_nodups$cohort <- as.factor(all_gen_ped_nodups$cohort)
all_gen_ped_nodups$KIN.gen2 <- all_gen_ped_nodups$KIN.gen
all_gen_ped_nodups[all_gen_ped_nodups$KIN.gen <= 0.0441,"KIN.gen2"] <- 0

source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
pops_c <- c("CEU","TSI","CARL","VBI","FVG-E","FVG-I","FVG-R","FVG-S")
pop_colors <- col_pop(pops_c)
require(ggplot2)
require(reshape2)

pl <- ggplot(all_gen_ped_nodups)
pl <- pl + geom_point()
pl <- pl + aes(x = KIN.gen2, y = KIN.ped,colour=cohort,shape=rel)
pl <- pl + scale_shape(solid = FALSE)
pl <- pl + scale_colour_manual("Cohorts", values=pop_colors[as.character(unique(all_gen_ped_nodups$cohort))])
pl <- pl + xlab("Genetic Kinship")
pl <- pl + ylab("Pedigree Kinship")
pl <- pl + facet_wrap(~cohort, ncol=2)
# pl <- pl + ggtitle(main)
pl <- pl + theme_bw()
pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
# pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))
   
base_plot <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/05162015/PLOTS"
# ggsave(filename=paste(base_plot,"/kinship_corr_all.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
ggsave(filename=paste(base_plot,"/kinship_corr_all_gen2.jpeg",sep=""),width=12, height=7,dpi=300,plot=pl)
