#Calculate correlation kinship gen-ped
rm(list=ls())
library(kinship2)

# fvr_ped <- read.table('/home/max/Work/FVG/pedigree_FVG/SEQ/resia_ped.csv.seq.list')
fvr_ped <- read.table('/home/max/Work/FVG/pedigree_FVG/VILLAGES/resia_ped_sorted_manmod.csv')
#had to manually fix the pedigree data
fvr_pedigree <- pedigree(id=fvr_ped$V1,dadid=fvr_ped$V2,momid=fvr_ped$V3,sex=fvr_ped$V4,miss=0)
fvr_kinship <- kinship(fvr_pedigree)

fve_ped <- read.table('/home/max/Work/FVG/pedigree_FVG/VILLAGES/erto_ped.csv',skip=1)
fve_pedigree <- pedigree(id=fve_ped$V1,dadid=fve_ped$V2,momid=fve_ped$V3,sex=fve_ped$V4,miss=0)
#I had to manually remove mothers from DUOS relationships (id with no father but mother info)
fve_kinship <- kinship(fve_pedigree)

fvi_ped <- read.table('/home/max/Work/FVG/pedigree_FVG/VILLAGES/illegio_ped_uniq.csv')
fvi_pedigree <- pedigree(id=fvi_ped$V1,dadid=fvi_ped$V2,momid=fvi_ped$V3,sex=fvi_ped$V4,miss=0)
#added manually founders
fvi_kinship <- kinship(fvi_pedigree)

fvs_ped <- read.table('/home/max/Work/FVG/pedigree_FVG/VILLAGES/sauris_ped_uniq.csv')
fvs_pedigree <- pedigree(id=fvs_ped$V1,dadid=fvs_ped$V2,momid=fvs_ped$V3,sex=fvs_ped$V4,miss=0)
#removed duplicates founders
fvs_kinship <- kinship(fvs_pedigree)

#now load genomic kinship done with king
fvg_seq_list <- read.table('/home/max/Work/FVG/SEQ_CALLING/FVG_seq_sorted.list',h=F)
colnames(fvg_seq_list) <- c("ID")
fvg_seq_list$ID <- as.character(fvg_seq_list$ID)

#split kinships in villages extracting only sequenced samples
fvr_kinship_seq <- fvr_kinship[colnames(fvr_kinship) %in% fvg_seq_list$ID, rownames(fvr_kinship) %in% fvg_seq_list$ID]
fve_kinship_seq <- fve_kinship[colnames(fve_kinship) %in% fvg_seq_list$ID, rownames(fve_kinship) %in% fvg_seq_list$ID]
fvi_kinship_seq <- fvi_kinship[colnames(fvi_kinship) %in% fvg_seq_list$ID, rownames(fvi_kinship) %in% fvg_seq_list$ID]
fvs_kinship_seq <- fvs_kinship[colnames(fvs_kinship) %in% fvg_seq_list$ID, rownames(fvs_kinship) %in% fvg_seq_list$ID]
# dim(fvr_kinship_seq)
# dim(fve_kinship_seq)
# dim(fvi_kinship_seq)
# dim(fvs_kinship_seq)

# dim(fvr_kinship)
# dim(fve_kinship)
# dim(fvi_kinship)
# dim(fvs_kinship)


#load genomic kinship previously calculated
base_folder <- "/home/max/Work/Analyses/PURGE_INBREEDING/KINSHIP"
# base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING"

fve_kinship_king_seq <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Erto.keeplist.ibs0.kinship.conv.kin",sep=""),header=F)
fvi_kinship_king_seq <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Illegio.keeplist.ibs0.kinship.conv.kin",sep=""),header=F)
fvr_kinship_king_seq <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Resia.keeplist.ibs0.kinship.conv.kin",sep=""),header=F)
fvs_kinship_king_seq <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Sauris.keeplist.ibs0.kinship.conv.kin",sep=""),header=F)


dim(fve_kinship_king_seq)
dim(fvi_kinship_king_seq)
dim(fvr_kinship_king_seq)
dim(fvs_kinship_king_seq)

#load conversion table data
conv_fvg <- read.table("/home/max/Work/Analyses/PURGE_INBREEDING/WGS_FVG_id_conversion_complete_sorted.list")
colnames(conv_fvg) <- c("SEQ","CLIN")
conv_fvg$SEQ <- as.character(conv_fvg$SEQ)

#clean tables from useless columns
kin_set <- c("fve_kinship_king_seq","fvi_kinship_king_seq","fvr_kinship_king_seq","fvs_kinship_king_seq")
current_set <- NULL
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
	current_merged <- merge(current_set,conv_fvg,by.x="ID1",by.y="SEQ",all.x)
	current_merged <- merge(current_merged,conv_fvg,by.x="ID2",by.y="SEQ",all.x)
	assign(set,current_merged)
}

library(reshape2)
kin_ped_set <- c("fve_kinship_seq","fvi_kinship_seq","fvr_kinship_seq","fvs_kinship_seq","fve_kinship","fvi_kinship","fvr_kinship","fvs_kinship")
# kin_ped_set <- c("fvr_kinship")
current_kin_ped_set <- NULL
for (ped_set in kin_ped_set){
	current_ped_set <- get(ped_set)
	#convert ids
	current_ped_set_melted <- melt(current_ped_set, varnames = c('ID1', 'ID2'), na.rm = TRUE)
	# current_ped_set_melted <- current_ped_set_melted[which(current_ped_set_melted$ID1 != current_ped_set_melted$ID2),]
	# current_ped_set_melted <- current_ped_set_melted[which(paste(current_ped_set_melted$ID1,current_ped_set_melted$ID2,sep="_") != paste(current_ped_set_melted$ID2,current_ped_set_melted$ID1,sep="_")),]
	assign(paste(ped_set,"_melted",sep=""),current_ped_set_melted)
}

kin_ped_set_melted <- c("fve_kinship_seq_melted","fvi_kinship_seq_melted","fvr_kinship_seq_melted","fvs_kinship_seq_melted","fve_kinship_melted","fvi_kinship_melted","fvr_kinship_melted","fvs_kinship_melted")
for (ped_set in kin_ped_set_melted){
	current_ped_set <- get(ped_set)
	#write tables
	write.table(current_ped_set,file=paste(getwd(),ped_set,sep="/"),quote=F,sep="\t",row.names=F,col.names=F)
}


dim(fve_kinship_king_seq)
dim(fve_kinship_seq_melted)
fve_kinship_seq_melted$key <- paste(apply(fve_kinship_seq_melted[,1:2],1,min),apply(fve_kinship_seq_melted[,1:2],1,max),sep="_")
fve_kinship_seq_melted <- fve_kinship_seq_melted[!duplicated(fve_kinship_seq_melted$key),]
dim(fve_kinship_seq_melted)

dim(fvi_kinship_king_seq)
dim(fvi_kinship_seq_melted)
fvi_kinship_seq_melted$key <- paste(apply(fvi_kinship_seq_melted[,1:2],1,min),apply(fvi_kinship_seq_melted[,1:2],1,max),sep="_")
fvi_kinship_seq_melted <- fvi_kinship_seq_melted[!duplicated(fvi_kinship_seq_melted$key),]
dim(fvi_kinship_seq_melted)

dim(fvr_kinship_king_seq)
dim(fvr_kinship_seq_melted)
fvr_kinship_seq_melted$key <- paste(apply(fvr_kinship_seq_melted[,1:2],1,min),apply(fvr_kinship_seq_melted[,1:2],1,max),sep="_")
fvr_kinship_seq_melted <- fvr_kinship_seq_melted[!duplicated(fvr_kinship_seq_melted$key),]
dim(fvr_kinship_seq_melted)

dim(fvs_kinship_king_seq)
dim(fvs_kinship_seq_melted)
fvs_kinship_seq_melted$key <- paste(apply(fvs_kinship_seq_melted[,1:2],1,min),apply(fvs_kinship_seq_melted[,1:2],1,max),sep="_")
fvs_kinship_seq_melted <- fvs_kinship_seq_melted[!duplicated(fvs_kinship_seq_melted$key),]
dim(fvs_kinship_seq_melted)

length(unique(fve_kinship_seq_melted$ID1))
length(unique(fve_kinship_king_seq$CLIN.x))

length(unique(fvi_kinship_seq_melted$ID1))
length(unique(fvi_kinship_king_seq$CLIN.x))

length(unique(fvr_kinship_seq_melted$ID1))
length(unique(fvr_kinship_king_seq$CLIN.x))

length(unique(fvs_kinship_seq_melted$ID1))
length(unique(fvs_kinship_king_seq$CLIN.x))

fve_kinship_king_seq2 <- (fve_kinship_king_seq[which(fve_kinship_king_seq$CLIN.x %in% fve_kinship_seq_melted$ID1),][which(fve_kinship_king_seq$CLIN.y %in% fve_kinship_seq_melted$ID2),])
fvi_kinship_king_seq2 <- (fvi_kinship_king_seq[which(fvi_kinship_king_seq$CLIN.x %in% fvi_kinship_seq_melted$ID1),][which(fvi_kinship_king_seq$CLIN.y %in% fvi_kinship_seq_melted$ID2),])
fvr_kinship_king_seq2 <- (fvr_kinship_king_seq[which(fvr_kinship_king_seq$CLIN.x %in% fvr_kinship_seq_melted$ID1),][which(fvr_kinship_king_seq$CLIN.y %in% fvr_kinship_seq_melted$ID2),])
fvs_kinship_king_seq2 <- (fvs_kinship_king_seq[which(fvs_kinship_king_seq$CLIN.x %in% fvs_kinship_seq_melted$ID1),][which(fvs_kinship_king_seq$CLIN.y %in% fvs_kinship_seq_melted$ID2),])

length(unique(fve_kinship_king_seq2$CLIN.x))
length(unique(fve_kinship_seq_melted$ID1))

length(unique(fvi_kinship_king_seq2$CLIN.x))
length(unique(fvi_kinship_seq_melted$ID1))

length(unique(fvr_kinship_king_seq2$CLIN.x))
length(unique(fvr_kinship_seq_melted$ID1))

length(unique(fvs_kinship_king_seq2$CLIN.x))
length(unique(fvs_kinship_seq_melted$ID1))


fve_kinship_king_seq2$key <- paste(apply(fve_kinship_king_seq2[,4:5],1,min),apply(fve_kinship_king_seq2[,4:5],1,max),sep="_")
fve_kinship_king_seq3 <- fve_kinship_king_seq2[!duplicated(fve_kinship_king_seq2$key),]

fvi_kinship_king_seq2$key <- paste(apply(fvi_kinship_king_seq2[,4:5],1,min),apply(fvi_kinship_king_seq2[,4:5],1,max),sep="_")
fvi_kinship_king_seq3 <- fvi_kinship_king_seq2[!duplicated(fvi_kinship_king_seq2$key),]

fvr_kinship_king_seq2$key <- paste(apply(fvr_kinship_king_seq2[,4:5],1,min),apply(fvr_kinship_king_seq2[,4:5],1,max),sep="_")
fvr_kinship_king_seq3 <- fvr_kinship_king_seq2[!duplicated(fvr_kinship_king_seq2$key),]

fvs_kinship_king_seq2$key <- paste(apply(fvs_kinship_king_seq2[,4:5],1,min),apply(fvs_kinship_king_seq2[,4:5],1,max),sep="_")
fvs_kinship_king_seq3 <- fvs_kinship_king_seq2[!duplicated(fvs_kinship_king_seq2$key),]


length(unique(fve_kinship_king_seq3$CLIN.x))
length(unique(fvi_kinship_king_seq3$CLIN.x))
length(unique(fvr_kinship_king_seq3$CLIN.x))
length(unique(fvs_kinship_king_seq3$CLIN.x))

length(unique(fve_kinship_seq_melted$ID1))
length(unique(fvi_kinship_seq_melted$ID1))
length(unique(fvr_kinship_seq_melted$ID1))
length(unique(fvs_kinship_seq_melted$ID1))

#sort dataframes
fve_kinship_king_seq3 <- fve_kinship_king_seq3[with(fve_kinship_king_seq3,order(key)),]
fve_kinship_seq_melted <- fve_kinship_seq_melted[with(fve_kinship_seq_melted,order(key)),]
fvi_kinship_king_seq3 <- fvi_kinship_king_seq3[with(fvi_kinship_king_seq3,order(key)),]
fvi_kinship_seq_melted <- fvi_kinship_seq_melted[with(fvi_kinship_seq_melted,order(key)),]
fvr_kinship_king_seq3 <- fvr_kinship_king_seq3[with(fvr_kinship_king_seq3,order(key)),]
fvr_kinship_seq_melted <- fvr_kinship_seq_melted[with(fvr_kinship_seq_melted,order(key)),]
fvs_kinship_king_seq3 <- fvs_kinship_king_seq3[with(fvs_kinship_king_seq3,order(key)),]
fvs_kinship_seq_melted <- fvs_kinship_seq_melted[with(fvs_kinship_seq_melted,order(key)),]

fve_kinship_seq_melted <- fve_kinship_seq_melted[which(fve_kinship_seq_melted$ID1 != fve_kinship_seq_melted$ID2),]
fvi_kinship_seq_melted <- fvi_kinship_seq_melted[which(fvi_kinship_seq_melted$ID1 != fvi_kinship_seq_melted$ID2),]
fvr_kinship_seq_melted <- fvr_kinship_seq_melted[which(fvr_kinship_seq_melted$ID1 != fvr_kinship_seq_melted$ID2),]
fvs_kinship_seq_melted <- fvs_kinship_seq_melted[which(fvs_kinship_seq_melted$ID1 != fvs_kinship_seq_melted$ID2),]

(length(fve_kinship_king_seq3$key))
(length(fve_kinship_seq_melted$key))
(length(fvi_kinship_king_seq3$key))
(length(fvi_kinship_seq_melted$key))
(length(fvr_kinship_king_seq3$key))
(length(fvr_kinship_seq_melted$key))
(length(fvs_kinship_king_seq3$key))
(length(fvs_kinship_seq_melted$key))


#####################################################################################
#Correlation between kinship

fve_kinship_gen <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Erto.keeplist.ibs0.kinship.conv.kin.over_ped",sep=""),header=F)
fvi_kinship_gen <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Illegio.keeplist.ibs0.kinship.conv.kin.over_ped",sep=""),header=F)
fvr_kinship_gen <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Resia.keeplist.ibs0.kinship.conv.kin.over_ped",sep=""),header=F)
fvs_kinship_gen <- read.table(paste(base_folder,"/KINSHIP/KIN_GEN/FVG_Sauris.keeplist.ibs0.kinship.conv.kin.over_ped",sep=""),header=F)

fve_kinship_ped <- read.table(paste(base_folder,"/KINSHIP/KIN_PED/fve_kinship_melted.sorted.pkin.over_gen",sep=""),header=F)
fvi_kinship_ped <- read.table(paste(base_folder,"/KINSHIP/KIN_PED/fvi_kinship_melted.sorted.pkin.over_gen",sep=""),header=F)
fvr_kinship_ped <- read.table(paste(base_folder,"/KINSHIP/KIN_PED/fvr_kinship_melted.sorted.pkin.over_gen",sep=""),header=F)
fvs_kinship_ped <- read.table(paste(base_folder,"/KINSHIP/KIN_PED/fvs_kinship_melted.sorted.pkin.over_gen",sep=""),header=F)

colnames(fve_kinship_gen) <- c("ID1","ID2","KIN","KEY")
colnames(fve_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fve_gen_ped_merged <- merge(fve_kinship_gen,fve_kinship_ped,by.x="KEY",by.y="KEY",all.x)


colnames(fvi_kinship_gen) <- c("ID1","ID2","KIN","KEY")
colnames(fvi_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fvi_gen_ped_merged <- merge(fvi_kinship_gen,fvi_kinship_ped,by.x="KEY",by.y="KEY",all.x)

colnames(fvr_kinship_gen) <- c("ID1","ID2","KIN","KEY")
colnames(fvr_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fvr_gen_ped_merged <- merge(fvr_kinship_gen,fvr_kinship_ped,by.x="KEY",by.y="KEY",all.x)

colnames(fvs_kinship_gen) <- c("ID1","ID2","KIN","KEY")
colnames(fvs_kinship_ped) <- c("ID1","ID2","KIN","KEY")
fvs_gen_ped_merged <- merge(fvs_kinship_gen,fvs_kinship_ped,by.x="KEY",by.y="KEY",all.x)

#conversion table kinship2 / KING
# MZT: 0.5 > 0.3535
# 1st degree: 0.25 (0.1767,0.3535)
# 2nd degree: 0.125 (0.0883,0.1767)
# 3rd degree: 0.0625 (0.0441,0.0883)
# Unrelated: 0 < 0.0441

villages <- c("fve","fvi","fvr","fvs","vbi","carl","fvg")

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

# Kinship correlation values from cor.test(),Pearson's coefficient
# Fvg-E: r= 0.6517779, pval = <2.2e-16
# Fvg-I: r= 0.775517, pval = <2.2e-16
# Fvg-R: r= NA, pval = NA
# Fvg-S: r= 0.5232983, pval = <2.2e-16

# > cor_fve

# 	Pearson's product-moment correlation
# data:  fve_gen_ped_merged$KIN.x and fve_gen_ped_merged$KIN.y
# t = 32.4872, df = 1429, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.6260273 1.0000000
# sample estimates:
#       cor 
# 0.6517779 
# > cor_fvi
# 	Pearson's product-moment correlation

# data:  fvi_gen_ped_merged$KIN.x and fvi_gen_ped_merged$KIN.y
# t = 47.3051, df = 1483, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.7579145 1.0000000
# sample estimates:
#      cor 
# 0.775517 
# > cor_fvr
# 	Pearson's product-moment correlation
# data:  fvr_gen_ped_merged$KIN.x and fvr_gen_ped_merged$KIN.y
# t = NA, df = 2483, p-value = NA
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  NA  1
# sample estimates:
# cor 
#  NA 
# > cor_fvs
# 	Pearson's product-moment correlation
# data:  fvs_gen_ped_merged$KIN.x and fvs_gen_ped_merged$KIN.y
# t = 22.7794, df = 1376, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.4903437 1.0000000
# sample estimates:
#       cor 
# 0.5232983 

cor_fve <- plot(fve_gen_ped_merged$KIN.gen,fve_gen_ped_merged$KIN.ped)
cor_fvi <- plot(fvi_gen_ped_merged$KIN.gen,fvi_gen_ped_merged$KIN.ped)
cor_fvr <- plot(fvr_gen_ped_merged$KIN.gen,fvr_gen_ped_merged$KIN.ped)
cor_fvs <- plot(fvs_gen_ped_merged$KIN.gen,fvs_gen_ped_merged$KIN.ped)


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


vbi_kinship_ped <- read.table("/home/max/Work/Analyses/PURGE_INBREEDING/KINSHIP/KIN_PED/vbi_kinship_melted.sorted.pkin.over_gen",header=F)
vbi_kinship_gen <- read.table("/home/max/Work/Analyses/PURGE_INBREEDING/KINSHIP/KIN_GEN/VBI.keeplist.ibs0.kinship.conv.kin.over_ped",header=F)

colnames(vbi_kinship_gen) <- c("ID1","ID2","KIN","KEY")
colnames(vbi_kinship_ped) <- c("ID1","ID2","KIN","KEY")
vbi_gen_ped_merged <- merge(vbi_kinship_gen,vbi_kinship_ped,by.x="KEY",by.y="KEY",all.x)


cor_vbi <- cor.test(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped,alternative="g")
cor.test(vbi_gen_ped_merged$KIN.cast,vbi_gen_ped_merged$KIN.ped,alternative="g")
corrgram(vbi_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)


# Kinship CARL
carl_kinship_ped <- read.table("/home/max/Work/CARLANTINO/pedigree_CARL/carlantino_sorted_nofam.ped",header=F)
carl_pedigree <- pedigree(id=carl_kinship_ped$V1,dadid=carl_kinship_ped$V2,momid=carl_kinship_ped$V3,sex=carl_kinship_ped$V4,miss=0)
carl_kinship <- kinship(carl_pedigree)

library(reshape2)
carl_kinship_melted <- melt(carl_kinship, varnames = c('ID1', 'ID2'), na.rm = TRUE)

write.table(carl_kinship_melted,file=paste(getwd(),"carl_kinship_melted",sep="/"),quote=F,sep="\t",row.names=F,col.names=F)

carl_kinship_ped <- read.table("/home/max/Work/Analyses/PURGE_INBREEDING/KINSHIP/KIN_PED/carl_kinship_melted.sorted.pkin.over_gen",header=F)
carl_kinship_gen <- read.table("/home/max/Work/Analyses/PURGE_INBREEDING/KINSHIP/KIN_GEN/CARL.keeplist.ibs0.kinship.conv.kin.over_ped",header=F)

colnames(carl_kinship_gen) <- c("ID1","ID2","KIN","KEY")
colnames(carl_kinship_ped) <- c("ID1","ID2","KIN","KEY")
carl_gen_ped_merged <- merge(carl_kinship_gen,carl_kinship_ped,by.x="KEY",by.y="KEY",all.x)
cor_carl <- cor.test(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped,alternative="g")
cor.test(carl_gen_ped_merged$KIN.cast,carl_gen_ped_merged$KIN.ped,alternative="g")
corrgram(carl_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)

rcor_fvg <- rcorr(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped)
rcor_vbi <- rcorr(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped)
rcor_carl <- rcorr(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped)

cov_fvg <- cov(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped)
cov_vbi <- cov(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped)
cov_carl <- cov(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped)

# Kinship correlation values from cor.test(),Pearson's coefficient
# Fvg: r= 0.611911,pval = <2.2e-16
# Vbi: r= 0.5021468,pval = <2.2e-16
# Carl: r= 0.8717032,pval = <2.2e-16

# > cor_fvg

# 	Pearson's product-moment correlation

# data:  fvg_gen_ped_merged$KIN.x and fvg_gen_ped_merged$KIN.y
# t = 63.6899, df = 6777, p-value < 2.2e-16
# alternative hypothesis: true correlation is greater than 0
# 95 percent confidence interval:
#  0.5992579 1.0000000
# sample estimates:
#      cor 
# 0.611911 

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
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvg_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
	corrgram(fvg_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/vbi_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
	corrgram(vbi_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/carl_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
	corrgram(carl_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fve_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
	corrgram(fve_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvi_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
	corrgram(fvi_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvr_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
	corrgram(fvr_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvs_gen_ped_merged.jpg",width=800, height=800,pointsize = 20)
	corrgram(fvs_gen_ped_merged[,c(4,7,8)],order=TRUE, lower.panel=panel.shade,upper.panel=panel.pts, text.panel=panel.txt)
dev.off()


jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvg_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fvg_gen_ped_merged$KIN.gen,fvg_gen_ped_merged$KIN.ped)
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/vbi_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(vbi_gen_ped_merged$KIN.gen,vbi_gen_ped_merged$KIN.ped)
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/carl_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(carl_gen_ped_merged$KIN.gen,carl_gen_ped_merged$KIN.ped)
dev.off()

jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fve_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fve_gen_ped_merged$KIN.gen,fve_gen_ped_merged$KIN.ped)
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvi_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fvi_gen_ped_merged$KIN.gen,fvi_gen_ped_merged$KIN.ped)
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvr_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fvr_gen_ped_merged$KIN.gen,fvr_gen_ped_merged$KIN.ped)
dev.off()
jpeg("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/PURGE_INBREEDING/KINSHIP/fvs_gen_ped_merged_ncorr.jpg",width=800, height=800,pointsize = 20)
	plot(fvs_gen_ped_merged$KIN.gen,fvs_gen_ped_merged$KIN.ped)
dev.off()
