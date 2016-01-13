###############################################################################################################
##### 22/09/2015 
##### calculate means Fst pairwise
###### PLOT with ggplot a heaatmap
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")
require(plotrix)
require(ggplot2)
require(reshape2)
base_folder <- getwd()
#WE DONT't have CARL pop anymore
pops <- c("CEU","TSI","VBI","FVG-E","FVG-I","FVG-R","FVG-S")
pops_for_tables <- c("CEU","TSI","VBI","Erto","Illegio","Resia","Sauris")
pops_p <- c("CEU_p","TSI_p","VBI_p","FVG-E_p","FVG-I_p","FVG-R_p","FVG-S_p")
pops_s <- c("CEU_s","TSI_s","VBI_s","FVG-E_s","FVG-I_s","FVG-R_s","FVG-S_s")
pops_n <- c("CEU","TSI","VBI","FVG_E","FVG_I","FVG_R","FVG_S")
pops_n_ingi <- c("VBI","FVG_E","FVG_I","FVG_R","FVG_S")
pops_ingi_novel <- c("VBI_n","FVG_n")
pops_ingi_class <- c("VBI_p","FVG_p","VBI_s","FVG_s")

hwe_mean_all <- NULL
for (pop1 in pops_for_tables) {
			current_hwe_file <- paste(base_folder,"/",pop1,"_hwe.hwe.clean.gz",sep="")
			current_hwe_table <- read.table(current_hwe_file,header=T)
			current_hwe_summary <- summary(current_hwe_table[,3])
			assign(paste(pop1,"hwe_summary",sep="_"),current_hwe_summary)
			assign(paste(pop1,"hwe",sep="_"),current_hwe_table)
}

all_hwe_info <- NULL

for (pop in pops_for_tables) {
	# current_hwe <- get(paste(pop,"hwe"))
	# assign(paste(pop,"hwe",sep="_"),current_hwe)	
	current_hwe <- get(paste(pop,"hwe",sep="_"))
	save(current_hwe,file=paste(pop,"hwe.RData",sep="_"))
	# current_hwe_summary <- summary(current_hwe[,3])
	# current_hwe_stde <- sd(current_hwe[,3])/sqrt(length(current_hwe[,3]))
	# current_hwe_info <- as.data.frame(t(c(cohort=pop,current_hwe_summary,sd=sd(current_hwe[,3]),stderr=current_hwe_stde)))
	# all_hwe_info <- rbind(all_hwe_info,current_hwe_info)
}
	
# save(fst_mean_all,file="fst_mean_all.RData")
load("fst_mean_all.RData")
#plot a heatmap with ggplot
fst_mean_all.m <- acast(fst_mean_all,pop1~pop2,value.var="Fst")
colnames(fst_mean_all) <- c("pop1","pop2","Fst")
fst_mean_all$pop1 <- as.character(fst_mean_all$pop1)
fst_mean_all$pop2 <- as.character(fst_mean_all$pop2)
fst_mean_all[grep("Erto",fst_mean_all$pop1),]$pop1 <- "FVG-E"
fst_mean_all[grep("Erto",fst_mean_all$pop2),]$pop2 <- "FVG-E"
fst_mean_all[grep("Illegio",fst_mean_all$pop1),]$pop1 <- "FVG-I"
fst_mean_all[grep("Illegio",fst_mean_all$pop2),]$pop2 <- "FVG-I"
fst_mean_all[grep("Resia",fst_mean_all$pop1),]$pop1 <- "FVG-R"
fst_mean_all[grep("Resia",fst_mean_all$pop2),]$pop2 <- "FVG-R"
fst_mean_all[grep("Sauris",fst_mean_all$pop1),]$pop1 <- "FVG-S"
fst_mean_all[grep("Sauris",fst_mean_all$pop2),]$pop2 <- "FVG-S"
fst_mean_all$pop1 <- as.factor(fst_mean_all$pop1)
fst_mean_all$pop2 <- as.factor(fst_mean_all$pop2)

# fst_mean_all$pop1 <- factor(fst_mean_all$pop1,level=names(pops))
# fst_mean_all$pop2 <- factor(fst_mean_all$pop2,level=names(pops))

pl <- ggplot(fst_mean_all, aes(pop1, pop2))
pl <- pl + geom_tile(aes(fill = Fst), colour = "white")
pl <- pl + scale_fill_gradient(low = "lightgray", high = "steelblue")

pl <- pl + xlab("")
pl <- pl + ylab("")
pl <- pl + guides(colour = guide_legend(override.aes = list(shape = 2)))
pl <- pl + theme_bw()
pl <- pl + theme(strip.text.x = element_text(size = 20))
pl <- pl + theme(axis.text.x=element_text(size = rel(1.2)))
pl <- pl + theme(axis.text.y=element_text(size = rel(1.2)))
pl <- pl + theme(axis.title= element_text(size=rel(1.2)))
pl <- pl + theme(legend.text= element_text(size = rel(1.2)), legend.title = element_text(size = rel(1.2)))

ggsave(filename=paste(base_folder,"/1_fst.jpeg",sep=""),width=8, height=4,units="in",dpi=300,plot=pl)
# png(paste(base_folder,"/1_fst.png",sep=""),width=700, height=200)
# print(pl)
# dev.off()

