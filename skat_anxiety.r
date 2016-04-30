#script for analyses on seq data using SKAT

#get boundaries for genes from bioMart:

library(biomaRt)
listMarts(host="www.ensembl.org")
ensemble <- useEnsembl("ensembl", GRCh=37)
#hgnc <- useEnsembl("HGNC", GRCh=37)
ensemble <- useDataset("hsapiens_gene_ensembl",mart=ensemble)

attributes <- listAttributes(ensemble)
filters <- listFilters(ensemble)

attributes_to_get <- c("chromosome_name","start_position","end_position","hgnc_symbol")
filters_to_use <- c("hgnc_symbol")
base_folder <- "/home/cocca/analyses/AnxietySheila/LISTS/"

genes <- read.table(paste(base_folder,"gene_names_list.txt",sep=""))
gene_ids <- as.character(genes$V1)

gene_list <- getBM(attributes=attributes_to_get,filters=filters_to_use,values=gene_ids,mart=ensemble)

write.table(gene_list,file=paste(base_folder,'gene_list_anxiety.list',sep=""),row.names=F,quote=F,col.names=T,sep="\t")

#now we need to read a pheno file, with samples and so on, in skat format
pheno_list <- read.table(paste(base_folder,"panico_ansia_list.phen",sep=""),header=T)

#intersect this list with the sequenced samples list
wgs_samples <- read.table(paste(base_folder,"cases_controls_sample_WGS_exact.list",sep=""),header=F)
colnames(wgs_samples) <- c("WGS_ID","CLINIC_ID")

wgs_samples_pheno <- merge(wgs_samples,pheno_list,by.x="CLINIC_ID",by.y="CODPAZ")

#now we will run analyses on:
# anx: ansia
# mdd: depressione 
# ap: panico
#using the gene list

# would be good to add sex and age as covariate, in a second round



