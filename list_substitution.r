#The matching list MUST BE in the format: NEW_ID OLD_ID

pheno <- read.table('PLT_trait.txt', header=TRUE,sep='\t')
pheno$id <- as.character(pheno$id)
pheno[1:10,]
pheno$id

matching <- read.table('../sanger_seq_list_109.txt', header=FALSE,sep='\t',col.names=c("new","old"))
matching$old <- as.character(matching$old)
matching$new <- as.character(matching$new)
matching[1:10,]

for(i in 1:length(matching$old)){
	print(i)
	print(matching[i,]$new)
	pheno$id <- gsub(as.character(matching[i,]$old),matching[i,]$new,pheno$id)
}

write.table(pheno,file='PLT_trait.txt.tmp',sep="\t",quote=F,row.names=F,col.names=F)

system("awk \'BEGIN{FS=\"\t\";OFS=\"\t\"}{print $0,\"n/a\"}\' PLT_trait.txt.tmp > PLT_trait_seq.txt")
