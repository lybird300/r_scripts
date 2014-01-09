#The matching list MUST BE in the format: NEW_ID OLD_ID

#region_up <- read.table('region_file_chr4_up.txt', header=F,sep='\t',col.names=c("chr","start","end","name"))
#region_down <- read.table('region_file_chr4_down.txt', header=F,sep='\t',col.names=c("chr","start","end","name"))

region <- read.table('region_file_chr12.txt', header=F,sep='\t',col.names=c("chr","start","end","name"))
new_region <- NULL

#FIXME:take in account the fact that the file has odd or even lines, so we want to use a step factor
# that let us use all lines

for(j in seq(10,length(region$chr), by=10)){
	#print(i)
	i=j-9
	#print(j)
	new_region_line <- cbind(as.character(region[i,]$chr),as.character(region[i,]$start),as.character(region[j,]$start),paste(as.character(region[i,]$name),as.character(region[j,]$name),sep="_"))
	new_region <-rbind(new_region,new_region_line)
}

colnames(new_region) <- c("chr","start","end","name")

write.table(new_region,file='region_file.bed',sep="\t",quote=F,row.names=F,col.names=F,append=T)

