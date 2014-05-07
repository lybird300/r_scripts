#r script to plot for enza

#plot 1 DAF spectrun
#upload data for each population:

chr <- "22"
pops <- c("CEU","TSI","VBI","FVG")

for (pop in pops) {

pop_table_file <- paste("/lustre/scratch113/projects/esgi-vbseq/20140430_purging/TABLES/",pop,".chr",chr,".tab",sep="")

pop_table <- read.table(pop_table_file,header=F)
colnames(pop_table) <- c("CHROM","POZ","POS","ID","REF","ALT","INFO","REC","ALC","DAC","MAC","DAF","MAF")
pop_table$DAF <- as.numeric(as.character(pop_table$DAF))

jpeg(paste(pop,"DAF.jpg",sep="_"),width=800, height=800,pointsize = 20)
  hist(pop_table$DAF,main=paste("DAF in ",pop,sep=""),xlab="DAF")
dev.off()
}