#r script to plot for enza
rm(list=ls())

base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/"
base_folder <- "/lustre/scratch113/projects/esgi-vbseq/20140430_purging/UNRELATED/FIVE_POPS"

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
