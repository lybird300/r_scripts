#calculate corrected p-value
pvalcheck <- function(current_table,gz_flag,trait,chrom=NULL){

	if(gz_flag == T){
		mycon<- gzcon(gzfile(current_table, open="r"))
			# pvals<-read.table(textConnection(readLines(mycon)),header=T,comment.char="",colClasses=c("numeric","character","numeric","numeric","numeric"))
			# pvals<-read.table(textConnection(readLines(mycon)),header=T,colClasses=c("numeric","character","numeric","NULL","numeric","numeric","NULL","NULL","numeric","NULL","NULL"),stringsAsFactors=F,comment.char="")
			pvals<-read.table(textConnection(readLines(mycon)),header=T,colClasses=c("character","numeric"),stringsAsFactors=F,comment.char="")
			#pvals<-read.table(textConnection(readLines(mycon)),header=T,colClasses=c("numeric","character","numeric","numeric","numeric","numeric","numeric"),stringsAsFactors=F,comment.char="")
			# pvals<-read.table(textConnection(readLines(mycon)),header=T,comment.char="",colClasses=c("numeric","numeric","numeric"))
		close(mycon)

	}else{
		# pvals<-read.table(current_table,header=T,comment.char="",colClasses=c("numeric","character","numeric","numeric","numeric"))
		#coltypes<-read.table(current_table,header=T,colClasses=c("numeric","character","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"),stringsAsFactors=F,comment.char="",nrows=5)
		coltypes<-read.table(current_table,header=T,colClasses=c("character","numeric"),stringsAsFactors=F,comment.char="",nrows=5)
		classes <- as.list(sapply(coltypes, class))
		pvals <- as.data.frame(scan(current_table,what=as.list(classes),sep="\t",skip=1),stringsAsFactors=F,comment.char="")
		# pvals<-read.table(current_table,header=T,comment.char="",colClasses=c("numeric","numeric","numeric"))
	}

	#rename key columns for plotting
	# pvals$n_miss <- NULL
	# pvals$l_remle <- NULL
	# pvals$l_mle <- NULL
	# pvals$p_lrt <- NULL
	# pvals$p_score <- NULL
	#colnames(pvals) <- c("CHR","SNP","BP","beta","se","P","CHI2")
	colnames(pvals) <- c("CHR","CHI2")
	#pvals$CHR <- as.numeric(as.character(pvals$CHR))
	#pvals$BP <- as.numeric(as.character(pvals$BP))
	#pvals$beta <- as.numeric(as.character(pvals$beta))
	#pvals$se <- as.numeric(as.character(pvals$se))
	#pvals$P <- as.numeric(as.character(pvals$P))
	#pvals$CHI2 <- as.numeric(as.character(pvals$CHI2))
	#calc chi2 
	# chi2 <- (pvals$beta/pvals$se)^2
	# calc lambda 
	# qchisq(0.5,1) = 0.4549364
	lambda <- median(pvals$CHI2,na.rm=T)/qchisq(0.5,1)
	lambda
	# # calc pvalue 
	# p <- pchisq(chi2, df=1, ncp=0, lower.tail = F, log.p = FALSE)
	# calc pgc if lamdba > 1 
	if(lambda <= 1){ 
		#pgc <- pvals$P
	}else{ 
		chi2_gc <- (pvals$CHI2/lambda)
		pgc <- pchisq(chi2_gc, df=1, ncp=0, lower.tail = F, log.p = FALSE) 
	}

	pvals_new <- as.data.frame(cbind(pvals,PGC=pgc))

	#Create a manhattan plot for non corrrected pvalues
	if(!is.null(chrom)){
		man_name <- paste(trait, "_", chrom, ".manhattan.unfiltered.jpg", sep="")
		qq_name <- paste(trait, "_", chrom, ".qq.unfiltered.jpg", sep="")
		plot_main <- paste(trait, "_", chrom, " unfiltered", sep="")

		#create a plot for corrected pval
		man_name_c <- paste(trait, "_", chrom, ".manhattan.pgc.jpg", sep="")
		qq_name_c <- paste(trait, "_", chrom, ".qq.pgc.jpg", sep="")
		plot_main_c <- paste(trait, "_", chrom, " PGC unfiltered", sep="")

	}else{
		man_name <- paste(trait,".manhattan.unfiltered.jpg", sep="")
		qq_name <- paste(trait,".qq.unfiltered.jpg", sep="")
		plot_main <- paste(trait,"unfiltered", sep=" ")

		#create a plot for corrected pval
		qq_name_c <- paste(trait,".qq.pgc.jpg", sep="")
		man_name_c <- paste(trait,".manhattan.pgc.jpg", sep="")
		plot_main_c <- paste(trait,"PGC unfiltered", sep=" ")
	}

	#jpeg(man_name, width = 1300, height = 600, pointsize = 16)
	#manhattan(pvals, pch=16, main=plot_main)
	#dev.off()

	#create a qqplot for non corrected pval
	#jpeg(qq_name, width = 1300, height = 600, pointsize = 16)
	#qq(pvals$P, main=plot_main)
	#dev.off()

	#create a qqplot for corrected pval
	#jpeg(qq_name_c, width = 1300, height = 600, pointsize = 16)
	#qq(pvals_new$PGC, main=plot_main_c)
	#dev.off()

	#Create a manhattan plot for corrected pvalues
	#jpeg(man_name_c, width = 1300, height = 600, pointsize = 16)
	#manhattan(pvals_new, pch=16, main=plot_main_c)
	#dev.off()

	return (list(pvals=pvals_new,lambda=lambda))
}
