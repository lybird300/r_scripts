#function to read a phenotype table and add a couple of columns
pheno_read <- function(el1){
	pheno_all <- read.table(el1,sep="\t",header=T)
	
	#if there is a age column, remove under 18 people
	if("age" %in% colnames(pheno_all)){
		pheno_all$age <- as.numeric(as.character(pheno_all$age))
		pheno <- pheno_all[which(pheno_all$age >= 18),]
		pheno$age2 <- ((pheno$age)^2)
	}

	pheno$id <- as.character(pheno$id)
	pheno$sex <- as.numeric(as.character(pheno$sex))
	
	pheno
}



#function to add a column based on values of other 2 columns

pheno_add <- function(el1,el2){
		if(is.na(el1) && is.na(el2)){
			all <- as.numeric(el2)
			}
		else if(!is.na(el1) && is.na(el2)){
			all <- el1
			}
		else if(is.na(el1) && !is.na(el2)){
			all <- el2
			}
			all
}

