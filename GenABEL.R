args=commandArgs(trailing=TRUE)

#gli args sono:
#args[[1]]: file fenotipo
#args[[2]]: tratto
#args[[3]]: covariate
#args[[4]]: Genomic kinship già calcolata precedentemente
#args[[5]]: Genotypes file
#args[[6]]: cohort

library(GenABEL)
# gkins <- as.matrix(read.table("/nfs/servizio/FVG.kinship"))
# gkins <- as.matrix(read.table(args[[4]]))

# colnames(gkins) <- rownames(gkins)
# chr <- load.gwaa.data(phenofile= "/home/cocca/analyses/MetabolicSyndrome/FVG/fvg_all_metabolic_ALL_MetS_score.csv", genofile = "/home/cocca/analyses/MetabolicSyndrome/FVG/FVG_out")
chr <- load.gwaa.data(phenofile= args[[1]], genofile = args[[5]])

# load("/nfs/servizio/FVG.kinship")
load(args[[4]])
gkins <- kin.matr[chr@phdata$id,chr@phdata$id]

if(args[[3]] == "none"){
   print("No covariates included in polygenic model")
	#qui provo ad utilizzare il nuovo modello poligenico
   model = try(polygenic_hglm(eval(parse(text=args[[2]])) ~ 1, gkins, chr, family=gaussian(link = "identity"), conv=1e-08, maxit=100))
   
   if(is(model,"try-error")){
   print("using old polygenic model")
	#se il nuovo modello da errore uso il vecchio
   model = polygenic(eval(parse(text=args[[2]])), gkins, chr, trait.type = "gaussian", opt.method = "optim", fglschecks=FALSE)
   }
}

if(args[[3]] != "none"){
   print("Covariates are included in polygenic model")
   model = try(polygenic_hglm(eval(parse(text=args[[2]])) ~ eval(parse(text=args[[3]])), gkins, chr, family=gaussian(link = "identity"), conv=1e-08, maxit=100))
   # model = try(polygenic_hglm(eval(parse(text="MetS_score")) ~ eval(parse(text="age")), gkins, chr, family=gaussian(link = "identity"), conv=1e-08, maxit=100))

   if(is(model,"try-error")){
   print("using old polygenic model")
   model = polygenic(eval(parse(text=args[[2]])) ~ eval(parse(text=args[[3]])), gkins, chr, trait.type = "gaussian", opt.method = "optim", fglschecks=FALSE)
   # model = polygenic(eval(parse(text="MetS_score")) ~ eval(parse(text="age")), gkins, chr, trait.type = "gaussian", opt.method = "optim", fglschecks=FALSE)
  }
}

#Get residuals from analysis, based on covariate effects only. 
trait_residuals <- model$residualY 

#Create table with two columns: id and trait of ALL individuals
#using residuals
pheno_residuals <- data.frame(id=chr@phdata$id, trait_residuals=trait_residuals)

#save it into the file. We will use this file in ProbABEL 
write.table(pheno_residuals, file="res.pheno", row.names=F, quote=F)

#get inverse of the variance-covariance matrix 
InvSigma <- model$InvSigma 

rownames(InvSigma) <- na.omit(pheno_residuals)$id

#save it to file. We'll use it in ProbABEL for mmscore 
write.table(InvSigma, file="varcovar.mat", col.names=F, quote=F) 

# Write unfiltered output
GWA <- formetascore(model, data = chr, stat = mmscore, transform = "no", verbosity = 2)
write.table(GWA, file = paste(args[[6]], "_", args[[2]], ".mmscore", sep=""), row.names = F, sep = " ", quote = F)
# write.table(GWA, file = paste("FVG", "_", "MetS_score", ".mmscore", sep=""), row.names = F, sep = " ", quote = F)

system(paste("grep -v name ", args[[6]], "_", args[[2]], ".mmscore | cut -f 1 -d \" \" > genotyped_SNPs.list", sep=""))
# system(paste("grep -v name ", "FVG", "_","MetS_score", ".mmscore | cut -f 1 -d \" \" > genotyped_SNPs.list", sep=""))


# Write filtered output
dim(GWA)

GWA <- GWA[which(GWA$effallelefreq >= 0.01),]
dim(GWA)

GWA <- GWA[which(GWA$pexhwe >= 0.00001),]
dim(GWA)

Quality <- rep(1,nrow(GWA))
Rsq <- rep(1,nrow(GWA))

GWA <- cbind(GWA[,c(1,8,5,9,9)], Quality, Rsq, GWA[,c(10,9,2,3,11,12)])
colnames(GWA) <- c("name","A1","A2","Freq1","MAF","Quality","Rsq","n","Mean_predictor_allele","chrom","position","beta_SNP_add","sebeta_SNP_add")

check_allele = as.character(GWA$A1) == as.character(GWA$A2)
if(length(which(check_allele == "FALSE") == nrow(GWA))){print("Alleles are ok")
  if(length(which(as.numeric(GWA$MAF) <= 0.5) == nrow(GWA))){print("Frequencies are ok")

      # write.table(GWA, file = paste(args[[6]], "_", args[[2]], ".mmscore_MOD", sep=""), row.names = F, sep = " ", quote = F)
      write.table(GWA, file = paste("FVG", "_", "MetS_score", ".mmscore_MOD", sep=""), row.names = F, sep = " ", quote = F)
   }
}   

q(save="no")

