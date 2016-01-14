rm(list=ls())
library(buRlo)
args=commandArgs(trailing=TRUE)

#gli args sono:
#args[[1]]: file fenotipo
#args[[2]]: tratto
#args[[3]]: covariate
#args[[4]]: Genomic kinship
#args[[5]]: Genotypes file
#args[[6]]: cohort
#args[[7]]: chromosome
#args[[8]]: Imputation path and prefix
# gkins <- as.matrix(read.table("/nfs/servizio/FVG.kinship"))
# gkins <- as.matrix(read.table(args[[4]]))

# colnames(gkins) <- rownames(gkins)
# chr <- load.gwaa.data(phenofile= "/home/cocca/analyses/MetabolicSyndrome/FVG/fvg_all_metabolic_ALL_MetS_score.csv", genofile = "/home/cocca/analyses/MetabolicSyndrome/FVG/FVG_out")
pheno <- args[[1]]
trait <- args[[2]]
covariates <- args[[3]]
kinship <- args[[4]]
geno <- args[[5]]
cohort <- args[[6]]
chromosome <- args[[7]]
imp_path <- args[[8]]

gwa_data <- load.gwaa.data(phenofile= pheno, genofile = geno)

# load("/nfs/servizio/FVG.kinship")
load(kinship)
gkins <- kin.matr[chr@phdata$id,chr@phdata$id]
# formula <- paste(trait,"~",covariates,sep=" ")
formula <- eval(parse(text=trait)) ~ eval(parse(text=covariates))

GWA(formula, gwa_data, gkins, chr=chromosome, trait.type="gaussian", plot=T, an.type="GRAMMAR", imp.extens=imp_path, model.SNP="additive", Rsq=0, maf=0.00005, extens=paste(cohort,trait,chromosome,sep="_"), p.annot=1e-5)

q(save="no")

