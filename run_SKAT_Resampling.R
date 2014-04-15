library(SKAT,lib.loc="/software/rarevar/r_libs/")
sessionInfo()

args = (commandArgs(TRUE));

OutType = as.character(args[1]);
Permutes = as.numeric(args[2]);
Method = as.character(args[3]);
File.Bed = args[4];
File.Bim = args[5];
File.Fam = args[6];
File.SetID = args[7];
File.SSD = args[8];
File.Info = args[9];

if (length(args)==10)
{
	File.cov = args[10];
}

#OutType = paste("\"",OutType,"\"",sep="")

cat("Input parameters:\n","out_type              = ",OutType,"\n","n.Resampling              = ",Permutes,"\n","File.Bed              = ",File.Bed,"\n","File.Bim            = ",File.Bim,"\n","File.Fam    = ",File.Fam,"\n","File.SetID   = ",File.SetID,"\n","File.SSD        = ",File.SSD,"\n","File.Info = ",File.Info,"\n",sep="");

output.file = sub(".SSD",".results.txt", File.SSD)

##############################################
# 	Generate SSD file

# To use binary ped files, you have to generate SSD file first.
# If you already have a SSD file, you do not need to call this function. 

if (file.info(File.SetID)$size == 0) {q(save="no")};

Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
FAM<-read.table(File.Fam,header=FALSE)

##############################################
# 	Open SSD, and Run SKAT
#
y<-FAM[,6]

#
# To use a SSD file, please open it first. 
# 
SSD.INFO<-Open_SSD(File.SSD, File.Info)

# Number of samples 
SSD.INFO$nSample 

# Number of Sets
SSD.INFO$nSets

output <-file(output.file,"w"); 

if (exists("File.cov") && file.exists(File.cov))
{	cat("Running with covariates..\n")
	X<-read.table(File.cov,header=FALSE)
	#print(X)
	obj <- SKAT_Null_Model(y ~ as.matrix(X), out_type = OutType, n.Resampling = Permutes, type.Resampling = "bootstrap")
} else {
	obj <- SKAT_Null_Model(y ~ 1, out_type = OutType, n.Resampling = Permutes, type.Resampling = "bootstrap")
}

#p.val1<-rep(0,SSD.INFO$nSets)
##p.val2<-rep(0,SSD.INFO$nSets)
#for(i in 1:SSD.INFO$nSets){

#	SetID<-SSD.INFO$SetInfo$SetID[i]
#	p.val1[i]<-SKAT.SSD.OneSet(SSD.INFO,SetID, obj)$p.value

#	cat(SetID,p.val1[i],sep="\t",file=output);
#	cat("\n",file=output);
	#or 
	#p.val2[i]<-SKAT.SSD.OneSet_SetIndex(SSD.INFO,i, y ,out_type = "C")$p.value

#}
#p.val1
#p.val2

############################################
# Use different types of weight

#p.val1<-rep(0,SSD.INFO$nSets)
#p.val3<-rep(0,SSD.INFO$nSets)

#names(p.val3)<-SSD.INFO$SetInfo$SetID

methods=unlist(strsplit(Method,","))

cat("SetID",paste(c("SKAT_"),rep(methods,each=2),c("","_perm"),sep=""),sep="\t",file=output)
cat("\n",file=output);
	
for(i in 1:SSD.INFO$nSets){
	
	SetID<-SSD.INFO$SetInfo$SetID[i]
	Z<-Get_Genotypes_SSD(SSD.INFO, i)
	
	pvalues <- SetID
	
	for (m in 1:length(methods))
	{
		re <- SKAT(Z, obj, kernel = "linear.weighted", method=methods[m], estimate_MAF=2)
		p.val1<-re$p.value
		Get_Resampling_Pvalue(re)
		p.val3<-Get_Resampling_Pvalue(re)$p.value
		
		pvalues <- rbind(pvalues,p.val1,p.val3)	
	}
	
	cat(pvalues,sep="\t")
	cat("\n")
	
	cat(pvalues,sep="\t",file=output);
	cat("\n",file=output);

}


close(output);	

warnings()

# You must close SSD after you finish to use it.
Close_SSD()

q(save="no");
