########################################
# 16/02/2016
# Post analysis QC

rm(list=ls())

args=commandArgs(trailing=TRUE)
require(GWAtoolbox)
cohort <- args[[1]]
chr <- args[[2]]
outpath <- args[[3]]
outfile <- args[[4]]

# build script
to_GWAtoolbox <- paste(outpath,"/results/",chr,"/",cohort,"_GWA_toolbox_chr",chr,".conf",sep="")
write("# Description of input data columns", file=to_GWAtoolbox)

#we need to read from a result file with the following header:
#  SNP,Chromosome,Position,A0,A1,NoMeasured,CallRate,Pexact,MarkerType,Rsq,p,beta,sebeta,effallelefreq,MAF,strand
write("MARKER SNP", append=T, file=to_GWAtoolbox)
write("CHR Chromosome", append=T, file=to_GWAtoolbox)
write("POSITION Position", append=T, file=to_GWAtoolbox)
write("N NoMeasured", append=T, file=to_GWAtoolbox)
write("ALLELE A0 A1", append=T, file=to_GWAtoolbox)
write("STRAND strand", append=T, file=to_GWAtoolbox)
write("EFFECT beta", append=T, file=to_GWAtoolbox)
write("STDERR sebeta", append=T, file=to_GWAtoolbox)
write("PVALUE p", append=T, file=to_GWAtoolbox)
write("FREQLABEL effallelefreq", append=T, file=to_GWAtoolbox)
write("IMPUTED MarkerType", append=T, file=to_GWAtoolbox)
write("IMP_QUALITY Rsq", append=T, file=to_GWAtoolbox)

write("# High quality filters", append=T, file=to_GWAtoolbox)
write("HQ_SNP 0.01 0.3", append=T, file=to_GWAtoolbox)

write("# Plotting filters", append=T, file=to_GWAtoolbox)
write("MAF 0.01 0.05", append=T, file=to_GWAtoolbox)
write("IMP 0.4 0.5", append=T, file=to_GWAtoolbox)

write("# Prefix for output files", append=T, file=to_GWAtoolbox)
write("PREFIX Result_", append=T, file=to_GWAtoolbox)

write("# Input file with GWA data", append=T, file=to_GWAtoolbox)
write("VERBOSITY 2", append=T, file=to_GWAtoolbox)
write("SEPARATOR COMMA", append=T, file=to_GWAtoolbox)
write(paste("PROCESS ",outfile,sep=""), append=T, file=to_GWAtoolbox)

gwasqc(to_GWAtoolbox)
##################################################################
