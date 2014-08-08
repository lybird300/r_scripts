#script to extract the ooverall concordance and discordace for a certain data set

concordance <- read.table("all_samples_all_chr_concordance_discordance_table_no_mismatch.txt",sep="\t", header=T)
overall <- cbind(sum(concordance$SRR_ORR),sum(concordance$SRR_ORA),sum(concordance$SRR_OAA),sum(concordance$SRR_ONA),sum(concordance$SRA_ORR),sum(concordance$SRA_ORA),sum(concordance$SRA_OAA),sum(concordance$SRA_ONA),sum(concordance$SAA_ORR),sum(concordance$SAA_ORA),sum(concordance$SAA_OAA),sum(concordance$SAA_ONA),sum(concordance$SNA_ORR),sum(concordance$SNA_ORA),sum(concordance$SNA_OAA),sum(concordance$SNA_ONA))

colnames(overall) <- colnames(concordance)[5:20]

overall <- as.data.frame(overall)

#calculate overall OGC 
overall_conc_num <- sum(overall['SRR_ORR'],overall['SRA_ORA'],overall['SAA_OAA'])
overall_conc_den <- sum(overall['SRR_ORR'],overall['SRR_ORA'],overall['SRR_OAA'],overall['SRA_ORR'],overall['SRA_ORA'],overall['SRA_OAA'],overall['SAA_ORR'],overall['SAA_ORA'],overall['SAA_OAA'])
OGC <- overall_conc_num/overall_conc_den


#calculate overall NRD
non_ref_num <- sum(overall['SRA_ORR'],overall['SRA_OAA'],overall['SAA_ORR'],overall['SAA_ORA'],overall['SRR_ORA'],overall['SRR_OAA'])
non_ref_den <- sum(overall['SRR_ORA'],overall['SRR_OAA'],overall['SRA_ORR'],overall['SRA_ORA'],overall['SRA_OAA'],overall['SAA_ORR'],overall['SAA_ORA'],overall['SAA_OAA'])
NRD <- non_ref_num/non_ref_den

overall_OGC_NRD <- as.data.frame(cbind(OGC,NRD))
colnames(overall_OGC_NRD) <- c("Concordance","NRdiscordance")

write.table(overall_OGC_NRD,file="overall_concordance_discordance_no_mismatch.txt",sep="\t",col.names=T,quote=F,row.names=F)
