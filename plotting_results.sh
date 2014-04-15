#!/usr/local/bin/bash
#ARGS:
#$1=trait
#$2=filepath

trait=$1
res_path=$2

#IFS=$'\n'
#set $(cat parameter_file.txt)

# RUN GenABEL.R
#bsub -J "G_$2_$4" -o GenABEL.log  -M7000000 -R"select[mem>7000] rusage[mem=7000]" \
#-q long R CMD BATCH '--args '$6' '$4' '$8' '${10}' '${12}' '$2'' /nfs/team151/GWAS/IMPUTATION-1KG/VALBORBERA/GWAS/GenABEL.R

# RUN ProbABEL
#for i in $(seq 1 22)
#do
#bsub -J "P_$2_$4_chr${i}" -w "ended(G_$2_$4)" -o chr${i}.log -M7000000 -R"select[mem>7000] \
#rusage[mem=7000]" -q basement /nfs/users/nfs_g/gp5/ProbABEL/bin/palinear --mmscore varcovar.mat --pheno res.pheno \
#--chrom ${i} --map ${14}/${i}/mach.${i}.machlegend \
#--info ${16}/${i}/mach.${i}.machinfo \
#--dose ${18}/output${i}.gen.dose.fvi --out chr${i}.palinear
#done

# RUN merge.R
#bsub -o merge.log -J "merge_$2_$4" -w "ended(P_$2_$4_chr*)" \
#-M10000000 -R"select[mem>10000] rusage[mem=10000]" -q hugemem \
#R CMD BATCH '--args '$2' '$4'' /nfs/team151/GWAS/IMPUTATION-1KG/VALBORBERA/GWAS/merge.R

# RUN exclude.pl
#bsub -o exclude.log -J "PERL_$2_$4" -w "ended(merge_$2_$4)" \
#-M6000000 -R"select[mem>6000] rusage[mem=6000]" -q long \
#perl /nfs/users/nfs_a/ar10/exclude.pl all.palinear.out all.palinear.filtered

# RUN FORMATTING AND PLOTTING STEPS
for chr in {1..22} X
do
bsub -J "plot_${trait}_${chr}" -o "%J_plot_${trait}_${chr}.log" -e "%J_plot_${trait}_${chr}.err" \
-M6000000 -R"select[mem>6000] rusage[mem=6000]" -q basement \
R CMD BATCH '--args '${trait}' '${chr}' '${res_path}/${trait}.chr${chr}.tab.assoc.txt'' /nfs/users/nfs_m/mc14/Work/r_scripts/formatting_plot.R

done

