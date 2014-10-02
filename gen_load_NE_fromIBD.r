######################################################################################################
# Plot 10: Ne from IBD
###############################################################################
# We calculated the IBD sharing segments using GERMLINE , than the sharing density for each chromosome.
# Now we'd like to: visualize this sharing density, select regions to remove because of too much sharing
# 
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")

# pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
pops <- c("ALL")
base_folder <- getwd()
minDens <- 75

ylab <- "Sharing density"
xlab <- "Chrom segment"

for (chr in seq(1,22)){
  # current_shared_path <- paste(base_folder,"/SHARE_D",sep="")
  current_shared_path <- paste(base_folder,"/CHR",chr,sep="")
  main <- paste("IBD sharing density across chr",chr)
  current_chr_all_pop <- NULL
  current_chr_all_pop_names <- NULL
  current_chr_all_pop_sd <- NULL
  print(chr)

  for (pop in pops){
    current_chr_current_pop_path <- paste(current_shared_path,"/",pop,".",chr,".non_missing.match.",minDens,".shareDens",sep="")
    # current_chr_current_pop_path <- paste(current_shared_path,"/",pop,".",chr,".non_missing.filtered.match.shareDens",sep="")
    current_chr_current_pop <- read.table(current_chr_current_pop_path,header=F)
    current_chr_current_pop$cohort <- pop
    colnames(current_chr_current_pop) <- c("bin","share","cohort")
    # print(pop)
    current_chr_all_pop <- rbind(current_chr_all_pop,current_chr_current_pop)
    current_sd <- sd(current_chr_current_pop$share)
    current_mean <- mean(current_chr_current_pop$share)
    current_chr_current_pop_sd <- data.frame(share_sd=current_sd,mean=current_mean,mean_sd=(current_mean+current_sd),
      sd2=(current_mean+(2*current_sd)),
      # sd_2=(current_mean-(2*current_sd)),
      sd5=(current_mean+(5*current_sd)),
      # sd_5=(current_mean-(5*current_sd)),
      cohort=pop)
    current_chr_all_pop_sd <- rbind(current_chr_all_pop_sd, current_chr_current_pop_sd)
    # we need to remove regions with more than 5sd sharing density
    # first extract those:
    current_chr_regions_to_keep <- current_chr_current_pop[which(current_chr_current_pop$share < (current_mean+(5*current_sd))),]
    len <- 0
    reg <- 0
    reg_to_keep <- NULL
    current_reg_to_keep <- NULL

    for (i in 1:(dim(current_chr_regions_to_keep)[1])){
      
      if(i!=(dim(current_chr_regions_to_keep)[1])){
        
        if ((current_chr_regions_to_keep$bin[i+1] - current_chr_regions_to_keep$bin[i]) == 0.5){
          len <- len+0.5
          current_reg_to_keep <- rbind(current_reg_to_keep,current_chr_regions_to_keep[(i),])
        }else{
          # print(len)
          # print(dim(current_reg_to_keep))
          
          if (len < 45){
            len <- 0
            current_reg_to_keep <- NULL
          }else{
            print(len)
            len <- 0
            reg <- reg+1
            write.table(current_reg_to_keep,file=paste(current_chr_current_pop_path,"_R",reg,".to_include",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
            reg_to_keep <- rbind(reg_to_keep,current_reg_to_keep)
            current_reg_to_keep <- NULL
          }
        }  
      }else{
        # print(chr)
        # print(dim(current_chr_regions_to_keep)[1])
        # print(dim(current_reg_to_keep))
        
        if (len < 45){
          len <- 0
          current_reg_to_keep <- NULL
        }else{
          print(len)
          current_reg_to_keep <- rbind(current_reg_to_keep,current_chr_regions_to_keep[(i),])
          reg_to_keep <- rbind(reg_to_keep,current_reg_to_keep)
          reg <- reg+1
          write.table(current_reg_to_keep,file=paste(current_chr_current_pop_path,"_R",reg,".to_include",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
           
        }
      }
    }
    reg_to_keep$cohort <- NULL
    write.table(reg_to_keep,file=paste(current_chr_current_pop_path,".to_include",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
  }

  pl <- ggplot(current_chr_all_pop) + aes(x = bin, y = share,col=cohort) + geom_point(size = 3) + xlab(xlab) + ylab(ylab) + ggtitle(main) + facet_wrap(~cohort)
  
  # # jpeg(paste(base_folder,"/",chr,"_point_dens.jpg",sep=""),width=1800, height=800)
  ggsave(filename=paste(base_folder,"/",chr,"_point_dens.",minDens,".jpeg",sep=""),width=10, height=7,dpi=300,plot=pl + geom_hline(aes(yintercept = mean_sd), data=current_chr_all_pop_sd,linetype= 2) + geom_hline(aes(yintercept = sd2), data=current_chr_all_pop_sd,linetype= 3) + geom_hline(aes(yintercept = sd5), data=current_chr_all_pop_sd,linetype= 4))
}


##################################################################################
# Same as above but reading regions's names and files from a list
rm(list=ls())
source("/nfs/users/nfs_m/mc14/Work/r_scripts/col_pop.r")

# pops <- c("CEU","TSI","VBI","FVG","CARL","Erto","Illegio","Resia","Sauris")
pops <- c("ALL")
base_folder <- getwd()
minDens <- 0

ylab <- "Sharing density"
xlab <- "Chrom segment"

share_list <- read.table(paste(base_folder,"share_region.list",sep="/"),header=F)
share_list$V1 <- as.character(share_list$V1)
share_list$V2 <- as.character(share_list$V2)

for (j in 1:length(share_list$V1)) {
  current_shared_path <- paste(base_folder,sep="")
  main <- paste("IBD sharing density across region",share_list[j,2])
  current_chr_all_pop <- NULL
  current_chr_all_pop_names <- NULL
  current_chr_all_pop_sd <- NULL

  for (pop in pops){
    current_chr_current_pop_path <- paste(current_shared_path,"/",share_list[j,1],sep="")
    # current_chr_current_pop_path <- paste(current_shared_path,"/",pop,".",chr,".non_missing.filtered.match.shareDens",sep="")
    current_chr_current_pop <- read.table(current_chr_current_pop_path,header=F)
    current_chr_current_pop$cohort <- pop
    colnames(current_chr_current_pop) <- c("bin","share","cohort")
    # print(pop)
    current_chr_all_pop <- rbind(current_chr_all_pop,current_chr_current_pop)
    current_sd <- sd(current_chr_current_pop$share)
    current_mean <- mean(current_chr_current_pop$share)
    current_chr_current_pop_sd <- data.frame(share_sd=current_sd,mean=current_mean,mean_sd=(current_mean+current_sd),
      sd2=(current_mean+(2*current_sd)),
      # sd_2=(current_mean-(2*current_sd)),
      sd5=(current_mean+(5*current_sd)),
      # sd_5=(current_mean-(5*current_sd)),
      cohort=pop)
    current_chr_all_pop_sd <- rbind(current_chr_all_pop_sd, current_chr_current_pop_sd)
    # we need to remove regions with more than 5sd sharing density
    # first extract those:
    current_chr_regions_to_keep <- current_chr_current_pop[which(current_chr_current_pop$share < (current_mean+(5*current_sd))),]
    len <- 0
    reg <- 0
    reg_to_keep <- NULL
    current_reg_to_keep <- NULL

    for (i in 1:(dim(current_chr_regions_to_keep)[1])){
      
      if(i!=(dim(current_chr_regions_to_keep)[1])){
        
        if ((current_chr_regions_to_keep$bin[i+1] - current_chr_regions_to_keep$bin[i]) == 0.5){
          len <- len+0.5
          current_reg_to_keep <- rbind(current_reg_to_keep,current_chr_regions_to_keep[(i),])
        }else{
          # print(len)
          # print(dim(current_reg_to_keep))
          
          if (len < 45){
            len <- 0
            current_reg_to_keep <- NULL
          }else{
            len <- 0
            reg <- reg+1
            write.table(current_reg_to_keep,file=paste(current_chr_current_pop_path,"_R",reg,".to_include",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
            reg_to_keep <- rbind(reg_to_keep,current_reg_to_keep)
            current_reg_to_keep <- NULL
          }
        }  
      }else{
        print(len)
        print(dim(current_reg_to_keep))
        
        if (len < 45){
          len <- 0
          current_reg_to_keep <- NULL
        }else{
          current_reg_to_keep <- rbind(current_reg_to_keep,current_chr_regions_to_keep[(i),])
          reg_to_keep <- rbind(reg_to_keep,current_reg_to_keep)
          reg <- reg+1
          write.table(current_reg_to_keep,file=paste(current_chr_current_pop_path,"_R",reg,".to_include",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
           
        }
      }
    }
    reg_to_keep$cohort <- NULL
    write.table(reg_to_keep,file=paste(current_chr_current_pop_path,".to_include",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
  }

  pl <- ggplot(current_chr_all_pop) + aes(x = bin, y = share,col=cohort) + geom_point(size = 3) + xlab(xlab) + ylab(ylab) + ggtitle(main) + facet_wrap(~cohort)
  
  # # jpeg(paste(base_folder,"/",chr,"_point_dens.jpg",sep=""),width=1800, height=800)
  ggsave(filename=paste(base_folder,"/REG_",share_list[j,2],"_point_dens.",minDens,".jpeg",sep=""),width=10, height=7,dpi=300,plot=pl + geom_hline(aes(yintercept = mean_sd), data=current_chr_all_pop_sd,linetype= 2) + geom_hline(aes(yintercept = sd2), data=current_chr_all_pop_sd,linetype= 3) + geom_hline(aes(yintercept = sd5), data=current_chr_all_pop_sd,linetype= 4))
}

#resolve IBD sharing plateau problem
whole_reg_file <- paste(getwd(),"ALL.8.non_missing_rs62485438-rs117143950.ped.filtered.match.snp_dens",sep="/")
plateau_reg_file <- paste(getwd(),"ALL.8.non_missing_rs62485438-rs117143950.ped.filtered.match.plateau.snp_dens",sep="/")

whole_reg_file <- paste(getwd(),"ALL.9.non_missing_rs75838906-rs9792573.ped.filtered.match.snp_dens",sep="/")
plateau_reg_file <- paste(getwd(),"ALL.9.non_missing_rs75838906-rs9792573.ped.filtered.match.plateau.snp_dens",sep="/")


whole_reg <- read.table(whole_reg_file,header=F)
plateau_reg <- read.table(plateau_reg_file,header=F)
dim(whole_reg)
dim(plateau_reg)
summary(whole_reg$V4)
summary(plateau_reg$V4)

> summary(whole_reg$V4)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
9.123e-05 6.333e-04 1.378e-03 4.804e-03 3.467e-03 5.840e-02 
> summary(plateau_reg$V4)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0005622 0.0119700 0.0149300 0.0200600 0.0197800 0.0584000 
> 
