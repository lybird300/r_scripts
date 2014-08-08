#define function for population colors in plots
col_pop <- function(pops){
  all_cols <- NULL

  for(i in 1:length(pops)){
  if(pops[i] == "CEU"){
    cur_col <- "#E8D0A9"
    density <- 1
  }
  if(pops[i] == "TSI"){
    cur_col <- "#B7AFA3"
    density <- 1
  }
  if(pops[i] == "FVG"){
    cur_col <- "#5A5E61"
    density <- 1
  }
  if(pops[i] == "FVG_p"){
    cur_col <- "#5A5E61"
    density <- 2
  }
  if(pops[i] == "FVG_s"){
    cur_col <- "#5A5E61"
    density <- 3
  }
  if(pops[i] == "FVG_n"){
    cur_col <- "#5A5E61"
    density <- 4
  }
  if(pops[i] == "Erto"){
    cur_col <- "#96D4CA"
    density <- 1
  }
  if(pops[i] == "Illegio"){
    cur_col <- "#95C0BB"
    density <- 1
  }
  if(pops[i] == "Resia"){
    cur_col <- "#85C0E7"
    density <- 1
  }
  if(pops[i] == "Sauris"){
    cur_col <- "#66939E"
    density <- 1
  }
  if(pops[i] == "VBI"){
    cur_col <- "#DF5E5E"
    density <- 1
  }
  if(pops[i] == "VBI_p"){
    cur_col <- "#DF5E5E"
    density <- 2
  }
  if(pops[i] == "VBI_s"){
    cur_col <- "#DF5E5E"
    density <- 3
  }
  if(pops[i] == "VBI_n"){
    cur_col <- "#DF5E5E"
    density <- 4
  }
   if(pops[i] == "CARL"){
    cur_col <- "#F7A6A6"
    density <- 1
  }
  if(pops[i] == "CARL_p"){
    density <- 2
    cur_col <- "#F7A6A6"
  }
  if(pops[i] == "CARL_s"){
    density <- 3
    cur_col <- "#F7A6A6"
  }
  if(pops[i] == "CARL_n"){
    cur_col <- "#F7A6A6"
    density <- 4
  }
  pop_col <- cbind(cur_col,pops[i],density)
  all_cols <- rbind(all_cols,pop_col)
}
all_cols <- as.data.frame(all_cols)
colnames(all_cols) <- c("color","pop","density")
all_cols$color <- as.character(all_cols$color)
all_cols$pop <- as.character(all_cols$pop)
all_cols$density <- as.character(all_cols$density)


return (all_cols)
}
