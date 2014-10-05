#define function for population colors in plots

col_pop <- function(pops){
  all_cols <- NULL
  CAR_col = "#F54E4E"
  CAR_p_col = "#F77171"
  CAR_s_col = "#F99595"
  CAR_n_col = "#FBB8B8"
  CEU_col = "#13256F"
  Erto_col = "#10DFBC"
  FVG_col ="#505E69"
  FVG_p_col ="#737E87"
  FVG_s_col ="#969EA5"
  FVG_n_col ="#B9BFC3"
  Illegio_col = "#48867F"
  Resia_col = "#16C224"
  Sauris_col = "#00631F"
  TSI_col = "#619FE0"
  VBI_col = "#CA8A1A"
  VBI_p_col = "#D5A148"
  VBI_s_col = "#DFB976"
  VBI_n_col = "#EAD0A3"

  for(i in 1:length(pops)){
  if(pops[i] == "CEU"){
    cur_col <- CEU_col
    density <- 8
  }
  if(pops[i] == "TSI"){
    cur_col <- TSI_col
    density <- 8
  }
  if(pops[i] == "FVG"){
    cur_col <- FVG_col
    density <- 8
  }
  if(pops[i] == "FVG_p"){
    cur_col <- FVG_p_col
    density <- 3
  }
  if(pops[i] == "FVG_s"){
    cur_col <- FVG_s_col
    density <- 2
  }
  if(pops[i] == "FVG_n"){
    cur_col <- FVG_n_col
    density <- 1
  }
  if(pops[i] == "Erto" || pops[i] == "FVE"){
    cur_col <- Erto_col
    density <- 8
  }
  if(pops[i] == "Illegio" || pops[i] == "FVI"){
    cur_col <- Illegio_col
    density <- 8
  }
  if(pops[i] == "Resia" || pops[i] == "FVR"){
    cur_col <- Resia_col
    density <- 8
  }
  if(pops[i] == "Sauris" || pops[i] == "FVS"){
    cur_col <- Sauris_col
    density <- 8
  }
  if(pops[i] == "VBI"){
    cur_col <- VBI_col
    density <- 8
  }
  if(pops[i] == "VBI_p"){
    cur_col <- VBI_p_col
    density <- 3
  }
  if(pops[i] == "VBI_s"){
    cur_col <- VBI_s_col
    density <- 2
  }
  if(pops[i] == "VBI_n"){
    cur_col <- VBI_n_col
    density <- 1
  }
   if(pops[i] == "CARL" || pops[i] == "CAR" ){
    cur_col <- CAR_col
    density <- 8
  }
  if(pops[i] == "CARL_p" || pops[i] == "CAR_p"){
    density <- 3
    cur_col <- CAR_p_col
  }
  if(pops[i] == "CARL_s" || pops[i] == "CAR_s"){
    density <- 2
    cur_col <- CAR_s_col
  }
  if(pops[i] == "CARL_n" || pops[i] == "CAR_n"){
    cur_col <- CAR_n_col
    density <- 1
  }
  pop_col <- cbind(cur_col,pops[i],density)
  all_cols <- rbind(all_cols,pop_col)
}
all_cols <- as.data.frame(all_cols)
colnames(all_cols) <- c("color","pop","density")
all_cols$color <- as.character(all_cols$color)
all_cols$pop <- as.character(all_cols$pop)
all_cols$density <- as.numeric(as.character(all_cols$density))


return (all_cols)
}
