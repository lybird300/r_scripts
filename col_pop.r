#define function for population colors in plots

col_pop <- function(pops){
  all_cols <- NULL
  CAR_col = "#F54E4E"
  CAR_p_col = "#F77171"
  CAR_s_col = "#F99595"
  CAR_n_col = "#FBB8B8"
  CEU_col = "#13256F"
  FVG_col ="#505E69"
  FVG_p_col ="#737E87"
  FVG_s_col ="#969EA5"
  FVG_n_col ="#B9BFC3"
  Erto_col = "#10DFBC"
  Illegio_col = "#48867F"
  Resia_col = "#16C224"
  Sauris_col = "#00631F"
  Erto_n_col = "#10DFBC"
  Illegio_n_col = "#48867F"
  Resia_n_col = "#16C224"
  Sauris_n_col = "#00631F"
  Erto_p_col = "#10DFBC"
  Illegio_p_col = "#48867F"
  Resia_p_col = "#16C224"
  Sauris_p_col = "#00631F"
  Erto_s_col = "#10DFBC"
  Illegio_s_col = "#48867F"
  Resia_s_col = "#16C224"
  Sauris_s_col = "#00631F"
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
  if(pops[i] == "Erto" || pops[i] == "FVE" || pops[i] == "FVG-E" ){
    cur_col <- Erto_col
    density <- 8
  }
  if(pops[i] == "Illegio" || pops[i] == "FVI" || pops[i] == "FVG-I" ){
    cur_col <- Illegio_col
    density <- 8
  }
  if(pops[i] == "Resia" || pops[i] == "FVR" || pops[i] == "FVG-R" ){
    cur_col <- Resia_col
    density <- 8
  }
  if(pops[i] == "Sauris" || pops[i] == "FVS" || pops[i] == "FVG-S" ){
    cur_col <- Sauris_col
    density <- 8
  }
  if(pops[i] == "Erto_n" || pops[i] == "Erto_s" || pops[i] == "Erto_p"){
    cur_col <- Erto_col
    density <- 1
  }
  if(pops[i] == "Illegio_n" || pops[i] == "Illegio_s" || pops[i] == "Illegio_p"){
    cur_col <- Illegio_col
    density <- 1
  }
  if(pops[i] == "Resia_n" || pops[i] == "Resia_s" || pops[i] == "Resia_p"){
    cur_col <- Resia_col
    density <- 1
  }
  if(pops[i] == "Sauris_n" || pops[i] == "Sauris_s" || pops[i] == "Sauris_p"){
    cur_col <- Sauris_col
    density <- 1
  }
  if(pops[i] == "Erto_s" || pops[i] == "Erto_p"){
    cur_col <- Erto_col
    density <- 2
  }
  if(pops[i] == "Illegio_s" || pops[i] == "Illegio_p"){
    cur_col <- Illegio_col
    density <- 2
  }
  if(pops[i] == "Resia_s" || pops[i] == "Resia_p"){
    cur_col <- Resia_col
    density <- 2
  }
  if(pops[i] == "Sauris_s" || pops[i] == "Sauris_p"){
    cur_col <- Sauris_col
    density <- 2
  }
  if( pops[i] == "Erto_p"){
    cur_col <- Erto_col
    density <- 3
  }
if( pops[i] == "Illegio_p"){
    cur_col <- Illegio_col
    density <- 3
  }
if( pops[i] == "Resia_p"){
    cur_col <- Resia_col
    density <- 3
  }
if( pops[i] == "Sauris_p"){
    cur_col <- Sauris_col
    density <- 3
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
  pop_col <- cbind(cur_col,pops[i],density)
  all_cols <- rbind(all_cols,pop_col)
}
all_cols <- as.data.frame(all_cols)
colnames(all_cols) <- c("color","pop","density")
all_cols$color <- as.character(all_cols$color)
all_cols$pop <- as.character(all_cols$pop)
all_cols$density <- as.numeric(as.character(all_cols$density))

cols <- c(all_cols$color)
names(cols) <- all_cols$pop

return (cols)
}
