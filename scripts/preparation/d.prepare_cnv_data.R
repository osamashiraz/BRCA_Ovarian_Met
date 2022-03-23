process_CNV_data = function(){
  # LOAD COPY NUMBER DATA
  load(file = "../input/Nanostring_copynumber.Rdata")
  
  # SIMPLIFY CNV DATA
  mat_simple = as.matrix(CNV_mat)
  mat_simple[CNV_mat <= 1.5] = "LOH"
  mat_simple[CNV_mat <= 1] = "DEL"
  mat_simple[CNV_mat >= 3] = "GAIN"
  mat_simple[CNV_mat >= 5] = "AMP"
  mat_simple[mat_simple %notin% c("AMP","GAIN","LOH","DEL")] = "WT"
  
  return(mat_simple)
}

load(file = "../input/Nanostring_copynumber_coord_map.Rdata")
