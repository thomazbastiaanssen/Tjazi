impute_zeroes = function(vec, method = "unif"){
  if(! method %in% c("unif", "const")){stop("`method` must be exactly `unif` or `const`")}

  #Find detection limit
  DL = min(vec[vec != 0])
  
  if(method == "unif"){
     vec[vec == 0] = runif(n = sum(vec == 0), min =  0.1 * DL, max = DL)
  }
  else if(method == "const"){
    vec[vec == 0] = 0.65 * DL
}
  return(vec)
}
