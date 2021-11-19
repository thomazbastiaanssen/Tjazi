impute_zeroes = function(vec, method = "logunif"){
      if(! method %in% c("logunif", "unif", "const")){stop("`method` must be exactly `logunif`, `unif` or `const`")}
      
      #Find detection limit
      DL = min(vec[vec != 0])
      if(method == "logunif"){
        vec[vec == 0] = DL/(10^(runif(n = sum(vec == 0), min =  0, max = 1)))
      }
      else if(method == "unif"){
        vec[vec == 0] = runif(n = sum(vec == 0), min =  0, max = DL)
      }
      else if(method == "const"){
        vec[vec == 0] = 0.65 * DL
      }
      return(vec)
    } 
