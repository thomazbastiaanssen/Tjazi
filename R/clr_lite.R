clr_lite = function(counts, samples_are = "cols", method = "unif"){
  
  if(! samples_are %in% c("rows", "cols")){stop("`samples_are` must be exactly `rows` or `cols`")}
  
  if(samples_are == "rows"){counts = t(counts)}
  #impute zeroes
    for(col in 1:ncol(counts)){
      counts[,col] = impute_zeroes(vec = counts[,col], method = method)
    }
  if(samples_are == "rows"){counts = t(counts)}
  #Apply CLR transformation  
  if(samples_are == "cols"){counts = t(counts)}
      clr_counts = data.frame(compositions::clr(counts))
  if(samples_are == "cols"){
    clr_counts = t(clr_counts)
    row.names(clr_counts) = row.names(counts)
  }
    
    return(data.frame(clr_counts))
    }
