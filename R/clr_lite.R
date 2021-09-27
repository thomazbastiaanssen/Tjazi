clr_lite = function(counts, samples_are = "cols", method = "unif"){
  temp_counts = counts
  if(! samples_are %in% c("rows", "cols")){stop("`samples_are` must be exactly `rows` or `cols`")}
  
  if(samples_are == "rows"){temp_counts = t(temp_counts)}
  #impute zeroes
    for(col in 1:ncol(temp_counts)){
      temp_counts[,col] = impute_zeroes(vec = temp_counts[,col], method = method)
    }
  if(samples_are == "rows"){temp_counts = t(temp_counts)}
  #Apply CLR transformation  
  if(samples_are == "cols"){temp_counts = t(temp_counts)}
      clr_counts = data.frame(compositions::clr(temp_counts))
  if(samples_are == "cols"){clr_counts = t(clr_counts)}
    clr_counts = data.frame(clr_counts)
  
    rownames(clr_counts) = rownames(counts)
    colnames(clr_counts) = colnames(counts)

    return(clr_counts)
    }
