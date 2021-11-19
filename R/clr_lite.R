clr_lite = function(counts, samples_are = "cols", method = "logunif", replicates = 1000) 
{
  temp_counts = counts
  
  if(! method %in% c("logunif", "unif", "const"))
  {stop("`method` must be exactly `logunif`, `unif` or `const`")}
  
  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }
  
  temp_counts = apply(X          = temp_counts, 
                      MARGIN     = 2, 
                      FUN        = clr_imputed, 
                      method     = method, 
                      replicates = replicates)
  
  if(samples_are == "rows"){
    temp_counts = data.frame(t(temp_counts))
  }
  
  clr_counts = data.frame(temp_counts)
  rownames(clr_counts) = rownames(counts)
  colnames(clr_counts) = colnames(counts)
  return(clr_counts)
}
