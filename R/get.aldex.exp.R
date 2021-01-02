get.aldex.exp = function(clr, useMC = TRUE, verbose = TRUE){
  is.multicore = FALSE
  if ("BiocParallel" %in% rownames(installed.packages()) & 
      useMC) {
    message("multicore environment is OK -- using the BiocParallel package")
    is.multicore = TRUE
  }
  else {
    if (verbose == TRUE) 
      message("operating in serial mode")
  }
  
  
  nr <- ALDEx2::numFeatures(clr)
  rn <- ALDEx2::getFeatureNames(clr)
  cn <- ALDEx2::getSampleIDs(clr)
  
  if (is.multicore == TRUE) 
    clr.list <- BiocParallel::bplapply(ALDEx2::getMonteCarloInstances(clr), function(m) {
      t(apply(m, 1, median))
    })
  if (is.multicore == FALSE) 
    clr.list <- lapply(ALDEx2::getMonteCarloInstances(clr), function(m) {
      t(apply(m, 1, median))
    })
 
  
  clr.exp = data.frame(matrix(unlist(clr.list), nrow = nr, byrow=T),stringsAsFactors=FALSE)
  row.names(clr.exp) = rn
  colnames(clr.exp)  = cn
  
  return(clr.exp)
}
