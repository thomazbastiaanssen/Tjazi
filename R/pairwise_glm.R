pairwise_glm <- function(clr, 
                         model = "~ .", 
                         metadata, 
                         posthoc.method = "BH", 
                         features.as.rownames = FALSE){
  #Generate output data.frame
  out_df = rbind()
  
  #Run the models
  for(feature in 1:nrow(clr)){
    temp_df = metadata
    temp_df$depend = unlist(clr[feature,])
    fit = summary(glm(formula = paste("depend", model, sep = " "),
                      data    = temp_df))
    
    #Collect the stats
    out_df = rbind(out_df,  t(coefficients(fit))[1:length(coefficients(fit))])
    
  }
  
  #Reformat and name the stats part
  out_df = data.frame(out_df)
  colnames(out_df) = paste(rep(row.names(coefficients(fit)), 
                                   each = ncol(coefficients(fit))), 
                               rep(colnames(coefficients(fit))))
  
  
  #Perform FDR
  pvals <- colnames(out_df)[grepl("Pr\\(>", colnames(out_df))]
  df.bh <- out_df[, pvals]
  if(posthoc.method == "Storey"){
    colnames(df.bh) <- paste0(colnames(df.bh), ".Storeys.Q")
    for (j in 1:ncol(df.bh)) {
      df.bh[, j] <- qvalue::qvalue(p = df.bh[, j])$qvalues
      }
    } else {
      colnames(df.bh) <- paste(colnames(df.bh), posthoc.method, sep = ".")
      for (j in 1:ncol(df.bh)) {
        df.bh[, j] <- p.adjust(df.bh[, j], method = posthoc.method)
        }
    }
  
  #Consolidate results
  out_df = cbind(out_df, df.bh)
  
  if(features.as.rownames){
    row.names(out_df) = rownames(clr)
    
  } else {
    out_df = cbind(data.frame(microbe = rownames(clr)), out_df)
    
  }
  
  return(out_df)
}
