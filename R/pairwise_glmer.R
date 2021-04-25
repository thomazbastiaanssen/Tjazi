pairwise_glmer <- function(clr,  
                           model = "~ . + (1|ID)", 
                           metadata, 
                           posthoc.method = "BH", 
                           features.as.rownames = FALSE, CI = TRUE){
  #Generate output data.frame
  out_df = rbind()
  
  #Run the models
  for(feature in 1:nrow(clr)){
    temp_df = metadata
    temp_df$depend = unlist(clr[feature,])
    fit = summary(lmerTest::lmer(formula = paste("depend", model, sep = " "),
                      data    = temp_df))
    
    #Remove the "df" column to standardize results
    fit_coef = coefficients(fit)[,!colnames(coefficients(fit)) == "df"]
    
    #Compute 95% Confidence interval and add to results
    if(CI){fit_coef = cbind(fit_coef, car::Confint(lmer(formula = paste("depend",  model, sep = " "), data = temp_df))[,2:3])}
    
    #Collect the stats
    out_df = rbind(out_df,  t(fit_coef)[1:length(fit_coef)])
    
  }
  
  #Reformat and name the stats part
  out_df = data.frame(out_df)
  colnames(out_df) = paste(rep(row.names(fit_coef), 
                               each = ncol(fit_coef)), 
                           rep(colnames(fit_coef)))
  
  
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
