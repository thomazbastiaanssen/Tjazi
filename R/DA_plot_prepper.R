DA_plot_prepper = function(DA_df, 
                           limit_p     = 0.1, 
                           limit_e     = 0.25, 
                           name_cutoff = ".*ales_", 
                           sig_symbols = c("*"    = 0.1, 
                                           "**"   = 0.01, 
                                           "***"  = 0.001)
                           ){
  
  #Distinguish p-values and effect sizes
  cols_containing_pvals <-  1:((ncol(DA_df))/2)*2
  cols_containing_evals <- (1:((ncol(DA_df))/2)*2) +1
  
  #Change missing values to show no effects
  change_NA_to_1      <- function(x) { replace(x, is.na(x), 1) }
  change_NA_to_0      <- function(x) { replace(x, is.na(x), 0) }
  
  DA_df[,cols_containing_pvals] <- change_NA_to_1(DA_df[,cols_containing_pvals])
  DA_df[,cols_containing_evals] <- change_NA_to_0(DA_df[,cols_containing_evals])
  
  
  #Select features based on p-values and effect sizes
  signig_filter = rep(FALSE, nrow(DA_df))
  for(number in (1:nrow(DA_df))){
    for(pval in (cols_containing_pvals)){
      if(DA_df[number, pval] <= limit_p & abs(DA_df[number, pval + 1]) >= limit_e){
        signig_filter[number] <- TRUE
      }
    }
  }
  
  DA_df <- DA_df[signig_filter,]
  #Shorten names of microbes while keeping taxonomic order
  if(nchar(name_cutoff) > 0){
    DA_df$microbe <- sub(".*ales_", "", DA_df$microbe)
  }
  
  #Reshape data to long format for easy plotting
  plot.m <- tidyr::pivot_longer(DA_df[,c(1, cols_containing_evals)], cols = !microbe)
  pval   <- tidyr::pivot_longer(DA_df[,c(1, cols_containing_pvals)], cols = !microbe)
  
  #Deterimine significance characters
  stars <- rep("", length(pval$value))
  for(number in 1:nrow(pval)){
    if(pval$value[number] <= sig_symbols[1]){
      stars[number] <- "*"
    }
    if(pval$value[number] <= sig_symbols[2]){
      stars[number] <- "**"
    }
    if(pval$value[number] <= sig_symbols[3]){
      stars[number] <- "***"
    }
  }
  
  #Collect output into tibble for easy transfer to ggplot2
  plot.m$stars <- stars
  
  return(plot.m)
}
