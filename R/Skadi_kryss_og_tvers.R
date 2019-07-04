skadi_kryss <- function(x_vector, y_metric, method = "spearman"){
  res_df_cor = data.frame(p.value   = rep(NA,   nrow(x_vector)), 
                          statistic = rep(NA,   nrow(x_vector)), 
                          out.index = rep(TRUE, nrow(x_vector)))
  row.names(res_df_cor) <- rownames(x_vector)
  
  for(microbe in 1:(nrow(x_vector))){
    skresp = skadi(x = unlist(x_vector[microbe,]), 
                   y = y_metric, 
                   method = method, 
                   euclid.outlier.check = T, diagnostic.plot = F,
                   max.distance = 1, give.uncorrected.p.value = F, 
    )
    skres = skadi(x = unlist(x_vector[microbe,]), 
                  y = y_metric, 
                  method = method, 
                  euclid.outlier.check = T,
                  max.distance = 1, give.uncorrected.p.value = T, 
                  xlab=paste(strsplit(row.names(x_vector[microbe,]), split = ".*D_4__" ),
                             "\n p = ", skresp$p.value))
    
    res_df_cor[microbe,]$p.value   = skres$p.value
    res_df_cor[microbe,]$statistic = skres$estimate
    res_df_cor[microbe,]$out.index = skres$outlier
  }
  res_df_cor$q.value = qvalue(p = res_df_cor$p.value)$qvalues
  
  return(res_df_cor)
}

skadi_tvers <- function(x_vector, y_vector, method = "spearman"){
  q_df <- data.frame(names = row.names(x_vector))
  r_df <- data.frame(names = row.names(x_vector))
  
  for(column in 1:ncol(y_vector)){
    skadi_kryss_output <-skadi_kryss(x_vector = x_vector, 
                                     y_metric = unlist(y_vector[,column]), 
                                     method = method)
    
    q_df[,column] <- skadi_kryss_output$q.value
    r_df[,column] <- skadi_kryss_output$statistic
  }
  colnames(q_df) <- colnames(y_vector)
  row.names(q_df)<- row.names(x_vector)
  colnames(r_df) <- colnames(y_vector)
  row.names(r_df)<- row.names(x_vector)
  
  resultList <- list("q_values" = q_df, 
                     "correlation_statistics" = r_df)
  
  return(resultList)
}



testList    <- skadi_tvers(x_vector = species.exp, 
                           y_vector = metadata)
