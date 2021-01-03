pairwise_DA_tester = function (clr, groups, comparisons, 
                             verbose = TRUE, parametric = T, 
                             ignore.posthoc = F, 
                             paired.test = FALSE){
  
  #Define colnames for output
  pval_name = "p.value"
  if(!ignore.posthoc){pval_name = "BH.adjusted.p.value"}
  
  #Create output data.frame
  out_df = data.frame(microbe = rownames(clr))
  
  for(comparison in 1:nrow(comparisons)){
    
    #Name the comparison for output data.frame
    comp_name <- paste(comparisons[comparison, 1], 
                       comparisons[comparison, 2], 
                       sep = " vs ")

    #Perform the appropriate t-test
    if(parametric){
    pvals = lapply(1:nrow(clr), function(i) t.test(
      as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,1]]))), 
      as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,2]]))), 
      paired = paired.test)[c("p.value")])
    } else {
      pvals = lapply(1:nrow(clr), function(i) wilcox.test(
        as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,1]]))), 
        as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,2]]))), 
        paired = paired.test)[c("p.value")])
    }
    pvals = unlist(pvals)
    

    #Calculate appropriate effect size
    if(sum(c(groups == comparisons[comparison,1]),
             groups == comparisons[comparison,2]) > 40 ){
    evals = lapply(1:nrow(clr), function(i) effsize::cohen.d(
      d = as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,1]]))),
      f = as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,2]]))),
      paired = paired.test,
      hedges.correction = FALSE)[c("estimate")])
    } else {
      evals = lapply(1:nrow(clr), function(i) effsize::cohen.d(
        d = as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,1]]))),
        f = as.numeric(as.character(unlist(clr[i,groups == comparisons[comparison,2]]))),
        paired = paired.test,
        hedges.correction = TRUE)[c("estimate")])
    }
    evals = unlist(evals)
    
    
    #Perform FDR procedure
    if(!ignore.posthoc){
      pvals = p.adjust(p = pvals, method = "BH")
    }
    
    #Gather results per comparison
    comp_res = data.frame(a = unlist(pvals), 
                          b = unlist(evals)) 
    colnames(comp_res) = c(paste(comp_name, pval_name), 
                           comp_name)
    
    #Add results to output
    out_df = cbind(out_df, comp_res)
  }
  
  return(out_df)
  
}
