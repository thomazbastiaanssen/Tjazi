#' @export
pairwise_glm <- function (clr, y = "microbe", model = "~ microbe",
                          metadata, posthoc.method = "BH", family = gaussian(link = "identity"),
                          features.as.rownames = FALSE, CI = TRUE, verbose = TRUE)
{
  if(verbose){print(family)}
  out_df = rbind()
  if(y == "microbe"){
    if(verbose){print(paste0("Running models with formula: microbe ", model))}
    for (feature in 1:nrow(clr)) {
      temp_df = metadata
      temp_df$microbe = unlist(clr[feature, ])
      temp_df$depend  = unlist(clr[feature, ])
      fit = summary(glm(formula = paste("depend", model, sep = " "),
                        data = temp_df, family = family))
      fit_coef = coefficients(fit)[, !colnames(coefficients(fit)) == "df"]
      if(CI){fit_coef = cbind(fit_coef, car::Confint(glm(formula = paste("depend", model, sep = " "),
                                                         data = temp_df, family = family))[,2:3])}
      out_df = rbind(out_df, t(fit_coef)[1:length(fit_coef)])
    }
  }
  else if(y %in% colnames(metadata) ){
    if(verbose){
      print(paste("Running models with formula:",  y , model,  sep = " "))
    }
    for (feature in 1:nrow(clr)) {
      temp_df = metadata
      temp_df$microbe = unlist(clr[feature, ])
      colnames(temp_df)[colnames(temp_df) == y] = "depend"
      fit = summary(glm(formula = paste("depend", model, sep = " "),
                        data = temp_df, family = family))
      fit_coef = coefficients(fit)[, !colnames(coefficients(fit)) == "df"]
      if(CI){fit_coef = cbind(fit_coef, car::Confint(glm(formula = paste("depend", model, sep = " "),
                                                         data = temp_df, family = family))[,2:3])}
      out_df = rbind(out_df, t(fit_coef)[1:length(fit_coef)])

    }
  }
  out_df = data.frame(out_df)

  colnames(out_df) = paste(rep(row.names(fit_coef), each = ncol(fit_coef)),
                           rep(colnames(fit_coef)))
  pvals <- colnames(out_df)[grepl("Pr\\(>", colnames(out_df))]
  df.bh <- out_df[, pvals]
  if (posthoc.method == "Storey") {
    if(verbose){print("Adjusting for FDR using Storey's q-value procedure.")}
    colnames(df.bh) <- paste0(colnames(df.bh), ".Storeys.Q")
    for (j in 1:ncol(df.bh)) {
      pvals = df.bh[, j]
      df.bh[, j] <- qvalue::qvalue(p = pvals, lambda = seq(0, max(pvals), 0.05))$qvalues
    }
  }
  else {
    colnames(df.bh) <- paste(colnames(df.bh), posthoc.method,
                             sep = ".")
    if(verbose){print(paste0("Adjusting for FDR using ", posthoc.method, " procedure."))}
    for (j in 1:ncol(df.bh)) {
      df.bh[, j] <- p.adjust(df.bh[, j], method = posthoc.method)
    }
  }
  out_df = cbind(out_df, df.bh)
  if (features.as.rownames) {
    row.names(out_df) = rownames(clr)
  }
  else {
    out_df = cbind(data.frame(microbe = rownames(clr)), out_df)
  }
  return(out_df)
}
