#' run a model on every row of your data. Made with microbiome data in mind.
#' @export
#' @examples
#'
#' metadata = data.frame(a = sample(letters[1:3], ncol(mtcars), replace=T),
#'                       b = sample(letters[4:6], ncol(mtcars), replace=T),
#'                       c = rnorm(ncol(mtcars)))
#'
#'  fw_fit(x = mtcars, f = x ~ a + b, metadata = metadata, format = "fits", model = "lm", get_CI = T, order = "fact")
#'
fw_fit <- function(x, f, metadata, verbose = T, get_CI = F, format = "wide", model = "lm", order = "fact", adjust.method = "BH", ...){

  stopifnot("Model type not recognised. shoudl be `lm` or `lmer`." = model %in% c("lm", "lmer"))

  stopifnot("Format unrecognised. should be one of `fits`, 'list', `wide` or `brief`. " = format %in% c("fits", "list", "wide", "brief"))

  order = interpret_order(order, model = model, verbose = verbose)

  f = as.formula(f)

  f = update.formula(x ~ 1, f)

  if(verbose){print(paste("Using the following formula:", deparse(f)))
  }

  if(model == "lm"){
    fits = fit_glm(x = x, f = f, metadata = metadata, ...)
  }

  if(model == "lmer"){
    fits    = fit_glmer(x = x, f = f, metadata = metadata, ...)
    if(order$full){
    f.null  = fit_null(x = x, f = f, metadata = metadata, ...)
    if(verbose){
      message("Be extra critical at this point, overall model fitting for mixed effects models is a contested subject. ")
      print(paste("Using the following null formula for conditional tests:", deparse(f.null[[1]]$m.null)))
      print(paste("Using the following null formula for marginal tests:", deparse(f.null[[1]]$m.ran)))
    }
    }
  }

  if(format == "fits"){
    return(fits)
  }

  out_list = make_output(order)

  if(order$full){
    if(model == "lm"){
    out_list$full   = lapply(X = fits, FUN = gather_full)
    }
    if(model == "lmer"){
    out_list$full   = compare_to_null(f.null = f.null, f.full = fits)
    }
  }

  if(order$anovas){
    out_list$anovas = lapply(X = fits, FUN = gather_anova)
  }

  if(order$coefs){
    out_list$coefs  = lapply(X = fits, FUN = gather_coefs, get_CI = get_CI)
  }
  if(order$tukeys){
    out_list$tukeys = lapply(X = fits, FUN = gather_tukeys)
  }

  if(format %in% c("wide", "brief")){

    out_list <- lapply(out_list, coef_to_wide)

    out_df = do.call(data.frame, list(out_list, check.names = F))

    #Organize colnames alphabetically
    out_df = out_df[,order(colnames(out_df))]

    #Add feature names as columns
    out_df = cbind(feature = row.names(out_df), out_df)
    row.names(out_df) <- NULL

    #Adjust for FDR
    out_df = cbind(out_df, adjust_fdr(out_df, method = adjust.method, verbose = verbose))

    if(format == "brief"){
      out_df = out_df[,!grepl("Intercept", colnames(out_df))]

      out_df = out_df[, grepl("feature|p.value|Pr\\(>|Estimate|F value|squared", colnames(out_df))]
    }

    if(sum(unlist(order)) == 1){
      colnames(out_df) <- gsub("full\\.|anovas\\.|coefs\\.|tukeys\\.", "", colnames(out_df))
    }
    return(out_df)
  }

  return(out_list)
}

#' lm wrapper around fw_fit.
#'
#' @export
#'
fw_glm <- function(x, f, metadata, verbose = T, get_CI = T, format = "wide", order = "c", adjust.method = adjust.method, ...){
  fw_fit(x = x, f = f, metadata = metadata, verbose = verbose, get_CI = get_CI, format = format, model = "lm", order = order, ...)
}

#' lmer wrapper around fw_fit.
#'
#' @export
#'
fw_glmer <- function(x, f, metadata, verbose = T, get_CI = T, format = "wide", order = "c", adjust.method = adjust.method, ...){
  fw_fit(x = x, f = f, metadata = metadata, verbose = verbose, get_CI = get_CI, format = format, model = "lmer", order = order, ...)
}

#' Fit a linear model
#'
fit_glm <- function(x, f, metadata, ...)
{apply(X = x, MARGIN = 1, FUN = function(y){
  metadata = cbind(x = c(y), metadata);
  return(lm(formula = f, data = metadata, ...))})}

#' Fit a linear mixed effect model
#'
fit_glmer <- function(x, f, metadata, ...)
{apply(X = x, MARGIN = 1, FUN = function(y){

  metadata = cbind(x = c(y), metadata);
  f.full = suppressMessages(lmerTest::lmer(formula = f, data = metadata, REML = F, ...))
  return(f.full)
  })}

#' Fit a null model for a linear mixed effect model
#'
fit_null <- function(x, f, metadata, ...)

{apply(X = x, MARGIN = 1, FUN = function(y){

  metadata = cbind(x = c(y), metadata);

  m.null <- update.formula(f, ~ 1 )
  f.null = lm(formula = m.null, data = metadata, ...)

  t = as.character(terms(f)[[3]])

  reff = t[grep("\\|", x = t)]
  reff = gsub(")", "", gsub(".*\\|", "", reff))

  m.ran <- update.formula(f, formula(paste("~ 1", paste(paste0("factor(", reff, ")"), collapse = " + "), sep = " + ")))
  fr.null = lm(formula = m.ran, data = metadata, ...)
  return(list(f.null = f.null, fr.null = fr.null, m.null = m.null, m.ran = m.ran))
})}

#' Fit a null model for a linear mixed effect model
#'
compare_to_null <- function(f.null, f.full){

  mapply(function(f.null, f.full) {
    out_vec <- vector("numeric", length = 4)
    out_vec[1] <- cor(residuals(f.null$f.null), fitted(f.full))^2
    out_vec[2] <- anova(f.full, f.null$f.null)["f.full", "Pr(>Chisq)"]
    out_vec[3] <- cor(residuals(f.null$fr.null), fitted(f.full))^2
    out_vec[4] <- anova(f.full, f.null$fr.null)["f.full", "Pr(>Chisq)"]


    names(out_vec) <- c("conditional.r.squared", "conditional.p.value", "marginal.r.squared", "marginal.p.value")

    return(out_vec)
  }, f.null = f.null, f.full = f.full, SIMPLIFY = F)
}

#' Calculate p-value for entire model fit
#'
get_model_p <- function(x){
  fstat = summary(x)$fstatistic
  p.val = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  return(p.val)
}

#' Extract and format full model results
#'
gather_full <- function(x){
  coef_out <- vector(mode = "numeric", length = 3)

  coef_out[1:2] <- unlist(summary(x)[c("r.squared", "adj.r.squared")])
  coef_out[3]   <- get_model_p(x)

  names(coef_out) <- c("r.squared", "adj.r.squared", "p.value")

  return(coef_out)
}

#' Extract and format coefficients.
#'
gather_coefs <- function(x, get_CI = T){
  coef_out <- coef(summary(x))
  if(get_CI){
    CI = car::Confint(x)
    CI = CI[row.names(coef_out),]
    coef_out <- cbind(CI, coef_out[,-1])
  }
  return(coef_out)
}

#' Extract and format anova results
#'
gather_anova <- function(x){
  #return values for non-residual coeficients
  coef_out <- anova(x)[1:sum(!is.na(anova(x)$`Pr(>F)`)),]
  return(coef_out)
}


#' Extract and format TukeyHSD results
#'
gather_tukeys <- function(x){
  coef_list = TukeyHSD(aov(x))
  coef_out  = do.call(rbind, coef_list)
  return(coef_out)
}

#' Take all coefficients and prepare them for wide format. Called by `coef_to_wide`.
#'
unroll_names <- function(x){
  if(length(dim(x)) == 2){
    suppressWarnings(paste(row.names(x), rep(colnames(x), each = nrow(x))))
  }
  else if (is.null(dim(x))){
    names(x)
  }
}

#' Take all coefficients and prepare them for wide format. Called by `coef_to_wide`.
#'
unroll_coefs <- function(x){
  out <- unlist(c(x))
  names(out) = unroll_names(x)
  return(out)
}

#' Take all coefficients and prepare them for wide format.
#'
coef_to_wide <- function(x){

  coef_list <- lapply(X = x, FUN = unroll_coefs)
  do.call(rbind, coef_list)

}

#' Go from vector to individual booleans
#' @param x A string of input flags
#' @return A list of orders
#'
interpret_order <- function(x, model, verbose){
  out_list = vector("list", length = 4)
  names(out_list) = c("full", "anovas", "coefs", "tukeys")

  out_list["full"  ] <- grepl("f", x)
  out_list["anovas"] <- grepl("a", x)
  out_list["coefs" ] <- grepl("c", x)
  out_list["tukeys"] <- grepl("t", x)

  if(model == "lmer" & out_list$tukeys){
    if(verbose){print("lmer is not compatible with TukeysHSD, so skipping that order.")};
    out_list["tukeys"] <- F
  }

  return(out_list)

}

#' Make a container list based on the order
#' @param x the order list output from `interpret_order`.
#'
make_output <- function(x){
  ord_vec <- unlist(x)
  ord_vec = ord_vec[ord_vec]

  out_list = vector("list", length = length(ord_vec))
  names(out_list) = names(ord_vec)

  return(out_list)
}

#' Adjust p.values for wide format
#'
adjust_fdr <- function(x, method = method, verbose = verbose){
  y = x[,grep(pattern = "p.value|Pr\\(>", x = colnames(x))]

  if(method %in% c("Storeys.Q", "Storey", "storey")){
    if(verbose){print("Adjusting for FDR using Storey's q-value procedure.")}

    out_df <- apply(X = y,
                    MARGIN = 2,
                    FUN = function(x){qvalue::qvalue(p = x, lambda = seq(0, max(x), 0.05))$qvalues})

  } else {
    if(verbose & method == "BH")
      {print("Adjusting for FDR using Benjamini & Hochberg's procedure.")}
  out_df <- apply(X = y,
                  MARGIN = 2,
                  FUN = function(x){p.adjust(x, method = method)})
  }

  colnames(out_df) <- paste(colnames(out_df), method, sep = ".")


  return(out_df)
}
