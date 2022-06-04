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
fw_fit <- function(x, f, metadata, verbose = T, get_CI = T, format = "wide", model = "lm", order = "fact", ...){

  stopifnot("Model type not recognised. shoudl be `lm` or `lmer`." = model %in% c("lm", "lmer"))

  stopifnot("Format unrecognised. should be one of `fits`, 'list', `long` or `wide`. " = format %in% c("fits", "list", "long", "wide"))

  order = interpret_order(order)

  f = as.formula(f)

  f = update.formula(x ~ 1, f)

  if(verbose){print(paste("Using the following formula:", deparse(f)))
  }

  if(model == "lm"){
    fits = fit_glm(x = x, f = f, metadata = metadata, ...)
  }
  if(model == "lmer"){
    message("Error: lmer version not implemented yet!")
    break()
  }
  if(format == "fits"){
    return(fits)
  }

  out_list = make_output(order)

  if(order$full){
    out_list$full   = lapply(X = fits, FUN = gather_full)
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

  if(format == "wide"){

    out_list <- lapply(out_list, coef_to_wide)

    out_df = do.call(data.frame, list(out_list, check.names = F))
    out_df = cbind(feature = row.names(out_df), out_df)
    return(out_df)
  }

  return(out_list)
}


#' Fit a linear model
#'
fit_glm <- function(x, f, metadata, ...)
{apply(X = x, MARGIN = 1, FUN = function(y){
  metadata = cbind(x = c(y), metadata);
  return(lm(formula = f, data = metadata, ...))})}

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
  coef_out <- anova(x)[-nrow(anova(x)),]
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
interpret_order <- function(x){
  out_list = vector("list", length = 4)
  names(out_list) = c("full", "anovas", "coefs", "tukeys")

  out_list["full"  ] <- grepl("f", x)
  out_list["anovas"] <- grepl("a", x)
  out_list["coefs" ] <- grepl("c", x)
  out_list["tukeys"] <- grepl("t", x)

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