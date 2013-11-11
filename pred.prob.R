get.stuff <- function(mod){
    var.used <- attr(mod$terms, "term.labels")
    var.classes <- attr(mod$terms, "dataClasses")
    var.classes <- var.classes[var.used]
  }

pred.prob <- function(model, test){

  var.classes <- get.stuff(model)
  coefs <- model$coef
  names(coefs) <- model$vars

  ### checking to see if we used factors or not in the model
  stopifnot(all(var.classes %in% c('numeric', 'logical')))
  coefs <- coef(model)
  have.log <- var.classes %in% "logical"
  if (any(have.log)){
    vars <- names(var.classes)[have.log]
    for (ivar in vars){
      names(coefs) <- gsub(paste0(ivar, "TRUE"), ivar, names(coefs))
    }
  }
  ## just getting the intercept - not needed in matrix mult
  intercept <- ifelse("(Intercept)" %in% names(coefs), 
    coefs["(Intercept)"], 
    0)
  coefs <- coefs[ !names(coefs) %in% "(Intercept)" ]
  tt <- test

  ### getting predictions - prob slower than matrix multiplication
  ### but more control
  if (length(coefs) > 0){
    if ("data.table" %in% class(tt)){
        expr <- paste0(coefs, "*", names(coefs), collapse="+" )
        expr <- parse(text=expr)
        tt <- tt[, eval(expr)]
    } else {
        for (icoef in names(coefs)){
          tt[, icoef] <- coefs[names(coefs) == icoef] * tt[, icoef, drop=FALSE]
          print(icoef)
        }
        tt<- rowSums(tt[, names(coefs), drop=FALSE])
    }
  } else {
    tt <- rep(0, nrow(tt))
  }
  pred <- tt + intercept
  pred <- exp(pred)/(1+exp(pred))
  return(pred)
}

### just get the coefficients (without the factor crap)
scrape.mod <- function(x){

  var.used <- attr(x$terms, "term.labels")
  add.int <- "(Intercept)" %in% x$names
  if (add.int) var.used <- c("(Intercept)", var.used)

  cn <- x$names
  if (!x$converged) {
    coefs <- rep(NA, length(cn))
  } else {
    coefs <- coef(mod)
  }
  names(coefs) <- var.used

  l <- list(coefficients=coefs, null.deviance=NA, deviance=NA, aic=NA, 
    converged=x$converged)
  return(l)
}

scrape.mod.glm <- function(x){
  cn <- names(x)
  cn <- cn[! cn %in% c("coefficients", 
    "null.deviance", 
    "deviance",
    "aic")]
  for (icn in cn) x[icn] <- NULL
  x
}

scrape.mod.sub <- function(x) {
    x$data <- NULL
    x$residuals <- NULL
    x$fitted.values <- NULL
    x$prior.weights <- NULL
    x$model <- NULL
    x$linear.predictors <- NULL
    x$weights <- NULL
    x$effects <- NULL
    x$qr <- NULL
    x$y <- NULL
    x
  }