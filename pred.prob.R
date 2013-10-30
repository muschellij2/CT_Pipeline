get.stuff <- function(mod){
    var.used <- attr(mod$terms, "term.labels")
    var.classes <- attr(mod$terms, "dataClasses")
    var.classes <- var.classes[var.used]
  }
pred.prob <- function(model, test){

  var.classes <- get.stuff(model)


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
  tt <- test[, names(coefs)]

  ### getting predictions - prob slower than matrix multiplication
  ### but more control
  for (icoef in names(coefs)){
    tt[, icoef] <- coefs[names(coefs) == icoef] * tt[, icoef]
    print(icoef)
  }
  pred <- rowSums(tt)
  pred <- pred + intercept
  pred <- exp(pred)/(1+exp(pred))
  return(pred)
}