dice <- function(prediction.obj){
    if (class(prediction.obj) != "prediction") {
        stop(paste(
            "Wrong argument types: First argument", 
            "must be of type", 
            "'prediction'"))    
    }
    argnames <- c()
    x.values <- list()
    y.values <- list()
    for (i in 1:length(prediction.obj@predictions)) {
      fp = prediction.obj@fp[[i]] 
        tp = prediction.obj@tp[[i]]
        fn = prediction.obj@fn[[i]]
        tn = prediction.obj@tn[[i]]
        cutoffs = prediction.obj@cutoffs[[i]]
        meas_dice = 2 * tp / (2*tp + fp + fn)
        x.values <- c(x.values, list(cutoffs))
        y.values <- c(y.values, list(meas_dice))
    }
    if (!(length(x.values) == 0 || 
            length(x.values) == length(y.values))) {
        stop("Consistency error.")
    }
    return(new("performance", x.name = "cutoff", 
        y.name = "dice", 
        alpha.name = "none", 
        x.values = x.values, y.values = y.values, 
        alpha.values = list())
    )
}


sub_samp = function(x, pct = .1, maxcut = 5e3){
    rr = range(x)
    breaks = seq(rr[1], rr[2], length.out = 10)
    cuts = cut(x, breaks= breaks, include.lowest=TRUE)
    levs = levels(cuts)
    l = lapply(levs, function(ilev){
        which(cuts == ilev)
    })
    n.inlev = sapply(l, length)
    n.inlev = ceiling(n.inlev*pct)
    n.inlev = pmin(n.inlev, maxcut)
    ind = sort(unlist(mapply(function(ind, n){
        sample(ind, size = n)
    }, l, n.inlev)))
    return(ind)
}



short_predict = function(object, newdata, 
    lthresh=  .Machine$double.eps^0.5,
    family = NULL){
    tt <- terms(object)
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, 
        na.action = na.pass, 
        xlev = object$xlevels)
    if (is.null(cl <- attr(Terms, "dataClasses"))) 
        stop("no dataclasses")
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    # p <- object$rank
    beta <- object$coefficients
    beta = beta[ !is.na(beta) ]
    predictor = drop(X[, names(beta), drop=FALSE ] %*% beta)

    if (is.null(family)) {
        family = family(object)
    }
    predictor <- family$linkinv(predictor)
    predictor[ predictor < lthresh] = 0
    predictor
}

opt.cut = function(perf, pred){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
}

get_max_dice = function(obj){
    dice.coef = dice(obj)
    ind = which.max(dice.coef@y.values[[1]])
    dice.cutoff = dice.coef@x.values[[1]][[ind]]
    dice.coef = dice.coef@y.values[[1]][ind]
    dice.coef = c(dice=dice.coef, cutoff= dice.cutoff)
    dice.coef = t(dice.coef)
    dice.coef
}

get_pauc = function(obj, fpr.stop){
    pauc = performance(obj, "auc", fpr.stop= fpr.stop)
    pauc = pauc@y.values[[1]] / fpr.stop
    pauc
}

get_acc = function(obj){
    acc = performance(obj, "acc")
    ind = which.max(acc@y.values[[1]])
    cutoff = acc@x.values[[1]][[ind]]
    acc = acc@y.values[[1]][ind]
    acc = t(as.matrix(c(accuracy=acc, cutoff= cutoff)))
    acc
}

get_senscut = function(obj, fpr.stop){
    perf <- performance(obj,"tpr","fpr")
    xind = perf@x.values[[1]] <= fpr.stop
    perf@x.values[[1]] = perf@x.values[[1]][xind]
    perf@y.values[[1]] = perf@y.values[[1]][xind]

    max.cut = obj@cutoffs[[1]][length(perf@y.values[[1]])]
    if (is.infinite(max.cut)){
        max.cut = 1
    }
    preds = obj@predictions[[1]]
    N = length(preds)
    Y = obj@labels[[1]]
    Y = Y == levels(Y)[2]
    res = extrantsr:::sim(
        dman = Y,
        dauto = preds > max.cut)
    sens.cut = matrix(c(res$sens, res$spec, res$accur, max.cut), 
        nrow=1)
    colnames(sens.cut) = c("sensitivity", "specificity", 
        "accuracy", "cutoff")
    return(sens.cut)
}

remove_lmparts = function(mod){
    keep = c("coefficients", "xlevels", 
        "contrasts", "family", "terms")
    mn = names(mod)
    mn = mn[ !( mn %in% keep)]
    for (iname in mn){
        mod[[iname]] = NULL
    }
    mod
}

remove_gamparts = function(mod){
    keep = c("assign", "cmX", "coefficients", "contrasts", 
        "family", 
        "formula", "model", "na.action", "nsdf", 
        "pred.formula", 
        "pterms", "smooth",  "Xcentre", "xlevels", "terms")
    # "Vp",
    mn = names(mod)
    mn = mn[ !( mn %in% keep)]
    for (iname in mn){
        mod[[iname]] = NULL
    }
    mod$model = mod$model[1,]
    mod
}   


get.acc = function(tab){
    sum(diag(tab)) / sum(tab)
}
get.sens = function(tab){
    tab["TRUE", "TRUE"] / sum(tab[, "TRUE"])
}
get.spec = function(tab){
    tab["FALSE", "FALSE"] / sum(tab[, "FALSE"])
}
get.dice = function(tab){
    tt = tab["TRUE", "TRUE"]
    ab <- tt
    estvol = sum(tab[, "TRUE"])
    truevol = sum(tab["TRUE", ])
    aplusb <- (estvol + truevol)
    aorb = sum(tab) - tab["FALSE", "FALSE"]
    dice <- 2 * ab/aplusb
}


limit_pauc = function(pred, fpr.stop){
    perf <- performance(pred,"tpr","fpr")
    for (iperf in seq_along(perf@x.values)){
        ind = which(perf@x.values[[iperf]] <= fpr.stop)
        perf@y.values[[iperf]] = perf@y.values[[iperf]][ind]
        perf@x.values[[iperf]] = perf@x.values[[iperf]][ind]
    }
    return(perf)
}


get_meas = function(xpred) {
    ord = order(xpred, decreasing = TRUE)
    predictions.sorted = xpred[ord]
    tp = cumsum(Y[ord])
    fp = CS - tp
    dups <- rev(duplicated(rev(predictions.sorted)))
    tp <- c(0, tp[!dups])
    fp <- c(0, fp[!dups])
    cutoffs <- c(Inf, predictions.sorted[!dups])
    fn = n.pos - tp
    tn = n.neg - fp

            #########
            # calculate accuracy
            #########
    acc = (tn + tp)/N 
    ind = which.max(acc)
    acc.cutoff = cutoffs[ind]
    acc = c(acc=acc[ind], 
        cutoff= acc.cutoff)
    acc = t(acc)
    colnames(acc) = c("acc", "cutoff")

            #########
            # calculate dice
            #########            
    dice = 2 * tp/(2 * tp + fp + fn)
    ind = which.max(dice)
    dice.cutoff = cutoffs[ind]
    dice = c(dice[ind], dice.cutoff)
    dice = t(dice)
    colnames(dice) = c("dice", "cutoff")


            #########
            # calculate senscut
            #########
    x <- fp/n.neg
    y <- tp/n.pos
            # finite.bool <- is.finite(x) & is.finite(y)
            # x <- x[finite.bool]
            # y <- y[finite.bool]

    ind <- max(which(x <= fpr.stop))
    sens.cut = matrix(c(y[ind], 
        1-x[ind], 
        (tn[ind]+tp[ind])/N, 
        cutoffs[ind]), 
    nrow=1)
    colnames(sens.cut) = c("sensitivity", "specificity", 
        "accuracy", "cutoff")

            #########
            # calculate pAUC
            #########
    tpr.stop <- approxfun(x[ind:(ind + 1)], 
        y[ind:(ind + 
            1)])(fpr.stop)
            x <- c(x[1:ind], fpr.stop)
            y <- c(y[1:ind], tpr.stop) 

            lag = 1
            n = length(x)
            dx = x[(1+lag):n] - x[1:(n-lag)]
            dy = y[(1+lag):n] + y[1:(n-lag)]
            adds = 0.5 * dx * dy
            pauc = sum(adds) / fpr.stop

            l = list(
                dice = dice,
                acc = acc,
                pauc = pauc,
                senscut = sens.cut 
                )              


            return(l)

}



strip_model = function(model){
    model$y = c()
  model$model = c()
  model$residuals = c()
  model$fitted.values = c()
  model$effects = c()
  model$qr$qr = c()  
  model$linear.predictors = c()
  model$weights = c()
  model$prior.weights = c()
  model$data = c()
  attr(model$terms,".Environment") = c()
  attr(model$formula,".Environment") = c()
  model
}

remake_img = function(vec, img, mask = NULL){
    if (is.null(mask)){
        mask = array(1, dim = dim(img))
    }
    img2 = niftiarr(img, 0)
    img2[mask == 1] = vec
    img2 = datatyper(img2)
    img2 = cal_img(img2)
    return(img2)
}


mod_func = function(test.pred, Y, fpr.stop){
    pred <- prediction( 
        test.pred, Y)
    perf <- performance(
        pred,"tpr","fpr")
    pauc.cut = t(opt.cut(
        perf, pred))

    sens.cut = get_senscut(pred, 
        fpr.stop=fpr.stop)

    print(sens.cut)
    pauc = get_pauc(pred, fpr.stop=fpr.stop)
    print(pauc)
    dice.coef = get_max_dice(pred)
    print(dice.coef)
    acc = get_acc(pred)
    print(acc) 

    L = list(
        mod.pauc = pauc, 
        mod.pauc.cut = pauc.cut, 
        mod.sens.cut = sens.cut, 
        mod.dice.coef = dice.coef, 
        mod.acc = acc,
        mod.perf = perf)
    return(L)
}

predict_logistic = function(mod, newdata){
    predict(mod, 
        newdata = newdata, 
        type = "response")
}

predict_gam = function(mod, newdata){
    predict(mod, 
        newdata = newdata, 
        newdata.guaranteed = TRUE,
        type="response", 
        block.size=1e5)
}

predict_rf = function(mod, newdata){
    p = predict(mod,
        newdata = newdata,
        type = "prob"
        )[, "1"]
}
predict_rf_reduced = predict_rf

predict_lasso = function(mod, newdata){
    coefs = rownames(coef(mod))
    coefs = coefs[ 
        !coefs %in% c("(Intercept)")
        ]
    mm = as.matrix(newdata[, coefs])
    predict(mod, 
    newx = mm,
    s = "lambda.1se", 
    type = "response")[, "1"]
}

trim_rf = function(object, printable = TRUE){
    keep_names = c("forest",  "type",  
        "na.action",  "forest$cutoff",  
        "classes",  "terms",  "importance",  
        "coefs")

    if (printable){
        keep_names = c(keep_names, 
        "err.rate", "confusion", 
        "ntree", "mtry", "call")
    }
    keep_names = unique(keep_names)

    nn = names(object)
    sd = setdiff(nn, keep_names)
    for (icn in seq_along(sd)){
        object[[sd[icn]]] = NULL
    }
    object
}
