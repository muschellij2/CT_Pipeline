myperf <- function (predictions, labels, cutoffs, fp, tp, fn, tn, n.pos, 
    n.neg, n.pos.pred, n.neg.pred, fpr.stop) 
{
    x <- fp/n.neg
    y <- tp/n.pos
    finite.bool <- is.finite(x) & is.finite(y)
    x <- x[finite.bool]
    y <- y[finite.bool]
    if (length(x) < 2) {
        stop(paste("Not enough distinct predictions to compute area", 
            "under the ROC curve."))
    }
    auc <- 0
    for (i in 2:length(x)) {
        auc <- auc + 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 
            1])
    }
    full.auc <- auc
    if (fpr.stop < 1) {
        ind <- max(which(x <= fpr.stop))
        tpr.stop <- approxfun(x[ind:(ind + 1)], y[ind:(ind + 
            1)])(fpr.stop)
        x <- c(x[1:ind], fpr.stop)
        y <- c(y[1:ind], tpr.stop)
    }
    auc <- 0
    for (i in 2:length(x)) {
        auc <- auc + 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 
            1])
    }
    xx <- x[2:length(x)] - x[1:(length(x)-1)]
    yy <- y[2:length(x)] - y[1:(length(x)-1)]
    newauc <- sum(xx*yy*0.5)
    ans <- c(full.auc, auc, newauc)
    names(ans) <- c("full.auc", paste0("pAUC", fpr.stop), "newauc")
    return(ans)
}
