rm(list=ls())
library(stargazer)
library(plyr)
# score = "NIHSS"
setwd('~/Dropbox/CTR/DHanley/CT_Registration/CT_Pipeline')

get.stuff = function(score){
  load(paste0("Regress_ROI_", score, "_Results.Rda"))
  best.mod = which.min(aics[, "With_Perc"])
  
  mods = mods[[best.mod]]
  co = mods$coefficients 
  co = rename(co, c("SexMale"="Sex: Male vs. Female", 
               "Base_ICH_10"="TICHVol per 10 cc",
               "perc_ROI" = "ROI Coverage per 10%"))
  mods$coefficients = co
  
  co = keep.cmod$coefficients 
  co = rename(co, c("SexMale"="Sex: Male vs. Female", 
                    "Base_ICH_10"="TICHVol per 10 cc", 
                    "LOCGlobus Pallidus"= 
                      paste0("Reader-Based Location",
                             "&&&& RRRRRRRRRR",
                             "\nRRRRR;RRRRR;Globus Pallidus")))
  names(co) = gsub("^LOC", "RRRRR;RRRRR;", names(co))
  keep.cmod$coefficients = co
  
  res = list(mods, keep.cmod)
}

nihss.res = get.stuff(score = "NIHSS")
gcs.res = get.stuff(score = "GCS")
# colnames(gcs.res[[1]]$model)[1] = "BLAH"
# colnames(gcs.res[[2]]$model)[1] = "BLAH"

gcs.res[[1]]$call = call("lm", 
                    list(formula = Y ~ Age + Sex + Base_ICH_10 + LOC))
gcs.res[[2]]$call = call("lm", 
                    list(formula = Y ~ Age + Sex + Base_ICH_10 + LOC))

cov.name = "ROI Coverage per 10%"
make.coef = function(mod, cov.name){
  nihss.coef = abs(coef(mod)[cov.name])
  nihss.coef = sprintf("%02.1f", nihss.coef)
  nihss.ci  = abs(confint(mod)[cov.name, ])
  nihss.ci = paste(sprintf("%02.1f", nihss.ci), collapse = ", ")
  nihss.ci = paste0("(95\\% CI: ", nihss.ci, ")")
  return(list(coef = nihss.coef, ci = nihss.ci))
}

nihss= make.coef(nihss.res[[1]], cov.name = cov.name)
nihss.coef = nihss$coef
nihss.ci = nihss$ci

nihss= make.coef(nihss.res[[2]], cov.name = "RRRRR;RRRRR;Putamen")
put.nihss.coef = nihss$coef
put.nihss.ci = nihss$ci

cap = paste0("Regression Models for ROI-Based Analysis. The models for ",
             "ROI coverage represent the best model based on the ",
             "model-fit measures. ",
             "We see that after adjusting for age, sex, and ",
            "total baseline ICH volume, increasing 10\\%",
            " coverage is expected to increase NIHSS score by ",
            nihss.coef, " ", nihss.ci, " points. ",
            " We see that all locations, compared to lobar hemorrhages",
            " have higher estimated NIHSS scores, but putaminal ",
            "hemorrhages were significantly higher by ",
            put.nihss.coef, " ", put.nihss.ci, " points.")

gcs= make.coef(gcs.res[[1]], cov.name = cov.name)
gcs.coef = gcs$coef
gcs.ci = gcs$ci

gcs= make.coef(gcs.res[[2]], cov.name = "RRRRR;RRRRR;Putamen")
put.gcs.coef = gcs$coef
put.gcs.ci = gcs$ci

cap2 = paste0("Adjusting for other covariates, increasing 10\\%",
             " coverage is expected to decrease GCS score by ",
             gcs.coef, " ", gcs.ci, " points. ",
             " We see that all locations, compared to lobar hemorrhages",
             " have lower estimated GCS scores, but none were ",
            "statistically different.")
cap = paste(cap, cap2 )

rr = stargazer(nihss.res, gcs.res, type = "latex", 
               title = cap,
               t.auto=FALSE, p.auto=FALSE,
               ci= TRUE, omit.stat="all", 
               single.row = TRUE, star.char="", notes="",
               notes.append=FALSE, 
               notes.label = "", 
               no.space = TRUE, 
               dep.var.caption = "", 
               label = "f:beta",
               omit.table.layout = "n",
               dep.var.labels = c("\\textbf{NIHSS Score}", 
                                  "\\textbf{GCS Score}"),
               column.labels = rep(c("\\textbf{ROI Coverage}", 
                                     "\\textbf{Reader-Location}"), 2),
               model.names = FALSE, model.numbers = FALSE, 
               digits = 1)
rr = gsub("RRRRR", "\\", rr, fixed=TRUE)
move = grep("\\caption", rr)
rr = rr[c(1:(move-1), (move+2):(length(rr)-1), move:(move+1), length(rr))]
empty.hline = grep("\\hline \\\\[-1.8ex]", rr, fixed=TRUE)
l = length(empty.hline)
empty.hline = empty.hline[c(1, (l-1):l)]
rr = rr[-empty.hline]
writeLines(rr, con="Beta_Table.tex")
# }