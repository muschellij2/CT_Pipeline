####################################
## LIME
####################################
#####################################
rm(list=ls())
library(methods)
library(fslr)
library(ichseg)
library(ROCR)
library(matrixStats)
library(randomForest)
library(lime)


rootdir = file.path("/Volumes/DATA_LOCAL", 
                    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  rootdir = Sys.getenv("dex")
}
basedir = file.path(rootdir, 
                    "PITCH_reconverted", 
                    "processed_data")
resdir = file.path(rootdir, 
                   "PITCH_reconverted", 
                   "results")
progdir = file.path(rootdir,
                    "PITCH_reconverted",
                    "code")
source(file.path(progdir, 
                 "iso_performance_functions.R"))

filename = file.path(resdir, 
                     "iso_filename.rds")
fdf = readRDS(filename)


##############################
# Keeping files where predictors exist
##############################
fname = file.path(resdir, 
                  "iso_candidate_aggregate_data.rda")
x = load(fname)

train$gr_medztemp = NULL
train$mode = fdf$mode[match(train.img, 
                            fdf$img)]
test$mode = fdf$mode[match(test.img, 
                           fdf$img)]

runnames = colnames(train)
nosmooth = c("any_zero_neighbor",
             "thresh", "pct_zero_neighbor")
runnames = runnames[ !(runnames %in% 
                         c("mask", "Y", "img", nosmooth))]


sds = colSds(as.matrix(train))
names(sds) = colnames(train)
novar = sds == 0
novar = colnames(train[novar])
novar = novar[!novar %in% c("mask")]
novar = c("", novar)
novar = paste0(novar, collapse=" - ")
formstr = paste0("Y ~ . - mask", novar) 

test$multiplier = test.mult.df[, 
                               "multiplier"]

stub_fname = file.path(resdir, 
                       "iso_aggregate_models")
fname = paste0(stub_fname, 
               "_rf.Rda")
what = load(fname)

model = modlist$mod
predict_model.randomForest <- function(x, newdata, type, ...) {
  res <- predict(x, newdata = newdata, type = ifelse(type == "raw", "response", type))
  switch(type,
         raw = data.frame(Response = res, stringsAsFactors = FALSE),
         prob = as.data.frame(res, check.names = FALSE)
  )
}

model_type.randomForest <- function(x, ...) "classification"

train = train[ , !colnames(train) %in% "mask"]
explainer = lime(train, model = model)

fname = paste0(stub_fname, 
               "_rf_lime.rda")
save(explainer, predict_model.randomForest,
     model_type.randomForest,
     file = fname)

