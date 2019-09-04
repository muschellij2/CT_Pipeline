####################################
## LIME
####################################
#####################################
rm(list=ls())
library(methods)
library(dplyr)
library(fslr)
library(ichseg)
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
mod_type = "rf"

stub_fname = file.path(resdir, 
                       "iso_aggregate_models")
fname = paste0(stub_fname, 
               "_rf_lime.rda")
load(fname)

fdf = fdf %>% 
  mutate(
    pimg = file.path(outdir,
                     paste0(mod_type, "_prob.nii.gz")),
    sm_pimg = sub("_prob", "_sm_prob", 
                  pimg),
    pred_img = sub("prob[.]", "pred.",
                   pimg),
    sm_pred_img = sub("prob[.]", "sm_pred.",
                      pimg),    
  )


##############################
# Keeping files where predictors exist
##############################
imod = as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(imod)) {
  imod = 1
}
L = nrow(fdf)

print(imod)
idf = fdf[imod,]

img = idf$img
mask_fname = sub("[.]nii",
                 "_mask.nii", img)

df = readRDS(idf$outfile)
df = df$df

pimg = readnii(idf$pred_img)
df$pred = c(pimg)

bad = df %>% 
  filter(Y != pred)

bad_run = bad[, names(explainer$feature_distribution)]
bad_run = bad_run[ sample(nrow(bad_run), 10), ]
bad_lime = lime::explain(bad_run, explainer = explainer, n_labels = 1,
                         n_features = 5)



