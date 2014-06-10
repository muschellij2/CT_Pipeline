rm(list=ls())
library(cttools)
library(ROCR)
library(fslr)
library(AnalyzeFMRI)

options(fsl.path='/usr/local/fsl')
setwd("~/Desktop/scratch")
# stub = "100-318_20070723_0957_CT_3_CT_Head-"
stub = "100-362_20100126_1926_CT_2_CT_ROUTINE"
img = readNIfTI(paste0(stub, ".nii.gz"), 
                reorient=FALSE)
mask = readNIfTI(paste0(stub, "_SS_Mask_0.01.nii.gz"), 
                 reorient=FALSE)
roi = readNIfTI(paste0(stub, "ROI.nii.gz"), 
                reorient=FALSE)
# erode the mask
mask = fslerode(file=mask, kopts = "-kernel box 1x1x1", 
                  reorient=FALSE, retimg = TRUE)

nvoxels = 1

preds = make_predictors(img, mask, nvoxels = 1)

df = data.frame(preds[ mask > 0, ])
df$Y = roi[mask > 0]

samps = seq(nrow(df)) %in% sample(nrow(df), size= 1e4)
train = df[samps,]
test = df[!samps,]

mod = glm(Y ~ . - zero_100 - gauss, data=train, family=binomial())

test.pred = predict(mod, newdata=test, type="response")
train.pred = predict(mod, newdata=train, type="response")

pred <- prediction( test.pred, test$Y)
perf <- performance(pred,"tpr","fpr", fdr.stop= .1)

auc = performance(pred, "auc")@y.values[[1]]
# plot(perf)
fpr.stop = .1
pauc = performance(pred, "auc", fpr.stop= fpr.stop)
pauc = pauc@y.values[[1]] / fpr.stop
pauc




predimg = img
predimg[!is.na(predimg)] = 0
predimg[is.na(predimg)] = 0
predimg[mask > 0][samps] = train.pred
predimg[mask > 0][!samps] = test.pred
predimg = window_img(predimg, window = c(0, 1))

diff = abs((roi - predimg))
dimg = predimg
dimg@.Data = diff

#########################################################
# Smoothing predictions
#########################################################
prob.img = mean_image(predimg, 1)
pimg = predimg
pimg[!is.na(pimg)] = 0
pimg[is.na(pimg)] = 0
pimg@.Data = prob.img
pimg[abs(pimg) < 1e-13] = 0

smooth.pred = pimg[mask > 0][!samps]
spred <- prediction( smooth.pred, test$Y)
pauc = performance(spred, "auc", fpr.stop= fpr.stop)
pauc = pauc@y.values[[1]] / fpr.stop
pauc

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.259e+01  8.884e-01 -14.173  < 2e-16 ***
#   moment1      1.914e-01  4.119e-02   4.647 3.37e-06 ***
#   moment2      1.603e+00  4.204e-01   3.813 0.000138 ***
#   moment3     -1.328e-01  2.779e-02  -4.778 1.77e-06 ***
#   moment4     -1.274e+00  3.790e-01  -3.362 0.000773 ***
#   pct.thresh   3.502e+00  1.110e+00   3.155 0.001607 ** 
#   pct.0       -3.638e+02  1.901e+04  -0.019 0.984729    
# value       -3.580e-02  2.596e-02  -1.379 0.167798    
# threshTRUE   6.171e-01  5.734e-01   1.076 0.281773 
# Second model
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -15.39703    0.62893 -24.481  < 2e-16 ***
#   moment1       0.19541    0.01942  10.061  < 2e-16 ***
#   moment2       1.10830    0.09911  11.182  < 2e-16 ***
#   moment3      -0.06771    0.01256  -5.391 7.01e-08 ***
#   moment4      -0.86272    0.08059 -10.705  < 2e-16 ***
#   pct.thresh    0.95659    0.75201   1.272  0.20336    
# pct.0       -16.41426    2.19523  -7.477 7.59e-14 ***
#   value         0.02185    0.01268   1.722  0.08501 .  
# threshTRUE    1.33921    0.36124   3.707  0.00021 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
