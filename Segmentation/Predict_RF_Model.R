###########################################################
## This code is for prediction of Image Segmentation of CT
##
## Author: John Muschelli
## Last updated: May 20, 2014
##########################################################
##########################################################
rm(list=ls())
library(plyr)
library(cttools)
library(fslr)
library(ROCR)
library(matrixStats)
library(randomForest)
library(getopt)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

segdir = file.path(progdir, "Segmentation")
source(file.path(segdir, "performance_functions.R"))

outdir = file.path(basedir, "results")

correct = "Rigid"
# options = c("none", "N3", "N4", "N3_SS", "N4_SS",
#         "SyN", "SyN_sinc", "Rigid", 
# "Affine", "Rigid_sinc", 
#         "Affine_sinc")
options = c("none", "N3_SS", "N4_SS",
    "Rigid",  "Rigid_sinc")
# options = c("none", "Rigid_sinc")
# options = "Rigid"

spec = matrix(c(
    'correct', 'c', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

correct = opt$correct
print(opt)

my.tab <- function(
  x, 
  y, 
  dnames=c("x", "y")) {
  x = as.numeric(x)
  y = as.numeric(y)
  stopifnot(all(unique(c(x,y)) %in% c(0, 1, NA)))
  tt = sum(x * y)
  t1=sum(x)
  t2=sum(y)
  tab = matrix(c(length(x)-t1-t2+tt, 
    t1-tt, t2-tt, tt), 2, 2)
  n = list(c("FALSE", "TRUE"), c("FALSE", "TRUE"))
  names(n) = dnames
  dimnames(tab) = n
  tab = as.table(tab)
  return(tab) 
}

types = c("_zval2", "_zval_all", "_zval2_medztemp")
# , "_zval2"ls
# "_include_all", 
type = types[1]


keep.obj = ls()

# for (correct in options){
    
    all.obj = ls()
    rm.obj = all.obj[!(all.obj %in% 
        c(keep.obj, "keep.obj"))]
    rm(list=rm.obj)
    for (i in 1:3) gc()
    correct = match.arg(correct, options)
    adder = switch(correct, 
        "none"= "",
        "N3"="_N3",
        "N4" = "_N4",
        "N3_SS" = "_N3_SS",
        "N4_SS" = "_N4_SS", 
        "SyN" = "_SyN",
        "SyN_sinc" = "_SyN_sinc",
        "Rigid" = "_Rigid",
        "Affine" = "_Affine",
        "Rigid_sinc" = "_Rigid_sinc",
        "Affine_sinc" = "_Affine_sinc")

    #### load voxel data
    outfile = file.path(outdir, "Voxel_Info.Rda")
    load(file=outfile )

    outfile = file.path(outdir, 
        "111_Filenames_with_volumes_stats.Rda")
    load(file = outfile)


    makedir = sapply( fdf$outdir, function(x) {
        if (!file.exists(x)){
            dir.create(x, showWarnings =FALSE)
        }
    })
    irow = 1
    x = fdf[irow,]

    ##############################
    # Keeping files where predictors exist
    ##############################
    outfiles = nii.stub(basename(fdf$img))
    outfiles = paste0(outfiles, "_predictors", 
        adder, ".Rda")
    outfiles = file.path(fdf$outdir, outfiles)
    stopifnot(file.exists(outfiles))

    # load(file = file.path(outdir, 
        # "Segmentation_Models.Rda"))
    ##############################
    # Run lmod number of models - 
    # not all the models - leave out
    ##############################
    mod.filename = file.path(outdir, 
        paste0("Collapsed_Models", adder, ".Rda"))
    xxxx = load(mod.filename)
    vol.sdatas = vol.datas = vol.data
    reses = sreses = res

    cut.vol.data = cut.vol.sdata = cut.vol.tsdata = vol.data
    pauc.cut.vol.data = pauc.cut.vol.sdata = cut.vol.data
    pauc.cut.vol.tsdata = pauc.cut.vol.sdata
    sens.cut.vol.data = sens.cut.vol.sdata = cut.vol.data
    dice.cut.vol.data = dice.cut.vol.sdata = cut.vol.data

    get.pred <- as.numeric(Sys.getenv("SGE_TASK_ID"))
    if (is.na(get.pred)) get.pred = 38
    x = fdf[get.pred,]
    print(get.pred)

# for (get.pred in runpreds){

    iddir = fdf$iddir[get.pred]
    id.outdir = fdf$outdir[get.pred]
    predname = nii.stub(basename(fdf$img[get.pred]))
    predname = file.path(id.outdir, 
        paste0(predname, "_predictors", adder, ".Rda"))
    load(predname)
    df = img.pred$df
    stopifnot(all(df$Y %in% c(0, 1)))

    nim = img.pred$nim
    keep.ind = img.pred$keep.ind
    rm(list="img.pred")
    for (i in 1:3) gc() 
    df$include = df$value >= 30 & df$value <= 100
    
    df$mode = fdf$mode[get.pred]

    fname = file.path(outdir, 
        paste0("Aggregate_data_cutoffs", adder, ".Rda"))

    load(file = fname)
    keepnames = colnames(est.cutoffs)
    keepnames = keepnames[!( keepnames %in% 
    c("dist_centroid", "smooth10", "smooth20", "moment2", 
        "moment4", "moment3"))]    

    med.ztemp = median(df$zscore_template[keep.ind])
    df$gr_medztemp = (df$zscore_template > med.ztemp)

    df$skew[is.nan(df$skew)] = 0
    df$kurtosis[is.nan(df$kurtosis)] = 0

    include = rep(TRUE, length=nrow(df))
    for (icut in keepnames){
        qcuts = est.cutoffs[, icut]
        colname = paste0(icut, ".cutoff")
        df[, colname] = df[, icut] >= qcuts[1] & 
            df[, icut] <= qcuts[2]
        include = include & df[, colname]
    }

    df$include.all = include

    df$zval = df[, "zscore3.cutoff"] & df$include &
        df$pct_thresh.cutoff
    df$zval2 = df[, "zscore2.cutoff"] & df$zval
    # df$dist_centroid <= 75
    df$zval_all = df[, "zscore_template.cutoff"] & df$zval2

    df$zval2_medztemp = df$zval2 &  df$gr_medztemp

    # df$subset = df$zval2
    df$subset = df$zval2

    df$in0100 = df$value >= 0 & df$value <= 100
    # df$in20_85 = df$value >= 20 & df$value <= 85
    df$mask = df$mask > 0
    
    ######################################
    # Get volume of ROI
    ######################################  
    vdim = voxdim(nim)
    vres = prod(vdim) / 1000

    Y =  df$Y

    ######################################
    # Get volume of ROI, using varslice if needed
    ######################################  
    vol.roi = sum(Y) * vres 

    if (correct %in% c("none", "N3_SS", "N4_SS")) {
        x$hdr = file.path(x$iddir, "Sorted",
            paste0(nii.stub(x$img, bn=TRUE), 
                "_Header_Info.Rda"))
        hdrload = load(x$hdr)
        img = nim
        img[is.na(img)] = 0
        img@.Data = array(df$Y, dim=dim(img))
        vol.roi = get_roi_vol(img, dcmtables)$truevol
    }


    not0100 = sum(df$Y[ !df$in0100 ])

    brain.vol = sum(df$mask) * vres
    if (correct %in% c("none", "N3_SS", "N4_SS")) {
        img = nim
        img[is.na(img)] = 0
        img@.Data[keep.ind] = 1
        brain.vol = get_roi_vol(img, dcmtables)$truevol
    }


    ######################################
    # Keep all ROI = 1, even if not inmask
    ######################################  
    roi.not.in = which(df$Y == 1)
    roi.not.in = roi.not.in[!(roi.not.in %in% keep.ind)]
    keep.ind = sort(c(keep.ind, roi.not.in))


    #### need this because the length of df has changed
    roi.not.in = which(keep.ind %in% roi.not.in)

    df$subset[ roi.not.in ] = FALSE


    Y =  df$Y
    benchmark = 1-mean(df$Y)

    df$img = fdf$img[1]

    fpr.stop = .01

    cc = complete.cases(df)
    stopifnot(all(cc))

    test.gam.pred = rep(0, length=nrow(df))
    test.rf.pred = test.gam.pred

    rf.pred = predict(rf.mod,
        newdata = df[ df$subset, ],
        type = "prob"
        )[, "1"]
    test.rf.pred[ df$subset ] = rf.pred

    preds = matrix(test.rf.pred, ncol = 1)
    colnames(preds) = "rf"

    cut.filename = file.path(outdir, 
    paste0("Model_Cutoffs", adder, ".Rda"))

    x = load(file=cut.filename)

    scut.filename = file.path(outdir, 
    paste0("Smooth_Model_Cutoffs", adder, ".Rda"))

    y = load(file=scut.filename)

    # stopifnot(all(names(all.cuts) == colnames(preds)))
    # stopifnot(all(names(all.pauc.cuts) == colnames(preds)))

    # stopifnot(all(names(all.scuts) == colnames(preds)))
    # stopifnot(all(names(all.spauc.cuts) == colnames(preds)))

    all.cuts = all.cuts[ colnames(preds)]
    all.scuts = all.scuts[ colnames(preds)]
    all.pauc.cuts = all.pauc.cuts[ colnames(preds)]
    all.spauc.cuts = all.spauc.cuts[ colnames(preds)]
    all.ssens.cuts = all.ssens.cuts[ colnames(preds)]
    all.sens.cuts = all.sens.cuts[ colnames(preds)]
    all.sdice.cuts = all.sdice.cuts[ colnames(preds)]
    all.dice.cuts = all.dice.cuts[ colnames(preds)]

    xpreds = preds

    for (type in types){
        if (type == ""){
            df$multiplier = 1
        }
        if (type == "_include"){
            df$multiplier = df$include
        }   
        if (type == "_include_all"){
            df$multiplier = df$include.all
        }
        if (type == "_zval"){
            df$multiplier = df$zval
        }       
        if (type == "_zval2"){
            df$multiplier = df$zval2
        }   
        if (type == "_zval_all"){
            df$multiplier = df$zval_all
        }               
        if (type == "_zval2_medztemp"){
          df$multiplier = df$zval2_medztemp
        }                   
        preds = xpreds
        #### only values 0 to 100 are likely to be ROI
        # preds = preds * df$in0100
        preds = preds * df$multiplier

        ################################
        # Running with variable slice thickness
        ################################

        run_vol = function(runpred){
            img = niftiarr(nim, runpred)
            img[is.na(img)] = 0
            # img@.Data[keep.ind] = runpred
            vol = get_roi_vol(img, dcmtables)$truevol 
        }

        # pred.imgs = vector(mode="list", 
            # length=ncol(preds))
        spreds = matrix(nrow=nrow(preds), ncol=ncol(preds))
        colnames(spreds) = colnames(preds)
        pb = txtProgressBar(max=ncol(preds), 
            style=3)
        ipred = 1
        for (ipred in seq(ncol(preds))){
            nd = preds[, ipred]
            img = nim
            nd[roi.not.in] = 0
            img = niftiarr(nim, nd)
            img = cal_img(img)
            img[is.na(img)]= 0
            img = datatyper(img, 
                datatype= convert.datatype()$FLOAT32, 
                bitpix = convert.bitpix()$FLOAT32)
            cn = colnames(preds)[ipred]
            # pred.imgs[[ipred]] = img
            outimg = nii.stub(fdf$img[get.pred], bn=TRUE)
            outimg = file.path(id.outdir, 
                paste0(outimg, "_", cn, adder)) 
            writeNIfTI(img, filename = outimg )

            outimg2 = paste0(outimg, 
                "_thresh")
            thresh = datatyper( 
                img > all.cuts[ipred]
                )
            writeNIfTI(thresh, 
                filename = outimg2 )              

            #### smooth the results
            sm.img  = mean_image(img, nvoxels = 1)
            # sm.img = local_moment(img, 
                # nvoxels = 1, moment =1,
            #   center = FALSE)[[1]]
            sm.img[abs(sm.img) < 
                .Machine$double.eps^0.5 ] = 0
            # sm.img = sm.img[keep.ind]
            spreds[, ipred] = sm.img
            # # spreds[roi.not.in, ipred] = 0

            sm.img = niftiarr(img, sm.img)
            # sm.img@.Data[keep.ind] = spreds[, ipred]
            # sm.img[is.na(sm.img)]= 0
            # sm.img[roi.not.in] = 0

            sm.img = cal_img(sm.img)
            sm.img  = datatyper(sm.img , 
                datatype= convert.datatype()$FLOAT32, 
                bitpix = convert.bitpix()$FLOAT32)

            outimg = paste0(outimg, "_smoothed", 
                adder, type)
            writeNIfTI(sm.img, filename = outimg )

            outimg2 = paste0(outimg, 
                "_thresh")
            thresh = datatyper( 
                sm.img > all.scuts[ipred]
                )
            writeNIfTI(thresh, 
                filename = outimg2 )  
            setTxtProgressBar(pb, ipred)

        }
        close(pb)
        
    }
    print(correct)
# }
# }
# }

# res = t(res)
# sres = t(sres)



