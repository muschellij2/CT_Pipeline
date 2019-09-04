#######################################
# Need to run ROI_Copying_Reorienting.R before
###################################
rm(list=ls())
library(fslr)
library(cttools)
library(plyr)
library(scales)
library(reshape2)
Sys.setenv(FSLOUTPUTTYPE = "NIFTI")

#### delete all ROI files
### find . -regextype posix-extended -regex "^./[0-9].*[0-9]$"
###  -exec rm -r {} \;

  username <- Sys.info()["user"][[1]]

  cluster=FALSE
  if (username == "muschellij2"){
    # rootdir <- "/Volumes/DATA/New_Age_Test"
    rootdir <- "~/CT_Registration"
  } else {
    rootdir <- "/dexter/disk2/smart/stroke_ct/ident"
    cluster =TRUE;
  }
    rootdir <- path.expand(rootdir)

  # ss <- as.numeric(strsplit(id, "-")[[1]][2])
  # if (ss > 4000){
  #   study <- "CLEAR_III"
  #   dpath <- file.path("CLEAR", "CLEAR III")
  # } else if (ss > 300 & ss < 500){
  #   dpath <- study <- "MISTIE"
  # } else if (ss > 500 & ss < 4000) {
  #   dpath <- study <- "ICES" 
  # }


  rootdir <<- path.expand(rootdir)

#progdir <- file.path(dirname(basedir), "programs")
  progdir <- file.path(rootdir, "programs")
  basedir <<- file.path(rootdir, "Registration")
  tempdir = file.path(rootdir, "Template")
  atlasdir = file.path(tempdir, "atlases")  



#### setting up if things are on the cluster or not

# load(file.path(rootdir, "Registration", 
#   "Registration_Image_Names.Rda"))
load(file.path(rootdir, "Registration", 
  "All_Registration_Image_Names.Rda"))  
# df = df[grep("101-307|101-308|102-317|102-322", df$roi.nii), ]


ids.111 = read.csv(file.path(basedir, "111_patients.csv"), 
  stringsAsFactors= FALSE)
uid = ids.111$patientName
all.ids = ids.111$id

############################
# Take out 111 already ran for paper - this is for testing
############################
xdf = df
reorient = TRUE
normalize = TRUE
subsamp = TRUE
rerun = TRUE

bad.imgs = c(233, 275, 293, 312, 345, 358, 359, 362, 365, 366, 367, 
  370, 
375, 388, 421, 424, 428, 431, 444, 595, 1055, 1148, 1149, 1150
)

# iimg <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# if (is.na(iimg)) iimg = 123
iimg = bad.imgs
df = df[iimg, , drop=FALSE]


imgs = mlply(.fun = function(outfile, roi.nii, raw, ss){
  c(outfile, roi.nii, raw, ss)
}, .data=df[, c("outfile", "roi.nii", "raw", "ss")])
attr(imgs, "split_labels") <- NULL

x = imgs[[1]]

if (reorient){
  rets = laply(.data=imgs, .fun = function(x){
    acpc_reorient(infiles = x)
  }, .progress="text")
  stopifnot(all(rets == 0))
  rois = df$roi.nii
  l_ply(.data=rois, .fun = function(x) {
    fslbin(x, outfile=x)
  }, .progress = "text")
  # acpc_reorient(infiles = x)
}


df$out.roi.nii = file.path(dirname(df$roi.nii), 
  paste0("bws", basename(df$roi.nii)))
df$out.raw = file.path(dirname(df$raw), 
  paste0("w", basename(df$raw)))

### create list of data
imgs = mlply(.fun = function( raw, roi.nii, out.raw, 
  out.roi.nii, ss){
  c(raw=raw, roi.nii=roi.nii, 
    out.raw=out.raw, out.roi.nii=out.roi.nii, ss = ss)
}, .data=df[, c("raw", "roi.nii", "out.raw", "out.roi.nii",
  "ss")])
attr(imgs, "split_labels") <- NULL

x = imgs[[1]]
## 131-310
if (normalize){
  ### matlab problems if I don't change directories
    # - because dicomread and spm_affreg in there
  # gwd = getwd()
  setwd("~/")
  ### delete previous iterations of normalization
  dim.equal = laply(.data=imgs, .fun = function(x){
      hd = fslhd(x['raw'])
      dims = fslhd.parse(hd)[paste0("dim", 1:3),]
      hd = fslhd(x['roi.nii'])
      roi.dims = fslhd.parse(hd)[paste0("dim", 1:3),]      
      all(dims == roi.dims)
    }, .progress="text")

  rets = laply(.data=imgs, .fun = function(x){
      suppressWarnings(y <- 
        file.remove(x[c("out.raw", "out.roi.nii")])
        )
      r = run_ctnorm(rawfile = x["raw"], roifile = x["roi.nii"], 
        deleteinter = TRUE)    
      return(r)
    }, .progress="text")
  ### after orientation
  print(which(rets != 0))
  cat("Normalizing the SS data \n")
  rets = laply(.data=imgs, .fun = function(x){
      rfile = x["raw"]
      fol = dirname(rfile)
      sn_mat = paste0('c', basename(nii.stub(rfile)), "_sn.mat")  
      sn_mat = file.path(fol, sn_mat)
      r2 = run_ctnorm_write(sn_mat=sn_mat, files= x["ss"])
  }, .progress = "text")
  # setwd(gwd)

}


