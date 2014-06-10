#################################
# Save Voxel information for reporting
#################################
rm(list=ls())
library(oro.nifti)
library(plyr)
library(limma)
library(microbenchmark)
library(abind)
library(ggplot2)
library(data.table)
homedir = "/Applications"
rootdir = "/Volumes/DATA_LOCAL/Image_Processing"
if (Sys.info()[["user"]] %in% "jmuschel") {
  homedir = "~"
  rootdir = "/dexter/disk2/smart/stroke_ct/ident"
}
progdir = file.path(rootdir, "programs")
source(file.path(progdir, "convert_DICOM.R"))
source(file.path(progdir, "fslhd.R"))
basedir = file.path(rootdir, "Registration")
tempdir = file.path(rootdir, "Template")
atlasdir = file.path(tempdir, "atlases")

outdir = file.path(basedir, "results")

whichdir = "reoriented"

#### load voxel data
outfile = file.path(outdir, "Voxel_Matrix.Rda")
load(file=outfile )

fnames = colnames(mat)

#### keeping if over 10 people have ICH in that locaiton
ncut = 10
all.nvox = sum(rs > 0)
mat = mat[rs > ncut, ]

#### getting unique rows and duplications
dups = duplicated(mat)
# cn = paste0("V", 1:ncol(mat))
cn = colnames(mat)
mm = cbind(mat, id = 1:nrow(mat))
dt = data.table(mm, key=cn)
urows = unique(dt)
urows[, c('group') := 1:nrow(urows)]
urows[, c('id') := NULL]
xx = dt[urows]

group = xx[, list(id, group)]
setkey(group, 'id')

nuniq.rows = nrow(urows)
nvox = nrow(mat)
where = "Save_Voxel_Info.R"

output = file.path(outdir, "Voxel_Info.Rda")
save(nvox, ncut, nuniq.rows, 
	dups, group, where, all.nvox,
	fnames, 
	file=output)
