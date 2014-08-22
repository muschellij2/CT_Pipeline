rm(list=ls())
library(cttools)
library(fslr)
img = readNIfTI("~/Dropbox/Temp_Lung/NIfTI/113352678_20100806_1357_CT_2_PET^PETCT_AC_WB_CPT_AC_CT.nii.gz")
newnifti = function(img, arr){
    x = img
    x@.Data = arr
    x = cal_img(x)
}


thresh1 = newnifti(img, img <= -150 & img >= -800)
orthographic(cal_img(thresh1))

smooth = fslsmooth(thresh1, sigma=10, retimg=TRUE)
# orthographic(cal_img(smooth > .3))

mask = smooth > 0.3
### could do smooth > 0.3 and img <= -150 
lung = img
lung[mask == 0] = -1024
orthographic(lung)

lung2 = lung
lung2[lung2 > -150 | lung2 < -800] = -1024
lung2 = cal_img(lung2)

mask = cal_img(lung != -1024)
xyz = round(cog(mask))
mask.overlay(img, mask, xyz=xyz)
