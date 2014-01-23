###############################################################
# Obsolete - jsut smooth with sigma = 3 on mask then binarize
###############################################################

library(oro.dicom)
library(plyr)
mask <- readNIfTI("./Skull_Stripped/205-509_20100418_0000_CT_2_HEADBRAIN_WO_CONTRAST_ROUTINE_SS_No1024_Mask_0.1.nii.gz", reorient=FALSE)

img <- readNIfTI("205-509_20100418_0000_CT_2_HEADBRAIN_WO_CONTRAST_ROUTINE.nii.gz", reorient=FALSE)

dilate <- function(mask){
	inds <- which(mask > 0, arr.ind=TRUE)
	inds <- data.frame(inds)
	for (iz in unique(inds$dim3)){
		ind <- inds[inds$dim3 == iz, ]
		for (iy in unique(ind$dim2)){
			ind2 <- ind[ind$dim2 == iy, ]
		}
	}

}