## get the basename for a file, for example getBase("blah.nii.gz", 2) = "blah"
getBase <- function(x, ind=1){
  sapply(strsplit(x, split="\\."), function(xx) 
    paste(xx[1:(length(xx)-ind)], collapse=".", sep=""))
}
