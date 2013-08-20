fslthresh <- function(image, autoname=TRUE, outdir=NULL, outfile = NULL, 
                      lower=NULL, upper=NULL, intern=FALSE, 
                      mask=FALSE, verbose=TRUE){
  stopifnot(file.exists(image))
  stopifnot(!( is.null(lower) & is.null(upper) ))
  
  ## check to see if installed on sytem
  fslcmd <- makeFSLcmd("fslmaths")
  
  ### make the command
  cmd <- sprintf('%s "%s"', fslcmd, image)
  
  if (!is.null(lower)) {
    stopifnot(is.numeric(lower))
    cmd <- paste0(cmd, " -thr ", lower)
  }
  if (!is.null(upper)) {
    stopifnot(is.numeric(upper))
    cmd <- paste0(cmd, " -uthr ", upper) 
  }
  
  ## autonaming !
  if (autoname & !is.null(outfile)){
    warning("Autoname and outfile specified, outfile used")
  } 
  
  if (autoname & is.null(outfile)){
    if (!is.null(lower) & is.null(upper)){
      addstub <- paste0("L", lower)
    } else if (is.null(lower) & !is.null(upper)) {
      addstub <- paste0("U", upper)
    } else if (!is.null(lower) & !is.null(upper)) {
      addstub <- paste0("L", lower, "_U", upper)
    }  
    ## strip off .nii.gz or .nii
    img <- gsub("\\.gz$", "", image, fixed=FALSE)
    img <- gsub("\\.nii$", "", img, fixed=FALSE)
    
    stub <- basename(img)
    if (mask) addstub <- paste0(addstub, "_mask")
    outdir <- ifelse(is.null(outdir), dirname(image), outdir)
    outfile <- file.path(outdir, paste0(stub, "_", addstub))
  }
  stopifnot(!is.null(outfile))
  outfile <- gsub("\\.gz$", "", outfile, fixed=FALSE)
  outfile <- gsub("\\.nii$", "", outfile, fixed=FALSE)
  
  if (mask) cmd <- paste0(cmd, " -bin")  
  cmd <- paste0(cmd, ' "', outfile, '"')
#   cat(cmd, sep="\n")
  x <- system(cmd, ignore.stdout = !verbose, intern=intern)
  return(x)
}

makeFSLcmd <- function(fslcmd){
  ## check to see if installed on sytem
  x <- system(paste0("which ", fslcmd), ignore.stdout=TRUE)
  if (x > 0){
    cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
                 "export FSLDIR FSLOUTPUTTYPE; echo $FSLDIR;", 
                 "sh ${FSLDIR}/etc/fslconf/fsl.sh;")
    fslcmd <- paste0(cmd, " /usr/local/fsl/bin/", fslcmd )
  } 
  return(fslcmd)
}

flirt <- function(image, ref, outfile, outmat = NULL, rigid = FALSE, cost="corratio", 
                  searchcost=cost, verbose=TRUE, run=TRUE, addopts=""){

  stopifnot(file.exists(image))
  stopifnot(file.exists(ref))

  outfile <- gsub("\\.gz$", "", outfile, fixed=FALSE)
  outfile <- gsub("\\.nii$", "", outfile, fixed=FALSE)
  
  cmd <- makeFSLcmd("flirt")
  cmd <- paste0(cmd, sprintf(' -ref "%s" -in "%s"', ref, image))
  cmd <- paste0(cmd, sprintf(' -out "%s"', outfile))
  if (rigid) cmd <- paste0(cmd, " -dof 6")
  cmd <- paste0(cmd, " -cost ", cost)
  cmd <- paste0(cmd, " -searchcost ", searchcost)
  if (verbose) cmd <- paste0(cmd, " -v")
  if (!is.null(outmat)) cmd <- paste0(cmd, sprintf(' -omat "%s"', outmat))
  cmd <- paste(cmd, addopts, sep= " ")
  if (run) {
    system(cmd)
  } else {
    cat(cmd)
  }
}

flirt.wrap <- function(image, rigid=TRUE, cost="mutualinfo", 
                       searchcost="mutualinfo", run=FALSE, mask=FALSE, 
                       ...){
  img <- gsub("\\.gz$", "", image, fixed=FALSE)
  img <- gsub("\\.nii$", "", img, fixed=FALSE)
  stub <- basename(img)
  path <- dirname(img)
  if (mask) stub <- paste0(stub, "_mask")
  
  outfile <- file.path(path, paste0("rigid_", stub))
  outmat <- paste0(outfile, ".mat")
  flirt(image = image, outfile=outfile, outmat = outmat, rigid=rigid, 
        cost=cost, searchcost=searchcost, run=run, ... )
}