
fslsmooth <- function(file, mask=NULL, outfile=NULL, 
	sigma=10, intern=TRUE, local=FALSE){
	
	cmd <- NULL
	fsldir <- system("echo $FSLDIR", intern=TRUE)
	if (fsldir == "" | local) {
	  if (local) cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
	                          "export FSLDIR FSLOUTPUTTYPE; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
	}


	cmd <- paste(cmd, sprintf('fslmaths "%s"', file))
	if (! is.null(mask)) cmd <- paste(cmd, sprintf(' -mas "%s"', mask))
	if (is.null(outfile)) {
		outfile <- gsub("\\.gz$", "", file)
		outfile <- gsub("\\.nii$", "", outfile)
		outfile <- sprintf("%s_%s", outfile, sigma)
	}
	cmd <- paste(cmd, sprintf(' -s %s "%s";', sigma, outfile))

	### tempfile for mask.stub
	mask.stub <- basename(mask)
	mask.stub <- gsub("\\.gz$", "", mask.stub)
	mask.stub <- gsub("\\.nii$", "", mask.stub)
	mask.stub <- file.path(dirname(mask), mask.stub)
	mask.blur <- sprintf("%s_%s", mask.stub, sigma)
 	cmd <- paste(cmd, sprintf('fslmaths "%s" -s %s "%s";', 
 		mask, sigma, mask.blur))
 	cmd <- paste(cmd, sprintf('fslmaths "%s" -div "%s" -mas "%s" "%s";', 
 		outfile, mask.blur, mask, outfile))

# fslmaths  'FLAIRnorm_blur1_10.nii.gz' -div 'csfmask_blur1_10.nii.gz'   -mas 'csfmask.nii.gz'  'FLAIR_blur1_10_div.nii.gz'
	system(cmd, intern=intern)
	file.remove(paste0(mask.blur, ".nii.gz"))
}


fslmask <- function(file, mask=NULL, outfile=NULL, 
	unzip = FALSE, 
	intern=TRUE, 
	local=FALSE){
	
	cmd <- NULL
	fsldir <- system("echo $FSLDIR", intern=TRUE)
	niiftype = "FSLOUTPUTTYPE=NIFTI_GZ;"
	if (unzip) niiftype = "FSLOUTPUTTYPE=NIFTI;"
	if (fsldir == "" | local) {
	  if (local) cmd <- paste("FSLDIR=/usr/local/fsl; ", 
	                          "export FSLDIR; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
	}
	cmd = paste(cmd, niiftype)

	cmd <- paste(cmd, sprintf('fslmaths "%s" -mas "%s" "%s"', 
		file, mask, outfile))
# fslmaths  'FLAIRnorm_blur1_10.nii.gz' -div 'csfmask_blur1_10.nii.gz'   -mas 'csfmask.nii.gz'  'FLAIR_blur1_10_div.nii.gz'
	system(cmd, intern=intern)
	# if (unzip) {
	# 	cmd = sprintf('gunzip --force "%s".gz', outfile)
	# 	system(cmd, intern=intern)	
	# }
}


fslhd <- function(file, intern=TRUE, local=FALSE){
	cmd <- NULL
	fsldir <- system("echo $FSLDIR", intern=TRUE)
	if (fsldir == "" | local) {
	  if (local) cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
	                          "export FSLDIR FSLOUTPUTTYPE; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
	}  
	cmd <- paste(cmd, sprintf('fslhd "%s"', file))
	system(cmd, intern=intern)
}

fslhd.parse <- function(hd){
  ss <- strsplit(hd, split=" ")
  ss <- lapply(ss, function(x) x[!x %in% ""])
  ss <- lapply(ss, function(x){
    if (grepl("_xyz", x[1])) 
      x <- c(x[1], paste(x[2:length(x)], sep="", collapse= " "))
    if (grepl("form_name", x[1])) 
      x <- c(x[1], paste(x[2:length(x)], sep="", collapse= "-"))
    if (grepl("descrip", x[1])) 
      x <- c(x[1], paste(x[2:length(x)], sep="", collapse= " "))  
    return(x)
  })
  ss.len <- sapply(ss, length)
  ss <- ss[ss.len > 0]
  ss.len <- sapply(ss, length)
  stopifnot(all(ss.len %in% c(1,2)))
  ss <- lapply(ss, function(x){
    if (length(x) == 1) x <- c(x, NA)
    x
  })
  ss <- do.call("rbind", ss)
  rownames(ss) <- ss[,1]
  ss <- data.frame(value=ss[,2, drop=FALSE], stringsAsFactors=FALSE)
  return(ss)
}


getForms <- function(file){
	x <- fslhd(file)
	convmat <- function(form){
		ss <- strsplit(form, " ")
		ss <- t(sapply(ss, function(x) x[x!=""]))
		ss <- ss[, -1]
		class(ss) <- "numeric"
		return(ss)
	}
	sform <- x[grepl("sto_xyz:", x)]
	sform <- convmat(sform)
	qform <- x[grepl("qto_xyz:", x)]
	qform <- convmat(qform)

	sor <- x[grepl("sform_(x|y|z)orient", x)]
	qor <- x[grepl("qform_(x|y|z)orient", x)]

	short_orient <- function(orient){
		ss <- strsplit(orient, " ")
		ss <- sapply(ss, function(x) x[x!=""])[2,]
		first <- substr(ss, 1,1)
		ss2 <- strsplit(ss, "-")
		ss2 <- sapply(ss2, function(x) x[length(x)])
		second <- substr(ss2, 1,1)
		paste(first, second, sep="")
	}
	ssor <- short_orient(sor)
	sqor <- short_orient(qor)
  
  sform_code <- as.numeric(
    gsub("sform_code", "", x[grepl("sform_code", x)])
  )

	qform_code <- as.numeric(
	  gsub("qform_code", "", x[grepl("qform_code", x)])
	)  
	return(list(qform=qform, sform=sform, sor=sor, qor=qor, 
		ssor=ssor, sqor=sqor, sform_code= sform_code, qform_code=qform_code ))
}

checkout <- function(hd){
	det.equal <- sign(det(hd$sform)) == sign(det(hd$qform))
	lr.equal <- hd$ssor[1] == hd$sqor[1]

	if ( (det.equal & !lr.equal) | (!det.equal & lr.equal)) {
		return(FALSE)
	}
	return(TRUE)
}

check_file <- function(file){
	hd <- getForms(file)
	checkout(hd)
}

check_sform <- function(hd, value=0){
	hd$sform_code == value
}


check_sform_file <- function(file, value=0){
	hd <- getForms(file)
	check_sform(hd, value=value)
}
## if sign(det(res$sform)) == sign(det(res$qform)) and res$ssor[1] == res$sqor[1] then all good



fslrange <- function(file, intern=TRUE, local=FALSE){
	cmd <- NULL
	fsldir <- system("echo $FSLDIR", intern=TRUE)
	if (fsldir == "" | local) {
	  if (local) cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
	                          "export FSLDIR FSLOUTPUTTYPE; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
	}
	cmd <- paste(cmd, sprintf('fslstats "%s"  -R', file))
  system(cmd, intern=intern)
}


fslfill = function(file, outfile = file, bin=TRUE, intern=TRUE, 
                   local=FALSE){
  cmd <- NULL
  fsldir <- system("echo $FSLDIR", intern=TRUE)
  if (fsldir == "" | local) {
    if (local) cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
                            "export FSLDIR FSLOUTPUTTYPE; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
  }
  runbin = ""
  if (bin) runbin = "-bin"
  cmd <- paste(cmd, sprintf('fslmaths "%s" %s -fillh "%s"', file, runbin, outfile))
  system(cmd, intern=intern)
}


fslthresh = function(file, outfile = file, thresh = 0, intern=TRUE, 
                   local=FALSE, unzip=FALSE){
  cmd <- NULL
  fsldir <- system("echo $FSLDIR", intern=TRUE)
	niiftype = "FSLOUTPUTTYPE=NIFTI_GZ;"
	if (unzip) niiftype = "FSLOUTPUTTYPE=NIFTI;"
	if (fsldir == "" | local) {
	  if (local) cmd <- paste("FSLDIR=/usr/local/fsl; ", 
	                          "export FSLDIR; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
	}
	cmd = paste(cmd, niiftype)
  cmd <- paste(cmd, sprintf('fslmaths "%s" -thr %f "%s"', 
  	file, thresh, outfile))
  system(cmd, intern=intern)
}

fslsub2 = function(file, outfile = file, intern=TRUE, 
                   local=FALSE, unzip=FALSE){
  cmd <- NULL
  fsldir <- system("echo $FSLDIR", intern=TRUE)
	niiftype = "FSLOUTPUTTYPE=NIFTI_GZ;"
	if (unzip) niiftype = "FSLOUTPUTTYPE=NIFTI;"
	if (fsldir == "" | local) {
	  if (local) cmd <- paste("FSLDIR=/usr/local/fsl; ", 
	                          "export FSLDIR; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
	}
	cmd = paste(cmd, niiftype)
  cmd <- paste(cmd, sprintf('fslmaths "%s" -subsamp2 "%s"', file, outfile))
  system(cmd, intern=intern)
}
