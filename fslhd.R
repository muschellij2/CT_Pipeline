fslhd <- function(file, intern=TRUE){
	cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
			"export FSLDIR FSLOUTPUTTYPE; echo $FSLDIR; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
	cmd <- paste(cmd, sprintf('/usr/local/fsl/bin/fslhd "%s"', file))
	system(cmd, intern=intern)
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


fslrange <- function(file, intern=TRUE){
  cmd <- paste("FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; ", 
               "export FSLDIR FSLOUTPUTTYPE; echo $FSLDIR; sh ${FSLDIR}/etc/fslconf/fsl.sh;")
  cmd <- paste(cmd, sprintf('/usr/local/fsl/bin/fslstats "%s"  -R', file))
  system(cmd, intern=intern)
}
