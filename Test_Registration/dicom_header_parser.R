
make.numeric <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  oldx <- x
  x[x == ""] <- NA
  if (all(is.na(x))) return(oldx)
  ## if they all are real - turn to numeric
  suppressWarnings(num_x <- as.numeric(x))
  keep <- !(is.na(x) & is.na(num_x))
  ### if all numeric 
  if (any(is.na(num_x[keep]))) {
    return(oldx)
  } else {
    return(num_x)
  }
}

dicom_header_parse = function(hdr, tryConvert=TRUE){
  nfiles = length(hdr)
  all.names = unique(unlist(lapply(hdr, function(x) x[, "name"])))
  n.names = length(all.names)
  mat = matrix(NA, nrow=nfiles, ncol=n.names)
  colnames(mat) = all.names

  for (ifile in seq(nfiles)){
  	x = hdr[[ifile]]
  	val = x[, c("value")] 
  	nn = names(val) = x[, "name"]
  	mat[ifile, nn] = val[nn] 
  }

  mat = data.frame(mat, stringsAsFactors=FALSE)
  if (tryConvert) {
  	for (icol in seq(n.names)){
  		mat[,icol] = make.numeric(mat[, icol]) 
  	}
  }
  rownames(mat) = names(hdr)
  return(mat)
}
