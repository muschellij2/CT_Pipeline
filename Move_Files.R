# dir <- "/Volumes/DATA/New_Age_Test/100-4001/"
dir <- "~/Desktop/101-307"
dir <- path.expand(dir)
setwd(dir)

movefiles <- function(files, indices, outdir, num=5, verb=TRUE){
	maxnum <- 10^num
	if (max(indices) > maxnum) stop("Need new format")
	for (i in 1:length(files)){
		ind <- indices[i]
		fmt <- paste0("%0", num, ".0f.dcm")
	  	new <- file.path(outdir, sprintf(fmt, ind))
	  	system(sprintf('mv "%s" "%s"', files[i], new))   
      if (verb) print(i)
	}
}

runformats <- function(fmts, startind= 0, indir, outdir){
	istart <- startind
	for (ifmt in 1:length(fmts)){
		fmt <- fmts[ifmt]
		print(paste0("Running format ", fmt))
		## get files
		x <- list.files(path=indir, pattern=fmt, recursive=TRUE, 
        	        full.names=TRUE)	
		if (length(x) > 0) {
			### get new indices
			runn <- (istart+1):(istart+1+length(x))
			movefiles(x, runn, outdir=outdir)
			## re-index the formats
			istart <- max(runn)
		}

	}
	return(istart)
}

fmts <- c(".dcm", "Image[0-9].*[0-9]$", "C[0-9].*[0-9]$", "C[0-9].*[A-Z]$")
lastind <- runformats(fmts, indir=dir, outdir=dir)


#sh ~/DHanley/CT_Registration/ICES/dcmsort_Final.sh -D /Volumes/DATA/New_Age_Test/265-389/ -o 
# /Volumes/DATA/New_Age_Test/265-389/Sorted -m
## THEN PIPELINE TEST

# x <- list.files(path=dir, pattern=".dcm", recursive=TRUE, 
#                 full.names=TRUE)
# runn <- 1:length(files)
# movefiles(x, runn, outdir=dir)

# ### oly crap
# xl <- max(runn)
# x <- list.files(path=dir, pattern="Image[0-9].*[0-9]$", recursive=TRUE, 
#                 full.names=TRUE)

# runn <- (xl+1):(xl+1+length(x))
# movefiles(x, runn, outdir=dir)

# xl <- max(runn)

# ## other format
# x <- list.files(path=dir, pattern="C[0-9].*[0-9]$", recursive=TRUE, 
#                 full.names=TRUE)

# runn <- (xl+1):(xl+1+length(x))
# movefiles(x, runn, outdir=dir)

