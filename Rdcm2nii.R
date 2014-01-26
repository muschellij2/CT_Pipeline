Rdcm2nii <- function(basedir, sortdir, verbose=TRUE, ...){
 
     gf <- getfiles(basedir)
   files <- gf$files
   paths <- gf$paths
 
     ## THEN PIPELINE TEST
     ipath <- 1
   lpath <- length(paths)
 
     ### need dcm2nii in your path! - from MRICRON
     if (lpath >= 1){
         if (verbose) cat("Converting to nii using R \n")
     
           for (ipath in 1:length(paths)){
               path <- paths[ipath]
               res <- system(sprintf('rm "%s"/*.nii.gz', path), ignore.stdout = !verbose)
               
                 dcm2nii_worker(path)
         
                 niis <- dir(path=path, pattern=".nii.gz")
               stub <- basename(path)
         
                 iddir <- file.path(basedir)
               name <- stub
               
                 ### copy dicom header info
                 make_txt(path)
               
                 if (length(niis) > 1){    
                     # stop("it")
                       # cmd <- 'FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; export FSLDIR FSLOUTPUTTYPE; '
                       # cmd <- paste0(cmd, 'echo $FSLDIR; sh ${FSLDIR}/etc/fslconf/fsl.sh; ', 
                       #   '/usr/local/fsl/bin/fslmerge -z "%s"/"%s".nii.gz "%s"/*.nii.gz')
                       cmd <- c('fslmerge -z "%s"/"%s".nii.gz "%s"/*.nii.gz')
                     system(sprintf(cmd, iddir, name, path), ignore.stdout = !verbose)
                       
                     }
               if (length(niis) == 1){
                   system(sprintf('mv "%s"/*.nii.gz "%s"/"%s".nii.gz', path, iddir, name), 
                                     ignore.stdout = !verbose)
                 }
               system(sprintf('rm "%s"/*.nii.gz', path), ignore.stdout = !verbose)
               setwd(dirname(path))
               #     stop("me")
                 newpath <- file.path(dirname(path), name)
               system(sprintf('mv "%s" "%s"', path, newpath), ignore.stdout = !verbose)
         
                 # setwd(basename(newpath))
                 # tarfile <- paste(newpath, ".tar.gz", sep="")
                 # tar(tarfile, compression="gzip")
                 system(sprintf('tar -czf "%s" ./"%s"', paste(newpath, ".tar.gz", sep=""), 
                                           basename(newpath)), ignore.stdout = !verbose)
               
                 system(sprintf('rm -R "%s"', newpath), ignore.stdout = !verbose)
             }
       } # end loop over paths
 
   } ## end Rdcm2nii

