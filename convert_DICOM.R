require(stringr)
require(plyr)
require(oro.nifti)
require(oro.dicom)



#### wrapper to convert an entire directory to sort/move/nifti
convert_DICOM <- function(basedir, progdir, verbose=TRUE, 
  untar=FALSE, useR = TRUE, id = NULL, ...){
  
  setwd(basedir)

  sortdir <- file.path(basedir, "Sorted")
  if (!file.exists(sortdir)) system(sprintf('mkdir -p "%s"', sortdir))

  if (untar){
    ### unzip all the tar.gz (if you have to redo)
    if (verbose) cat("Un tarballing data \n")
    files <- dir(path=basedir, pattern=".tar.gz$", full.names=TRUE, 
          recursive=TRUE)
    gants <- grepl("_ungantry", files)

    ## don't need to untarball the gantry corrected -these willbe over wrtiten
    if (any(gants)){
        gants <- files[gants]
        rm.files <- gsub("_ungantry", "", gants)
        files <- files[! files %in% rm.files]
    }
    ### untarball the files using R's command
    if (length(files) > 0){
      for (ifile in files){
        dname <- dirname(ifile)
        untar(ifile, exdir = dname, compressed='gzip')
        if (verbose) print(ifile)
      }
    }
    files <- dir(path=basedir, pattern=".nii.gz", 
          full.names=TRUE, recursive=TRUE)
    file.remove(files)
  }


  ### Moving files into one big directory
  ### watch out for / ending in basedir
  onedir(basedir, verbose, ...)

  ## putting into respective folders using dcmdump
  if (useR) {
    # if (verbose) cat("dcm sorting")
    dcmtables = Rdcmsort(basedir, sortdir, id = id)
  } else {
    dcmsort(basedir, progdir, sortdir, verbose, ...)
  }

  if (useR & inherits(dcmtables, "logical")){
    return(FALSE)
  }
  ## gantry tilt correction
  file_gc(basedir, progdir, verbose, ...)
  
  ## conversion
  if (useR) {
    ret = Rdcm2nii(basedir, sortdir, verbose=TRUE, ...)
  } else {
    dcm2nii(basedir, progdir, sortdir, verbose, ...)
  }
  
  gf <- getfiles(basedir)

  if (verbose) cat("Deleting Empty Directories \n")
  expaths <- list.dirs(basedir, recursive=TRUE, full.names=TRUE)
  expaths <- expaths[!(expaths %in% c(sortdir, basedir))]
  expaths <- expaths[!(expaths %in% gf$paths)]
  expaths <- expaths[!grepl("plots|Skull_Stripped|Registered|RawNIfTI|reoriented|FLIRT", 
    expaths)]
  for (ipath in expaths) system(sprintf('rmdir "%s"', ipath))
  for (ipath in expaths) system(sprintf('rmdir "%s"', ipath))

  return(TRUE)
} ###end function



movefiles <- function(files, indices, outdir, num=6, verbose=TRUE){
  maxnum <- 10^num
  if (max(indices) > maxnum) stop("Need new format")
  for (i in 1:length(files)){
    ind <- indices[i]
    fmt <- paste0("%0", num, ".0f.dcm")
      new <- file.path(outdir, sprintf(fmt, ind))
      system(sprintf('mv "%s" "%s"', files[i], new))   
      if (verbose) print(i)
  }
}

## find /thisdir -type f -name '*.ogg' -exec mv -i {} /somedir \;


runformats <- function(fmts, startind= 0, indir, outdir, verbose=TRUE){
  istart <- startind
  for (ifmt in 1:length(fmts)){
    fmt <- fmts[ifmt]
    if (verbose) print(paste0("Running format ", fmt))
    ## get files
    x <- list.files(path=indir, pattern=fmt, recursive=TRUE, 
                  full.names=TRUE)  
    if (length(x) > 0) {
      ### get new indices
      runn <- (istart+1):(istart+1+length(x))
      movefiles(x, runn, outdir=outdir, verbose=verbose)
      ## re-index the formats
      istart <- max(runn)
    }

  }
  return(istart)
}

  getfiles <- function(basedir){
    files <- dir(path=basedir, pattern="dcm|DCM", full.names=TRUE, 
          recursive=TRUE)
    paths <- unique(dirname(files))
    return(list(files=files, paths=paths))
  }

onedir <- function(basedir, verbose=TRUE, ...){
    gf <- getfiles(basedir)


  ### Moving files into one big directory
    ### watch out for / ending in basedir
  if (!all(gf$paths == basedir)){
    fmts <- c(".dcm", ".DCM", "Image[0-9].*[0-9]$", 
      "C[0-9].*[0-9]$", "C[0-9].*[A-Z]$")
    lastind <- runformats(fmts, indir=basedir, 
      outdir=basedir, verbose=verbose)
  }
}

dcmsort <- function(basedir, progdir, sortdir, dcmsortopt = '', verbose=TRUE, ...){
  gf <- getfiles(basedir)


  ### need dcmdump in your path - from DCMTK
  ### sort them out
  if (length(gf$paths) > 0){
    if (all(gf$paths == basedir)){
      if (verbose) cat("Sorting DICOMs \n")
      cmd <- sprintf('sh "%s"/dcmsort_Final.sh -D "%s" -o "%s" -m -x %s', progdir,
                     basedir, sortdir, dcmsortopt)
      system(cmd)
    } else {
      stop("something is off - need to all have dcms in one file")
    }
  }
}

### extract from readDICOM header
extract.from.hdr = function(hdr, key, numeric=FALSE){
  key <- gsub("(", "", key, fixed=TRUE)
  key <- gsub(")", "", key, fixed=TRUE)
  ss = strsplit(key, ",")
  groups = sapply(ss, function(x) x[1])
  elements = sapply(ss, function(x) x[2])

  crit = hdr$group == groups & hdr$element == elements
  values = hdr[crit, "value"]
  values = str_trim(values)
  values = gsub("\\s\\s+" , " ", values)
  if (numeric) values = as.numeric(values)
  return (values)
}



extract.desc = function(hdr, key){
  desc = extract.from.hdr(hdr, key)
  desc = gsub(" ", "_", desc)
  desc = ifelse(is.na(desc) | length(desc) == 0, "unnamed", desc)
  desc = gsub("/", "", desc); 
  desc = gsub("\47", "", desc)
  return(desc)
} 

extract.time = function(hdr, key){
  runtime     = extract.from.hdr(hdr, key, numeric=TRUE)  
  NUMBER      = runtime / 100
  NUMBER      = sprintf("%04.0f", NUMBER)
  NUMBER        = ifelse(!is.na(NUMBER) & !(NUMBER %in% "") & 
    length(NUMBER) > 0, NUMBER, NA)  
}

miss.data = function(info){
  if (is.null(info)) return(TRUE)
  if (is.na(info)) return(TRUE)
  if ( length(info) == 0 ) return(TRUE)
  if ( info %in% "" ) return(TRUE)
  return(FALSE)
}

name.file = function(hdr, id = NULL){
  # print(colnames(hdr))
  PID         = extract.from.hdr(hdr, '0010,0020)')
  # print(PID)
  # print(id)
  PID         = ifelse(miss.data(PID), id, PID)

  StudyDate   = extract.from.hdr(hdr, '(0008,0020)', numeric=TRUE)
  NUMBER      = extract.time(hdr, '(0008,0030)')

  StudyDesc  = extract.desc(hdr, '(0008,1030)')
  SeriesDesc  = extract.desc(hdr, '(0008,103E)')

  SeriesDate  = extract.from.hdr(hdr, '(0008,0021)', numeric=TRUE)
  SENUMBER      = extract.time(hdr, '(0008,0031)')

  ACQNUM       = extract.time(hdr, '(0008,0032)')
  CNUM         = extract.time(hdr, '(0008,0033)')

  SID         = extract.from.hdr(hdr, '(0020,000D)')
  #  echo $SeriesDesc
  Modality    = extract.from.hdr(hdr, '(0008,0060)')
  # (0020,000d)
  UUID        = extract.from.hdr(hdr, '(0008,0018)')
  ss          = strsplit(UUID, "\\.")
  UUID        = sapply(ss, function(x) paste(x[1:(length(x)-2)], sep="", collapse="."))
    
  ### series number
  SNUM        = extract.from.hdr(hdr, '(0020,0011)')

  FNUM        = ifelse(is.na(NUMBER), SENUMBER, NUMBER)
  FNUM        = ifelse( miss.data(FNUM), ACQNUM, FNUM)
  FNUM        = ifelse( miss.data(FNUM), CNUM, FNUM)

  DATE        = ifelse(!is.na(StudyDate) & !(StudyDate %in% "") & length(StudyDate) > 0, 
    StudyDate, SeriesDate)
  DATE        = ifelse(!is.na(DATE) & !(DATE %in% "") & length(DATE) > 0, DATE, NA)

  DATER       = paste0(DATE, "_", FNUM)

  PNAME       = paste(PID, DATER, Modality, SNUM, StudyDesc, SeriesDesc, sep="_")
  PNAME       = gsub(" ", "_", PNAME)
  return(PNAME)
}


Rdcmsort = function(basedir, sortdir, id = NULL, 
  writeFile=FALSE, verbose = TRUE){
  dcms = getfiles(basedir)$files

  if (length(dcms) > 0){
    # ifile = 1;
    ### read in EVERY HEADER from this 
    if (verbose) cat("Reading Headers \n")
    hdrl = rereadDICOMHeader(dcms)

    names(hdrl) = dcms


    # hdr = hdrl[[length(dcms)]]
    if (verbose) cat("Making filenames \n")
    filenames = llply(hdrl, name.file, id = id,
      .progress="text")
    flen = sapply(filenames, length)
    over1 = flen > 1
    if (any(over1)){
      ind = which(over1)
      print("Multi Names")
      print(filenames[ind])
      filenames = llply(filenames, function(x) x[1])
    }
    filenames = unlist(filenames)
    names(filenames)= NULL

    basenames = basename(dcms)

    new.dirs = file.path(sortdir, filenames)
    x = sapply(unique(new.dirs), dir.create, showWarnings=FALSE)
    
    new.fnames = file.path(new.dirs, basenames)

    x = file.rename(dcms, new.fnames)
    stopifnot(all(x))

    # hdrl = lapply(hdrl, function(df) {
    #   df$uval = paste(df$group, df$element, df$name, sep="-")
    #    df$duplicated = duplicated(df$uval)
    #    return(df)
    #   })
    
    dcmtables = dicomTable(hdrl)
    rownames(dcmtables) = new.fnames

    if (writeFile){
      save(dcmtables, 
        file=file.path(basedir, "All_Header_Info.Rda"))
    }

    return(dcmtables)
  } else {
    return(FALSE)
  }
  # names(hdrl) = new.fnames
}

dcm2nii <- function(basedir, progdir, sortdir, verbose=TRUE, 
                    dcm2niicmd = "dcm2nii", ...){

  gf <- getfiles(basedir)
  files <- gf$files
  paths <- gf$paths

  ## THEN PIPELINE TEST
  ipath <- 1
  lpath <- length(paths)

  ### need dcm2nii in your path! - from MRICRON
  if (lpath >= 1){
    if (verbose) cat("Converting to nii \n")

    for (ipath in 1:length(paths)){
      path <- paths[ipath]
      x = system(sprintf('rm "%s"/*.nii.gz', path), intern=TRUE)
      intern=TRUE
      res <- system(sprintf('%s -b "%s"/CT_dcm2nii.ini "%s"', dcm2niicmd,
          progdir, path), intern=TRUE)
      if (intern) {
        errs <- any(grepl("Error", res))
      } else {
        errs <- res != 0
      }
      stopifnot(length(errs) == 1)
      if (  errs ){
              system(sprintf('rm "%s"/*.nii.gz', path))
              print("Error in DCM2NII")
              next
      }        
      niis <- dir(path=path, pattern=".nii.gz")
      stub <- basename(path)

      iddir <- file.path(basedir)
      name <- stub
      
      ### copy dicom header info
      header_txt(path)
#       dcm <- dir(path=path, pattern="dcm|DCM", full.names=TRUE, recursive=FALSE)[1]    
#       txtfile <- file.path(dirname(path), paste0(stub, ".txt"))
#       system(sprintf('dcmdump "%s" > "%s"', dcm, txtfile))
      
      
      if (length(niis) > 1){    
        # stop("it")
  # cmd <- 'FSLDIR=/usr/local/fsl; FSLOUTPUTTYPE=NIFTI_GZ; export FSLDIR FSLOUTPUTTYPE; '
  # cmd <- paste0(cmd, 'echo $FSLDIR; sh ${FSLDIR}/etc/fslconf/fsl.sh; ', 
  #   '/usr/local/fsl/bin/fslmerge -z "%s"/"%s".nii.gz "%s"/*.nii.gz')
        cmd <- c('fslmerge -z "%s"/"%s".nii.gz "%s"/*.nii.gz')
        system(sprintf(cmd, iddir, name, path))
          
      }
      if (length(niis) == 1){
        system(sprintf('mv "%s"/*.nii.gz "%s"/"%s".nii.gz', 
          path, iddir, name))
      }
      system(sprintf('rm "%s"/*.nii.gz', path))
      setwd(dirname(path))
      #     stop("me")
      newpath <- file.path(dirname(path), name)
      system(sprintf('mv "%s" "%s"', path, newpath))

      # setwd(basename(newpath))
      # tarfile <- paste(newpath, ".tar.gz", sep="")
      # tar(tarfile, compression="gzip")
      system(sprintf('tar -czf "%s" ./"%s"', paste(newpath, ".tar.gz", sep=""), 
          basename(newpath)))
      
      system(sprintf('rm -R "%s"', newpath))
    }
  } # end loop over paths

} ## end dcm2nii





Skull_Strip <- function(basedir, progdir, CTonly=TRUE, 
                      dropstring=NULL, opts = "-f 0.1 -b", 
                      verbose=TRUE){
  outdir <- file.path(basedir, "Skull_Stripped")
  if (!file.exists(outdir)) system(sprintf('mkdir -p "%s"', outdir))

  niis <- dir(path=basedir, pattern=".nii.gz", 
    full.names=TRUE, recursive=FALSE)
  
  if (CTonly) niis <- niis[grepl("_CT_", niis)]
  ## drop out scans (like CTA)
  if (!is.null(dropstring) & length(dropstring) > 0){
    for (istring in seq_along(dropstring)){
      niis <- niis[!grepl(dropstring[istring], niis)]
    }
  }
  ### skull strip  data
  inii <- niis[1]
  for (inii in niis){
    Skull_Strip_file(img=inii, progdir=progdir, 
      outdir=outdir, opts=opts, verbose=verbose)
  }

}

Skull_Strip_file <- function(img, progdir, outdir, opts = "", 
  verbose=TRUE){
    cmd <- sprintf('sh "%s"/Brain_Seg_Function.sh -i "%s" -o "%s" %s', 
      progdir, img, outdir, opts) 
    res <- system(cmd)
    if (verbose) print(img)
    return(res)  
}



## get key from header 
extract <- function(file, key){
  key <- gsub("(", "\\(", key, fixed=TRUE)
  key <- gsub(")", "\\)", key, fixed=TRUE)
  ret=file[grep(key, file)]
  pkey <- paste0("^", key, ".*\\[(.*)\\].*")
  ret <- gsub(pkey, "\\1", ret)
}

### conglomerate info from txt files from dcmdump
getInfo <- function(txt){
  xx <- readLines(txt)
  SeriesDesc <- extract(xx, "(0008,103e)")
  StudyDesc <- extract(xx, "(0008,1030)")
  itype <- extract(xx, "(0008,0008)")
  modal <- extract(xx, "(0008,0060)")
  SeriesNum <- extract(xx, "(0020,0011)")
  StudyID <- extract(xx, "(0020,0010)")
  GantryDetectorTilt <- extract(xx, "(0018,1120)")
  RotationDirection <- extract(xx, "(0018,1140)")
  TableHeight <- extract(xx, "(0018,1130)")
  ConvolutionKernel <- extract(xx, "(0018,1210)")
  
  if (length(SeriesDesc) == 0) SeriesDesc <- NA
  if (length(StudyDesc) == 0) StudyDesc <- NA
  if (length(StudyID) == 0) StudyID <- NA
  if (length(itype) == 0) itype <- NA
  if (length(modal) == 0) modal <- NA
  if (length(SeriesNum) == 0) SeriesNum <- NA
  if (length(GantryDetectorTilt) == 0) GantryDetectorTilt <- NA
  if (length(RotationDirection) == 0) RotationDirection <- NA
  if (length(TableHeight) == 0) TableHeight <- NA
  if (length(ConvolutionKernel) == 0) ConvolutionKernel <- NA
  
  return(c(SeriesDesc=SeriesDesc, StudyDesc=StudyDesc, itype=itype, Modality=modal,
           SeriesNum= SeriesNum, GantryDetectorTilt=GantryDetectorTilt,
           RotationDirection= RotationDirection, TableHeight= TableHeight, 
           ConvolutionKernel= ConvolutionKernel, StudyID=StudyID))
}

## get the basename for a file, for example getBase("blah.nii.gz", 2) = "blah"
  getBase <- function(x, ind=1){
    sapply(strsplit(x, split="\\."), function(xx) 
      paste(xx[1:(length(xx)-ind)], collapse=".", sep=""))
  }



includeMatrix <- function(basedir, keepAll = FALSE, keepMR = TRUE, 
  dropstring = NULL, verbose=TRUE, error=TRUE){

  sortdir <- file.path(basedir, "Sorted")

  #### checking if all got converted and dropping unneded scans
  tars <- basename(dir(path=sortdir, pattern=".tar.gz", full.names=TRUE, recursive=FALSE))
  niis <- basename(dir(path=basedir, pattern=".nii.gz", full.names=TRUE, recursive=FALSE))


  tars <- getBase(tars, ind = 2)
  niis <- getBase(niis, ind=2)
  if (error) {
    checker <- (niis %in% tars)
    if (!all(checker)) print(niis[!checker])
    stopifnot(all(checker))
  }
  mis <- tars[!(tars %in% niis)]
  if (!is.null(dropstring)) mis <- mis[!grepl(dropstring, mis)]
  if (verbose) print(mis)


  #############################################
  ######  Get information from Scans ##########
  ###### Drop those that aren't needed ########
  #############################################

  stubs <- basename(tars)
  stubs <- file.path(basedir, "Sorted", stubs)
  txts <- paste0(stubs, ".txt")


  stubs <- gsub("(.*)\\.txt", "\\1", txts)
  outs <- sapply(txts, getInfo)
  outs <- t(outs)
  rownames(outs) <- NULL
  outs <- data.frame(outs, stringsAsFactors=FALSE)
  outs$fname <- txts
  outs$SeriesNum <- as.numeric(outs$SeriesNum)
  outs$GantryDetectorTilt <- as.numeric(outs$GantryDetectorTilt)
  outs$StudyID <- as.numeric(outs$StudyID)
  
  ## take out MRIs, localizers, dose reports, bone, cervical
  outs$Takeout <- FALSE
  if (!keepMR) outs$Takeout <- outs$Takeout | outs$Modality %in% "MR"
  outs$Takeout <- outs$Takeout | grepl("LOCALIZER", outs$itype)
  outs$Takeout <- outs$Takeout | grepl("SECONDARY", outs$itype)
  
  outs$Takeout <- outs$Takeout | grepl("Dose Report", outs$SeriesDesc)
  ## Bone window scans - or diff convolution filters
  outs$Takeout <- outs$Takeout | grepl("BONE", outs$SeriesDesc)
  outs$Takeout <- outs$Takeout | grepl("Bone", outs$SeriesDesc)
  outs$Takeout <- outs$Takeout | grepl("H60s", outs$ConvolutionKernel)

  outs$Takeout <- outs$Takeout | grepl("CERVICAL", outs$StudyDesc)
  outs$Takeout <- outs$Takeout | grepl("Patient Protocol", outs$SeriesDesc)
  outs$Takeout <- outs$Takeout | grepl("SCOUT", outs$SeriesDesc)
  outs$Takeout <- outs$Takeout | grepl("CENTERING", outs$SeriesDesc)
  outs$Takeout <- outs$Takeout | grepl("SINUS", outs$SeriesDesc)
  if (!is.null(dropstring)) 
    outs$Takeout <- outs$Takeout | grepl(dropstring, outs$SeriesDesc)
  
  outs$Takeout <- outs$Takeout | outs$SeriesNum > 20
  
  outs$Takeout[is.na(outs$Takeout)] <- FALSE
  if (keepAll) outs$Takeout <- TRUE
  kept <- outs[! outs$Takeout, ]
  
  return(list(outs=outs, mis=mis))
}

file_gc <- function(basedir, progdir, verbose=TRUE, ...){
  
  ###
  gf <- getfiles(basedir)
  files <- gf$files
  paths <- gf$paths
  
  ## THEN PIPELINE TEST
  ipath <- 1
  lpath <- length(paths)
  
  ## THEN PIPELINE TEST
  ipath <- 1
  lpath <- length(paths)
  
  ### read headers and then the gantry tilt correction
  if (lpath >= 1){
    if (verbose) cat("Gantry Tilt Correction \n")
    
    for (ipath in 1:length(paths)){
      path <- paths[ipath]
      files <- getfiles(path)$files
      hdr <- dicomInfo(fname=files[1], pixelData=FALSE)$hdr
      row <- grep("GantryDetectorTilt", hdr[, "name"])
      if (length(row) == 0) next;
      if (length(row) > 1) {
        print("Too many gantry tilts!")
        next
      }
      if (length(row) == 1){
        tilt <- as.numeric(hdr[row, "value"])
        if (tilt != 0) gantry_correct(indir=path, 
          progdir=progdir, verbose=verbose)
      }
    }
  }
} # end function


### copy dicom header info
header_txt <- function(path){
  stub <- basename(path)
  dcm = getfiles(path)$files[1]
  make_txt(dcm, ispath=FALSE)
}



### dump the header contents and 
make_txt = function(path, ispath=TRUE,
  readfile = FALSE, 
  rmfile = FALSE) {
  stopifnot(!readfile | rmfile)
  dcm = path
  if (ispath) {
    dcm = getfiles(path)$files
  }
    stub <- basename(dcm)
    stub = gsub("\\.(dcm|DCM)$", "", stub)
    txtfile <- file.path(dirname(dcm), paste0(stub, ".txt"))
  df = data.frame(dcm = dcm, txtfile=txtfile, 
    stringsAsFactors = FALSE)
  dcm = txtfile = NULL

  hdrs = mlply(df, function(dcm, txtfile){
    if (rmfile) {
      hdr = system(sprintf('dcmdump "%s"', dcm), intern=TRUE)
    } else {
      hdr = system(sprintf('dcmdump "%s" > "%s"', dcm, txtfile), 
        intern=TRUE)
      if (readfile) hdr = readLines(txtfile)
    }
    return(hdr)   
  }, .progress= "none")

  return(hdrs)
}


### parse a make_txt list into a list of header information
parse_txt <- function(hdrs, verbose=TRUE){
  tables = llply(hdrs, function(hdr){
    hdr = str_trim(hdr)
    hdr = hdr[ grepl("^\\(", hdr)]
    ss = strsplit(hdr, " ")
    gel = sapply(ss, function(x) x[1])
    code = sapply(ss, function(x) x[c(2)])
    vlen= sapply(ss, function(x) gsub(",", "", x[length(x) -2]))
    name = sapply(ss, function(x) x[length(x)])
    value = sapply(ss, function(x) x[c(3)])

    gel = gsub("\\(|\\)", "", gel)
    gel = toupper(gel)
    ss.gel = strsplit(gel, ",")
    group = sapply(ss.gel, function(x) x[1])
    el = sapply(ss.gel, function(x) x[2])

    df = data.frame(group=group, element=el, name=name,
      code=code, length=vlen,
      value=value,  stringsAsFactors=FALSE)
    df$value = gsub("^\\[", "", df$value)
    df$value = gsub("\\]$", "", df$value)

    df$value[grep("\\.\\.\\.$", df$value)] = "skipped"
    df$value = gsub("\\", " ", df$value, fixed=TRUE)

    df = df[ !grepl("[A-Z]", df$group) | 
      (df$group %in% "7FE0" & df$element %in% "0010"), ]
    df[(df$group %in% "7FE0" & df$element %in% "0010"), 
      c("length", "value")] = c(-1, "skipped")

    df$value = gsub("(no", "", df$value, fixed=TRUE)
    return(df)
  }, .progress=ifelse(verbose, "text", "none"))
  return(tables)
}

rereadDICOMHeader = function(dcm, pixelData=FALSE,...){
  hdr.list = make_txt(dcm, ispath=FALSE, rmfile= TRUE, readfile=TRUE)
  hdrs = parse_txt(hdr.list, ...)
  return(hdrs)
}



### gantry tilt correction - calls in matlab
### gantry2_edit script
gantry_correct <- function(indir, progdir, verbose=TRUE){
  if (verbose) print(paste0("Gantry correction ", indir))
  find.matlab <- system("which matlab")
  cmd <- 'matlab -nodesktop -nosplash -nodisplay -r '  
  if (find.matlab != 0){
    cmd <- paste0('/Applications/MATLAB_R2012b.app/bin/', cmd)
  }
  ### gantry tilt correction - make new folder
  ### ranem old folder - zip it and then run matlab script
  indir <- path.expand(indir)
  dname <- dirname(indir)
  folname <- basename(indir)
  outname <- paste0(folname, "_ungantry")
  outdir <- file.path(dname, outname)
  if (!file.exists(outdir)) dir.create(outdir)
  system(sprintf('rm "%s"/*', outdir))
  system(sprintf('cp "%s"/*.dcm "%s"', indir, outdir))
  cmd <- paste(cmd, '"try, ')
  cmd <- paste(cmd, sprintf("addpath('%s');", file.path(progdir, "gantry")))
  cmd <- paste(cmd, sprintf("DIRlist(1,1).path_in = '%s';", indir))
  cmd <- paste(cmd, sprintf("DIRlist(1,1).path_out = '%s';", indir))
  cmd <- paste(cmd, "DIRlist(1,1).status = 'PROGRESS 0';")
  cmd <- paste(cmd, sprintf("DIRlist(1,1).skip = false;"))
  cmd <- paste(cmd, sprintf("gantry2_edit(DIRlist);"))
  cmd <- paste0(cmd, 'end; quit"')
  x <- system(cmd, ignore.stdout = !verbose )
  lastchar <- substr(outdir, nchar(outdir), nchar(outdir))
  while (lastchar == "/"){
    outdir <- substr(outdir, 1, nchar(outdir)-1)
    lastchar <- substr(outdir, nchar(outdir), nchar(outdir))
  }
  stopifnot( x== 0)
  

  
  # header_txt(outdir)
  gd <- getwd()
  setwd(dirname(outdir))
  bn <- basename(outdir)

  dcms = getfiles(outdir)$files
  hdrl = rereadDICOMHeader(dcms)      
  names(hdrl) = dcms

  dcmtables = dicomTable(hdrl)

  save(dcmtables, 
    file=paste0(outdir, "_Header_Info.Rda"))  

  system(sprintf('tar -czf "%s" ./"%s" --remove-files', paste0(outdir, ".tar.gz"), bn))  
  system(sprintf('rmdir "%s"', outdir))  
  
  setwd(gd)
  # header_txt(indir)
  
  return(outdir)
}




dcm2nii_worker <- function(path, outfile="output", 
  pixelSpacing = "0.45 0.45", addPixel = FALSE, ...){
  
  # the workhorse of the Rdcm2nii
  res = tryCatch({ 
    dcm = readDICOM(path=path, recursive=FALSE, 
    verbose=verbose)
    TRUE
  }, warning = function(war) {
    print(paste("Warning was: ", war))
    return(FALSE)
  }, error = function(e) {
    print(paste("Had error, skipping: ", e))
    return(FALSE)
  }, finally = function(e){
    return(TRUE)
  })
  ### added line for readDICOM

  if (!res) return(FALSE)

  dcmtable = dicomTable(dcm$hdr)
  keepcols = grepl("RescaleIntercept|RescaleSlope|PixelSpacing", 
                   colnames(dcmtable))
  dcmtab = dcmtable[, 
    c("0028-1052-RescaleIntercept", 
      "0028-1053-RescaleSlope", 
      "0028-0030-PixelSpacing"),
    drop=FALSE]
  stopifnot(ncol(dcmtab) == 3)
  colnames(dcmtab) = c("intercept", "slope", "pixelspacing")
  # print(dcmtab$pixelspacing)

  dcmtab$pixelspacing[ dcmtab$pixelspacing %in% "" ] = NA
  if (all(is.na(dcmtab$pixelspacing))){
    ## Adding in Pixel Spacing
    if (addPixel){
      if (verbose) cat("Changing PixelSpacing to 0.45")
      for (iimg in seq_along(dcm$hdr)){
        hd = dcm$hdr[[iimg]]
        hd[ hd$name == "PixelSpacing", "value"] = pixelSpacing
        dcm$hdr[[img]] = hd
      }
    } else { 
      return(FALSE)
    }
  }

  for (iimg in 1:length(dcm$img)){
    inter = as.numeric(dcmtab$intercept[iimg])
    slope = as.numeric(dcmtab$slope[iimg])
    dcm$img[[iimg]] = (dcm$img[[iimg]] + inter)/slope
    x = dcm$img[[iimg]]
    # if (verbose) print(range(dcm$img[[iimg]]))
    dcm$img[[iimg]][x < -1024] = -1024
    dcm$img[[iimg]][x > 3000] = 3000
  }

  # res = try({ 
  #   dcmNifti <- dicom2nifti(dcm, rescale=FALSE, reslice=FALSE, 
  #     descrip = NULL)
  # })

  res = tryCatch({ 
    dcmNifti <- dicom2nifti(dcm, rescale=FALSE, reslice=FALSE, 
      descrip = NULL)
    TRUE
  }, warning = function(war) {
    print(paste("Warning was: ", war))
    return(FALSE)
  }, error = function(e) {
    print(paste("Had error, skipping: ", e))
    return(FALSE)
  }, finally = function(e){
    return(TRUE)
  })
  ### added line for readDICOM

  if (inherits(res, "try-error")) return(FALSE)

  dcmNifti@scl_slope = 1
  dcmNifti@scl_inter = 0        

  #       vals <- t(sapply(dcm$hdr, function(x){
  #         int <- x[x[, "name"] == "RescaleIntercept", "value"]
  #         slope <- x[x[, "name"] == "RescaleSlope", "value"]
  #         return(c(int=int, slope=slope))
  #       }))
  #       class(vals) <- "numeric"
  #       vals <- data.frame(vals)
  #       uniq <- sapply(vals, function(x) length(unique(x)))
  #       stopifnot(all(uniq==1))
  #       vals <- apply(vals, 2, unique)
  #       dcmNifti@scl_slope <- vals["slope"]
  #       dcmNifti@scl_inter <- vals["int"]
  descrip.string <- extractHeader(dcm$hdr, 
    "SeriesDescription", 
    FALSE)[1]
  descrip.string[ is.na(descrip.string) ] = ""
  
  dcmNifti@descrip <- paste0("written by R - ", dcmNifti@descrip)
  outfile <- gsub("\\.gz$", "", outfile)
  outfile <- gsub("\\.nii$", "", outfile)
  writeNIfTI(nim=dcmNifti, file=file.path(path, outfile))
  return(TRUE)

}


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
    
    for (ipath in 1:lpath){
      path <- paths[ipath]
      res <- system(sprintf('rm "%s"/*.nii.gz', path), ignore.stdout = !verbose)
      
      res = dcm2nii_worker(path)
      if (!res) next;
      
      niis <- dir(path=path, pattern=".nii.gz")
      stub <- basename(path)
      
      iddir <- file.path(basedir)
      name <- stub
      
      ### copy dicom header info
      # header_txt(path)
      dcms = getfiles(path)$files

      if (length(dcms) > 0){
        # ifile = 1;
        ### read in EVERY HEADER from this 
        if (verbose) cat("Reading Headers \n")
        hdrl = rereadDICOMHeader(dcms)      
        names(hdrl) = dcms
    
        dcmtables = dicomTable(hdrl)

          save(dcmtables, 
            file=paste0(path, "_Header_Info.Rda"))
      }


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




### nnee to write this function - failing at differnent spots
readCT = function (fname, hdr=NULL, endian = "little", 
  flipud = TRUE, DICM = TRUE, 
    skipSequence = FALSE, pixelData = TRUE, warn = -1, debug = FALSE) 
{
    oldwarn <- getOption("warn")
    options(warn = warn)
    fsize <- file.info(fname)$size
    fraw <- readBin(fname, "raw", n = as.integer(fsize), endian = endian)
    skip128 <- fraw[1:128]
    if (debug) {
        cat("#", "First 128 bytes of DICOM header =", fill = TRUE)
        print(skip128)
    }
    if (DICM) {
        if (rawToChar(fraw[129:132]) != "DICM") {
            stop("DICM != DICM")
        }
    }
    dicomHeader <- sequence <- NULL
    seq.txt <- ""
    if (is.null(hdr)) {
      dcm <- parseDICOMHeader(fraw[133:fsize], seq.txt, endian = endian, 
         verbose = debug)
      hdr <- as.data.frame(dcm$header, stringsAsFactors = FALSE)
      row.names(hdr) <- NULL
      names(hdr) <- c("group", "element", "name", "code", "length", 
          "value", "sequence")
    } else {
        #  don't know yet
    }


    if (dcm$pixel.data && pixelData) {
        if (debug) {
            cat("##### Reading PixelData (7FE0,0010) #####", 
                fill = TRUE)
        }
        img <- parsePixelData(fraw[(132 + dcm$data.seek + 1):fsize], 
            hdr, endian, flipud)
    }
    else {
        if (dcm$spectroscopy.data && pixelData) {
            if (debug) {
                cat("##### Reading SpectroscopyData (5600,0020) #####", 
                  fill = TRUE)
            }
            img <- parseSpectroscopyData(fraw[(132 + dcm$data.seek + 
                1):fsize], hdr, endian)
        }
        else {
            img <- NULL
        }
    }
    options(warn = oldwarn)
    list(hdr = hdr, img = img)
}
