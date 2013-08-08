# ids <- c("100-4001", "100-4003", "100-4034", "100-4054", "100-4066", "100-4087")
ids <- c("205-509", "205-517", "205-519", "225-502", "225-503", "225-504", 
         "225-505", "225-506", "225-507", "225-510", "225-511","225-515", "225-522", 
         "225-523", "225-524", "232-512", "232-513", "232-514", "232-516", "232-521", 
         "232-526", "289-518", "289-525")
# , "301-520" is weird

# ids <- "225-510"
# ids <- "100-4087"
rsync <- FALSE
func <- ifelse(rsync, "rsync", "scp")
iid <- ids[1]
for (iid in ids){
  login <- "jmuschel@enigma2.jhsph.edu"
  ss <- as.numeric(strsplit(iid, "-")[[1]][2])
  if (ss > 4000){
    study <- "CLEAR_III"
    dpath <- file.path("CLEAR", "CLEAR III")
  } else if (ss > 300 & ss < 500){
    dpath <- study <- "MISTIE"
  } else if (ss > 500 & ss < 4000) {
    dpath <- study <- "ICES" 
  }
  dexdir <- file.path("/dexter/disk2/smart/stroke_ct/ident", study)
  newid <- iid
  iddir <- file.path(dexdir, newid)
  
  cmd <- sprintf("ssh %s  mkdir -p %s\n", login, iddir)
  cat(cmd)
#   system(cmd)
  datadir <- file.path("/Volumes/DATA/Image Archive", dpath)
  diddir <- file.path(datadir, newid, "Anonymized")
  if (!file.exists(diddir)) {
    checkdir <- file.path(datadir, newid)
    dirs <- basename(list.dirs(path=checkdir, 
                      recursive=FALSE, full.names=FALSE))
    dirs <- dirs[!grepl("Source", dirs)]
    stopifnot(length(dirs)== 1)
    stopifnot(nchar(dirs) == 3)
    diddir <- file.path(checkdir, dirs[1])
  }
  cat(cmd <- sprintf('%s -rv "%s"/* %s:"%s"\n', func, diddir, login, iddir))
#   system(cmd)

}