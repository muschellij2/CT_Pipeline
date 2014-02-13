###################
#### read in baseline NIHSS so that we can cross ref with Image
####
options(stringsAsFactors=FALSE)

rootdir = "~/Dropbox/CTR/DHanley/MISTIE"
rootdir = path.expand(rootdir)
homedir = file.path(rootdir, "MISTIE DSMB Analysis")
datadir = file.path(homedir, "statadata")
nihss = read.csv(file.path(datadir, "baseline_NIHSS.csv"), 
                 na.strings = c("NA", ".", ""))