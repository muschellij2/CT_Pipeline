###################
#### read in baseline NIHSS so that we can cross ref with Image
####
options(stringsAsFactors=FALSE)

rootdir = "~/CT_Registration"
rootdir = path.expand(rootdir)
progdir = file.path(rootdir, "programs")
datadir = file.path(rootdir, "data")
nihss = read.csv(file.path(datadir, "baseline_NIHSS.csv"), 
                 na.strings = c("NA", ".", ""))
load(file.path(datadir, "Registration_Image_Names.Rda"))

ids = df$id
ids = as.numeric(gsub("-", "", ids))

nihss = nihss[ nihss$patientName %in% ids, ]
