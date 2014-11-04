rm(list=ls())
library(fslr)
dexdir = "/dexter/disk2/smart/stroke_ct/ident"
regdir = file.path(dexdir, "Registration")
basedir = "~/CT_Registration/Final_Brain_Seg/results"
x = list.files(pattern=".png", 
	path = basedir)

df = data.frame(png = x, stringsAsFactors = FALSE)
df$id = gsub("(.*)_\\d{8}.*", "\\1", df$png)
df$iddir = file.path(regdir, df$id)
df$img = paste0(nii.stub(df$png), ".nii.gz")

write.csv(df[,"img", drop=FALSE], 
	file = file.path(basedir, "Image_Filenames.csv"), 
	row.names=FALSE)