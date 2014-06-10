rm(list=ls())
library(xtable)
library(plyr)
# score = "NIHSS"
setwd('~/Dropbox/CTR/DHanley/CT_Registration/CT_Pipeline')

load("Top_0.01_Pvalues_df.Rda")
n = 10
nihss.tab = head(pvalimg.tab[["EVE_1"]], n)
nihss.tab$nvox = sprintf("%02.1f", nihss.tab$nvox)
names(nihss.tab) = c("NIHSS ROI", "Area")

load("GCS_Top_1000_Pvalues_df.Rda")
gcs.tab = head(pvalimg.tab[["EVE_1"]], n)
gcs.tab$nvox = sprintf("%02.1f", gcs.tab$nvox)
names(gcs.tab) = c("GCS ROI", "Area")

load("Population_Table_Engagement.Rda")
pop.tab = head(xtabs[["EVE_1"]], n)
pop.tab$EVE_1 = sprintf("%02.1f", pop.tab$EVE_1)

names(pop.tab) = c("Area", "Population Prevalence")

df = merge(pop.tab, nihss.tab, sort=FALSE, all = TRUE)
df = merge(df, gcs.tab, sort=FALSE, all = TRUE)

# names(jhut1.list)[names(jhut1.list) == "Background"] = "Ventricles"

proper = function(x){
  substring(x, 1, 1) <- toupper(substring(x, 1, 1))
  substring(x, 2) <- tolower(substring(x, 2))
  x
}

df$Area = proper(df$Area)
df$Area = revalue(df$Area, c("Background" = "CSF"))
df$Area = gsub("_", " ", df$Area)
df = df[order(as.numeric(df$"Population Prevalence"), 
              as.numeric(df[, "NIHSS ROI"]), 
              as.numeric(df[, "GCS ROI"]), 
              decreasing=TRUE),]
df[sapply(df, is.na)] = ""
rownames(df) = NULL

xtab = xtable(df)

xtab = xtable(df, 
              caption= 
                paste0("Distribution of the top 10 areas of engagement ",
                       "for population ",
                       " 3D histogram, the NIHSS ROI was based on a p-value ",
                       " threshold of 0.01, the GCS ROI was based on voxels ",
                       "with 1000 smallest p-values.  Each value represents ", 
                       "the percentage of the ROI engaged in this area.  ",
                       "The population-level areas are percentages ",
                       "weighted by proportion. Each distribution ",
                       " of areas is based on the Eve atlas.  ",
                       "We see that the population is engaged in ",
                       "areas of the CSF, such as the ventricles, and ",
                       " the insular and putaminal regions.",
                       " The ROI based on the NIHSS analysis engages ",
                       "primarily areas of the internal capsule ",
                       "and ventricular regions. ",                       
                       "The ROI based on the GCS analysis engages ",
                       "primarily the left thalamus and superior ",
                       "corona radiata.",
                       "The Eve", 
                       "atlas can be used to calculate area engagement",
                       "on a per-scan level as well."), 
              align=c("llccc"),
              label="t:breakdown")
print.xtable(xtab, file="breakdown.tex", include.rownames = FALSE)



