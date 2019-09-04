##############################################
## This code is for Creating Predictors 
# for each individual
##
## Author: John Muschelli
## Last updated: May 20, 2014
##############################################
##############################################
rm(list =  ls())
library(methods)
library(neurobase)
library(dplyr)
library(readr)
library(tidyr)
RNGversion("3.5.0")
rootdir = file.path("/Volumes/DATA_LOCAL", 
                    "Image_Processing")
if (Sys.info()[["user"]] %in% "jmuschel") {
  rootdir = Sys.getenv("dex")
}
basedir = file.path(rootdir, 
                    "PITCH_reconverted", 
                    "processed_data")
resdir = file.path(rootdir, 
                   "PITCH_reconverted", 
                   "results")

ids = list.dirs(recursive = FALSE,
                path = basedir,
                full.names = TRUE)
imgs = sapply(ids, function(x) {
  f = list.files(pattern = "ROI.nii.gz", 
                 path = x, 
                 recursive = FALSE,
                 full.names = TRUE)
})
df = data_frame(
  id = basename(ids),
  id_dir = ids,
  roi = imgs,
  img = sub("ROI[.]", ".", roi),
  stub = nii.stub(img, bn = TRUE))
df = df %>% 
  mutate(
    outdir = file.path(id_dir, 
                       "processed"),
    outfile = file.path(outdir,
                        "predictors.rds"))


ids = c("100-318", "100-362", "100-365", "101-306",   "101-307",   "101-308",   "102-317",   "102-322",   
        "102-323",   "102-324",   "102-331",   "102-360",   "102-367",   "102-374",   "102-391",   "102-393",   
        "102-403",   "102-406",   "120-376",   "131-310",   "131-316",   "131-334",   "131-354",   "133-409",  
        "133-417",   "134-304",   "134-305",   "134-320",   "134-327",   "134-345",   "134-380",   "134-381",   
        "134-382",   "134-392",   "134-408",   "134-412",   "134-416",   "152-302",   "152-303",   "152-353",   
        "157-328",   "157-329",   "157-332",   "157-335",   "157-336",   "157-370",   "157-372",   "157-399",   
        "157-410",   "161-413",   "173-312",   "173-313",   "173-325",   "173-341",   "173-361",   "173-364",   
        "173-368",   "173-384",   "173-396",   "173-404",   "175-387",   "175-397",   "175-405",   "179-373",   
        "179-383",   "179-386",   "179-394",   "179-395",   "179-402",   "184-388",   "191-301",   "191-311",   
        "191-314",   "191-315",   "191-319",   "191-321",   "191-333",   "191-375",   "191-400",   "205-509",   
        "205-517",   "205-519",   "216-390",   "219-350",   "222-337",   "222-357",   "222-358",   "223-355",   
        "223-369",   "223-407",   "225-502",   "225-503",   "225-504",   "225-505",   "225-506",   "225-507",   
        "225-510",   "225-511",   "225-515",   "225-524",   "230-356",   "230-363",   "230-366",   "230-371",   
        "230-377",   "232-514",   "232-516",   "234-385",   "265-389",   "265-398",   "289-518",   "289-525")

df = df[ df$id %in% ids, ]




###########################################
# Create the groupings
###########################################
set.seed(20150504)
non.aggmods = 10
vec = seq(length(ids))
ind = sample(1:length(vec), size=non.aggmods)
train.ind = sort(vec[ind])
vec = vec[-ind]
nvalid = ceiling(length(vec)/2)
ind = sample(1:length(vec), size=nvalid)
valid.ind = sort(vec[ind])
test.ind = sort(vec[-ind])

ddf = data.frame(id = ids)
ddf$group = NA
ddf$group[train.ind] = "Train"
ddf$group[valid.ind] = "Validation"
ddf$group[test.ind] = "Test"

df = merge(df, ddf, by="id", sort=FALSE)
if (any(is.na(df$group))){
  stop(paste0("Something went wrong with ", 
              "group designation"))
}

df = df %>% 
  mutate(
    usemask = sub("[.]nii", "_usemask.nii",
                  img),
    usemask = file.path(outdir, 
                        basename(usemask)))


filename = file.path(resdir, 
                     "iso_filename.rds")
saveRDS(df, file = filename)

##########################################
# CNN grouping
##########################################
fdf = df %>% 
  dplyr::as_tibble() %>% 
  mutate(scan = sub("-", "", stub))  %>% 
  separate(scan, into = c("iid", "date", "time", "other"),
           extra = "merge", sep = "_") %>% 
  select(-other) %>% 
  mutate(no_dash_id = gsub("-", "", id))

split_file = file.path(resdir, "dataset_split.csv")
split_df = read_csv(split_file, col_names = FALSE)
colnames(split_df) = c("scan", "cnn_group")
split_df = split_df %>% 
  separate(scan, into = c("no_dash_id", "date", "time", "other"),
           extra = "merge", sep = "_", remove = FALSE) %>% 
  select(-other)

aj = anti_join(fdf, split_df)
stopifnot(nrow(aj) == 0)
fdf = left_join(fdf, split_df)
stopifnot(!any(is.na(fdf$cnn_group)))

filename = file.path(resdir, 
                     "iso_filename_with_cnn.rds")
saveRDS(fdf, file = filename)
