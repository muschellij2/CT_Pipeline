# ids <- c("100-4001", "100-4003", "100-4034", "100-4054", "100-4066", "100-4087")
# ids <- c("205-509", "205-517", "205-519", "225-502", "225-503", "225-504", 
#          "225-505", "225-506", "225-507", "225-510", "225-511","225-515", "225-522", 
#          "225-523", "225-524", "232-512", "232-513", "232-514", "232-516", "232-521", 
#          "232-526", "289-518", "289-525", "301-520")

ids <- c("100-318", "100-362", "100-365", "101-306", "101-307", "101-308", "102-317", "102-322", "102-323", "102-324", "102-326", "102-331", 
  "102-347", "102-349", "102-351", "102-360", "102-367", "102-374", "102-391", "102-393", "102-403", "102-406", "111-415", "120-376", 
  "131-310", "131-316", "131-334", "131-354", "133-409", "133-417", "134-304", "134-305", "134-320", "134-327", "134-338", "134-343", 
  "134-345", "134-379", "134-380", "134-381", "134-382", "134-392", "134-408", "134-411", "134-412", "134-416", "152-302", "152-303", 
  "152-348", "152-353", "157-328", "157-329", "157-332", "157-335", "157-336", "157-370", "157-372", "157-399", "157-410", "161-413", 
  "173-312", "173-313", "173-325", "173-339", "173-341", "173-361", "173-364", "173-368", "173-384", "173-396", "173-404", "175-387", 
  "175-397", "175-405", "179-359", "179-373", "179-383", "179-386", "179-394", "179-395", "179-401", "179-402", "184-342", "184-388", 
  "191-301", "191-309", "191-311", "191-314", "191-315", "191-319", "191-321", "191-330", "191-333", "191-375", "191-378", "191-400", 
  "210-344", "216-390", "216-414", "219-340", "219-350", "222-337", "222-357", "222-358", "223-355", "223-369", "223-407", "230-346", 
  "230-352", "230-356", "230-363", "230-366", "230-371", "230-377", "234-385", "265-389", "265-398")

# ,  is weird

# ids <- "225-510"
# ids <- "100-4087"
rsync <- FALSE
makedir <- TRUE
func <- ifelse(rsync, "rsync", "scp")
iout <- 2
outdir <- c("/dexter/disk2/smart/stroke_ct/ident", "/Volumes/Seagate Backup Plus Drive/Image Archive")[iout]
cluster <- TRUE
if (grepl("Volumes", outdir)) {
  func <- "cp"
  cluster <- FALSE
}
iid <- ids[1]

for (iid in ids){
  ss <- as.numeric(strsplit(iid, "-")[[1]][2])
  if (ss > 4000){
    study <- "CLEAR_III"
    dpath <- file.path("CLEAR", "CLEAR III")
  } else if (ss > 300 & ss < 500){
    dpath <- study <- "MISTIE"
  } else if (ss > 500 & ss < 4000) {
    dpath <- study <- "ICES" 
  }
  dexdir <- file.path(outdir, study)
  newid <- iid
  iddir <- file.path(dexdir, newid)
  addcmd <- ""
  if (cluster) {
    login <- "jmuschel@enigma2.jhsph.edu"  
    addcmd <- sprintf("ssh %s", login)
  }
  cmd <- ""
  if (makedir) cmd <- sprintf('%s mkdir -p "%s";\n\n', addcmd, iddir)
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
  cat("\n\n")
#   system(cmd)

}