################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Do Skull stripping on the Test 5 cases,
# Output: Skull stripped images
# Use of output: intensity -based normalization of CT
################################
# rootdir="/Volumes/DATA_LOCAL/Image_Processing"
# rootdir="/dexter/disk2/smart/stroke_ct/ident"
rootdir="/Users/muschellij2/CT_Registration"
basedir="${rootdir}/SS_Validation"
progdir="${rootdir}/programs"

cd "$basedir"
# id="238-4160"
# for id in `find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \;`; 
 # do
	# refdir="$basedir/$id"
	refdir="$basedir"
	cd "$refdir"
	# if [ "$id" != "programs" ] && [ "$id" != "Sorted" ] && [ "$id" != "ROI_Images" ] && [ "$id" != "Original_Images" ]  && [ "$id" != "results" ]  
	# then 
	mkdir "$refdir/Skull_Stripped"

	###f="205-519_20110630_0633_3.nii.gz"
	#f="205-509_20100418_1312_2.nii.gz"
	# shopt -s extglob;
	## exclude ROIs
	for f in *.nii.gz; 
	do
	case "$f" in
 		*ROI* ) continue;;
 		*Zeroed* ) continue;;
		esac;
		# exit 1;
			# sh "${progdir}/Brain_Seg_Function.sh" -i "$f" -o "$refdir/Skull_Stripped" -f 0.35
			sh "${progdir}/Brain_Seg_Function.sh" -i "$f" -o "$refdir/Skull_Stripped" -f 0.1
			# sh "${progdir}/Brain_Seg_Function.sh" -i "$f" -o "$refdir/Skull_Stripped" -f 0.7
	done;
	# fi;
# done;


