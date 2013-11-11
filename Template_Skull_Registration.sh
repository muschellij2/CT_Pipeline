################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Do Skull stripping on the Test 5 cases,
# Output: Skull stripped images
# Use of output: intensity -based normalization of CT
################################
# rootdir="/Volumes/DATA_LOCAL/Image_Processing"
rootdir="/dexter/disk2/smart/stroke_ct/ident"
tempdir="${rootdir}/Template"
basedir="${rootdir}/MISTIE"
progdir="${rootdir}/programs"
tempfile="scct_unsmooth_skull_400_1000.nii.gz"
tempimg="scct_unsmooth.nii.gz"
tempimg="${tempdir}/${tempimg}"
template="${tempdir}/${tempfile}"

cd "$basedir"
id="100-318"
for id in `find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \;`; 
 do
	refdir="$basedir/$id"
	# refdir="$basedir"
	cd "$refdir"
	echo "$refdir"
	if [ "$id" != "programs" ] && [ "$id" != "Skull_Stripped" ] \
		&& [ "$id" != "All_Images" ] && [ "$id" != "results" ]  
	then

		OUTDIR="$refdir/Registered" 
		mkdir "$OUTDIR"

		for file in *.nii.gz; 
		do
		case "$file" in
			*ROI* ) continue;;
			*Zeroed* ) continue;;
			*_MR_* ) continue;;
		esac;

	    stub=`basename $file`
	    stub=`echo $stub | awk '{ sub(/\.nii\.gz/, ""); print }'`
	    raw="${stub}_Skull_400_1000"
	    outfile="${OUTDIR}/${stub}_Skull_Registered"
	    outimg="${OUTDIR}/${stub}_Registered"

	    fslmaths $file -thr 400 -uthr 1000 "${OUTDIR}/${raw}"
	    flirt -in "${OUTDIR}/${raw}" -ref "${template}" -omat "$outimg.txt" -o "$outfile"
	#    rm "${outfile}.nii.gz"


	    flirt -verbose 1 -in "$refdir/${file}" -ref "${template}" -applyxfm -init "$outimg.txt" -o "$outimg"
	    echo "$outfile"
	    rm "${OUTDIR}/${raw}.nii.gz"
	    rm "${outfile}.nii.gz"

		done;
	fi;


done;


