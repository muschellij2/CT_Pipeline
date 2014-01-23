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
basedir="${rootdir}/Test_Registration/"
imgdir="${basedir}/RawNIfTI"
progdir="${rootdir}/programs"
tempfile="scct_unsmooth_skull_400_1000.nii.gz"
tempimg="scct_unsmooth.nii.gz"
tempimg="${tempdir}/${tempimg}"
template="${tempdir}/${tempfile}"

cd "$imgdir"
OUTDIR="$basedir/FLIRT" 
# for id in `find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \;`; 
 # do
    # refdir="$basedir"
    cd "$imgdir"
    echo "$imgdir"
    if [ "$id" != "programs" ] && [ "$id" != "Skull_Stripped" ] \
        && [ "$id" != "All_Images" ] && [ "$id" != "results" ]  
    then

        file="102360_20100116_0042.nii"
        for file in *.nii; 
        do
        case "$file" in
            *ROI* ) continue;;
            *Zeroed* ) continue;;
            *_MR_* ) continue;;
        esac;

        lower=400
        upper=1000

        stub=`basename $file`
        roifile=`echo $stub | awk '{ sub(/_/, "_ROI_"); print }'`
        stub=`echo $stub | awk '{ sub(/\.nii/, ""); print }'`
        stub=`echo $stub | awk '{ sub(/\.gz/, ""); print }'`
        stubroi=`echo $stub | awk '{ sub(/_/, "_ROI_"); print }'`
        raw="${stub}_Skull_${lower}_${upper}"
        outimg="${OUTDIR}/affine12_${stub}"
        routimg="${OUTDIR}/affine12_${stubroi}"
        outimg2="${OUTDIR}/2mm_affine12_${stub}"
        routimg2="${OUTDIR}/2mm_affine12_${stubroi}"

        fslmaths $file -thr $lower -uthr $upper "${OUTDIR}/${raw}"
        # fslmaths "${OUTDIR}/${raw}" -bin "${OUTDIR}/${raw}"
        flirt -v -in "${OUTDIR}/${raw}" -ref "${template}" \
         -omat "$outimg.txt" -o "$outimg"
    #    rm "${outfile}.nii.gz"


        flirt -verbose 1 -in "$imgdir/${file}" -ref "${template}" \
        -applyxfm -init "$outimg.txt" -o "$outimg"
        ## subsampling to 2mm
        fslmaths "$outimg" -subsamp2 "$outimg2"
        echo "$outimg"


        flirt -verbose 1 -in "$imgdir/${roifile}" -ref "${template}" \
        -applyxfm -init "$outimg.txt" -o "$routimg"
        ## subsampling to 2mm
        fslmaths "$routimg" -subsamp2 "$routimg2"
        echo "$routimg"


### 9 DOF registration
        outimg="${OUTDIR}/affine9_${stub}"
        routimg="${OUTDIR}/affine9_${stubroi}"
        outimg2="${OUTDIR}/2mm_affine9_${stub}"
        routimg2="${OUTDIR}/2mm_affine9_${stubroi}"

        flirt -v -in "${OUTDIR}/${raw}" -ref "${template}" \
         -omat "$outimg.txt" -o "$outimg" -dof 9
    #    rm "${outfile}.nii.gz"


        flirt -verbose 1 -in "$imgdir/${file}" -ref "${template}" \
        -applyxfm -init "$outimg.txt" -o "$outimg"
        fslmaths "$outimg" -subsamp2 "$outimg2"        
        echo "$outimg"

        flirt -verbose 1 -in "$imgdir/${roifile}" -ref "${template}" \
        -applyxfm -init "$outimg.txt" -o "$routimg"
        fslmaths "$routimg" -subsamp2 "$routimg2"
        echo "$routimg"

        rm "${OUTDIR}/${raw}.nii.gz"        

        done;
    fi;


# done;


