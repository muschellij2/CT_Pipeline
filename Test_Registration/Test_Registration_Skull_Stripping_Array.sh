################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Do Skull stripping on the Test Registration cases,
# Output: Skull stripped images
# Use of output: intensity -based normalization of CT
################################
# rootdir="/Volumes/DATA_LOCAL/Image_Processing"
rootdir="/dexter/disk2/smart/stroke_ct/ident"
basedir="${rootdir}/Test_Registration"
progdir="${rootdir}/programs"
refdir="${basedir}/RawNIfTI"

outdir="${refdir}/Skull_Stripped"
mkdir -p "$outdir"

shopt -s extglob

cd "$refdir"
array=();
i=0;
for f in *.nii; 
do
case "$f" in
	*ROI* ) continue;;
	*Zeroed* ) continue;;
	bws* ) continue;;
	w* ) continue;;
	2mm* ) continue;;
	c* ) continue;;
	esac;
	array[ $i ]="$f";
	(( i++ ));
	echo $i;
done;
echo ${array[0]}
echo ${array[@]}

ind=`expr $SGE_TASK_ID - 1`
if [ -z "${ind}" ]; then
	ind=0;
fi
echo "ind is $ind"

f=${array[$ind]}
echo "f is $f"
FSLOUTPUTTYPE="NIFTI"
sh "${progdir}/Brain_Seg_Function.sh" -i "$f" -o \
    "$outdir" -f 0.1 -g -b
FSLOUTPUTTYPE="NIFTI_GZ"