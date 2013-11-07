################################
# Written 2013Oct24
# Author: John Muschelli
# Purpose: Do Skull stripping on the Test 5 cases,
# Output: Skull stripped images
# Use of output: intensity -based normalization of CT
################################
# rootdir="/Volumes/DATA_LOCAL/Image_Processing"
rootdir="/dexter/disk2/smart/stroke_ct/ident"
basedir="${rootdir}/Test_5"
progdir="${rootdir}/programs"

cd "$basedir"
refdir="$basedir"
cd "$refdir"

mkdir "$refdir/Skull_Stripped"

i=0;
for f in *.nii.gz; 
do
case "$f" in
	*ROI* ) continue;;
	*Zeroed* ) continue;;
	esac;
	array[ $i ]="$f";
	(( i++ ));
	echo $i;
done;
echo ${array[0]}

ind=`expr $SGE_TASK_ID - 1`
f=${array[$ind]}			
sh "${progdir}/Brain_Seg_Function.sh" -i "$f" -o "$refdir/Skull_Stripped" -f 0.1


