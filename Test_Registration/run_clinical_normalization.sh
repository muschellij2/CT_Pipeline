#!/bin/bash 
matlab -nodesktop -nosplash -nodisplay -nojvm < run_clinical_normalization.m


getdim(){

    dim[0]=`fslhd $1 | grep ^dim1 | awk '{print $2}'`
    dim[1]=`fslhd $1 | grep ^dim2 | awk '{print $2}'`
    dim[2]=`fslhd $1 | grep ^dim3 | awk '{print $2}'`

    echo "${dim[@]:0:3}"
}

getpixdim(){

    dim[0]=`fslhd $1 | grep ^pixdim1 | awk '{ rounded = sprintf("%.0f", $2); print rounded }'`
    dim[1]=`fslhd $1 | grep ^pixdim2 | awk '{ rounded = sprintf("%.0f", $2); print rounded }'`
    dim[2]=`fslhd $1 | grep ^pixdim3 | awk '{ rounded = sprintf("%.0f", $2); print rounded }'`
    
    echo "${dim[@]:0:3}"
    
}

checkpixdim(){
    dim=`getpixdim $1`
    for i in ${dim[@]};
    do
        # echo "$i"
    # echo "$i"
       if [ "$i" -ne "$2" ]
       then
           return 1;
       fi;
    done;
    return 0;
}


## need to set this for for loops with null returns for strings
shopt -s nullglob


#### subsample all the images, appending 2mm on them
rootdir="/dexter/disk2/smart/stroke_ct/ident"
basedir="${rootdir}/Test_Registration/"
imgdir="${basedir}/reoriented"
progdir="${rootdir}/programs"

cd "$imgdir"
# for id in `find . -maxdepth 1 -mindepth 1 -type d -exec basename '{}' \;`; 
 # do
    # refdir="$basedir"

### delete previous files
rm 2mm_*.nii.gz

### delete voi files 
rm *.voi    

for file in *.nii; 
do
case "$file" in
    w* | bws* ) 

    stub=`basename $file`
    outfile=`echo $stub | awk '{ sub(/^bws|^w/, "2mm_"); print }'`
    ### only subsample if it's not 2mm already
    checkpixdim $file 2
    if [ "$?" -eq "1" ];
    then
        fslmaths $file -subsamp2 $outfile
    fi;
    echo "$outfile"
esac;

done;



for file in c*_sn.mat; 
do
    stub=`basename $file`
    outfile=`echo $stub | awk '{ sub(/^c/, ""); print }'`
    mv $stub $outfile
done;


# done;


