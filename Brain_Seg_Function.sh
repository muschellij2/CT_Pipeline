
#!/bin/bash

# Author: John Muschelli, 2013
# This code is to extract the brain using FSL
# BET (Brain Extraction Tool, Steve Smith) from CT Scans 
# Converted from dcm2nii (Chris Rorden)



function usage {
  echo "Usage: $0 -i File to be skull stripped"
  echo "          -o Output directory"
  echo "          -f Fraction used in Skull stripping"
  echo "          -h This page"
#  echo "          NEED TO ADD opts for which to do"
}

while getopts "hi:o:f:" flag
do
  case "$flag" in
    i)
      file=$OPTARG
      ;;
    o)
      OUTDIR=$OPTARG
      ;;
    f)
      intensity=$OPTARG
      ;;      
    h|?)
      usage
      exit 2
      ;;
  esac
done

if [ -z "${file}" ]; then
  echo "File is required"
  usage
  exit 2
fi

if [ -z "${OUTDIR}" ]; then
  echo "OUTDIR is required"
  usage
  exit 3
else
  mkdir -p ${OUTDIR}
fi


if [ -z "${intensity}" ]; then
  echo "No intensity given, using 0.35"
  intensity=0.35;
fi




  
      
      stub=`basename $file`
      zeroed=`echo $file | awk '{ sub(/\.nii\.gz/, "_Zeroed_"'${intensity}'"\.nii\.gz"); print }'`

      ### need this because then you can reorient them if you need to
      sform=`fslorient -getsformcode $file`  

      ### no adding or whatever
      raw=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS_No1024_"'${intensity}'"\.nii\.gz"); print }'`
      echo "No 1024 Bet $file file..";
      fslmaths $file -thr 0 -uthr 100 "${OUTDIR}/${raw}"
   
      echo "Bet 1 Running $raw"
      bet2 "$OUTDIR/$raw" "$OUTDIR/$raw" -f ${intensity}
      
      rawmask=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS_No1024_Mask_"'${intensity}'"\.nii\.gz"); print }'`
      fslmaths "${OUTDIR}/${raw}" -bin "${OUTDIR}/${rawmask}"
      fslmaths "$OUTDIR/${rawmask}" -fillh "$OUTDIR/${rawmask}"            

      minmax=`fslstats $file -R`
      min=`echo "$minmax" | cut -d ' ' -f 1`
      max=`echo "$minmax" | cut -d ' ' -f 2`


      result=`echo "($min < 0)" | bc`
      
      if [ $result ] 
      then
        echo "Rescaling so histograms are 'zeroed'"
        ### need to add this so that rest works - can't have pre-subtracted data.  Could adapt this to be anything other than 1024
        fslmaths $file -add 1024 $zeroed
        file="$zeroed" 
      fi 

# This is important because of 2 things - constraining to brain and getting rid of FOV
      h=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS_First_Pass_"'${intensity}'"\.nii\.gz"); print }'`
      echo "No Human extraction $file file..";
      fslmaths $file -thr 1024 -uthr 1124 "${OUTDIR}/${h}"

      echo "Brain from Room extraction $h file..";
      bet2 "${OUTDIR}/${h}" "${OUTDIR}/${h}" -f ${intensity}

      fslmaths "${OUTDIR}/${h}" -sub 1024 "${OUTDIR}/${h}"

      fpmask=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS_First_Pass_Mask_"'${intensity}'"\.nii\.gz"); print }'`
      fslmaths "${OUTDIR}/${h}" -bin "${OUTDIR}/${fpmask}"
      fslmaths "$OUTDIR/${fpmask}" -fillh "$OUTDIR/${fpmask}"

      ## fslmaths "$OUTDIR/$h" -thr 1024 -bin "$OUTDIR/$h"
      
      ### just trying bet straight up
      human=`echo $stub | awk '{ sub(/\.nii\.gz/, "_Human_"'${intensity}'"\.nii\.gz"); print }'`      
      echo "Human Extraction $human file..";
      bet2 $file "$OUTDIR/${human}" -f ${intensity}
      
    #   j=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS_"'${intensity}'"\.nii\.gz"); print }'`
    #   jmask=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS_Mask_"'${intensity}'"\.nii\.gz"); print }'`

    #   bbet=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS2_"'${intensity}'"\.nii\.gz"); print }'`
    #   bbetmask=`echo $stub | awk '{ sub(/\.nii\.gz/, "_SS2_Mask_"'${intensity}'"\.nii\.gz"); print }'`

    #   echo "Thresholding to range of 1024-1124 $human file..";
    #   fslmaths "$OUTDIR/${human}" -thr 1024 -uthr 1124 "$OUTDIR/${j}"
      
    #   echo "Bet 1 Running ${j}"
    #   bet2 "$OUTDIR/${j}" "$OUTDIR/${j}" -f ${intensity}
      
    #   echo "Bet 2 Running ${j}"
    #   ### this is if we want meshes
    # # bet $j $j -f $intensity -A
    #   bet2 "$OUTDIR/${j}" "$OUTDIR/${bbet}" -f ${intensity}



    # # echo "Translating to 0-100 range"
    #   echo "Making Binary Image"
    #   fslmaths "$OUTDIR/${j}" -thr 1024 -bin "$OUTDIR/${jmask}"
    #   # filling the holes
    #   fslmaths "$OUTDIR/${jmask}" -fillh "$OUTDIR/${jmask}"
    # # echo "Making Binary Image"
    # # fslmaths "$j" -thr 0 "$j"

    #   fslmaths "$OUTDIR/${bbet}" -thr 1024 -bin "$OUTDIR/${bbetmask}"
    #   fslmaths "$OUTDIR/${bbetmask}" -fillh "$OUTDIR/${bbetmask}"
    # # fslmaths "$bbet" -thr 0 "$bbet"

    # # echo "Translating to 0-100 range"
    #   echo "Subtracting 1024 from Image"
    #   fslmaths "$OUTDIR/${j}" -sub 1024 "$OUTDIR/${j}"
    # # echo "Making Binary Image"
    # # fslmaths "$j" -thr 0 "$j"

    #   fslmaths "$OUTDIR/${bbet}" -sub 1024 "$OUTDIR/${bbet}"


      fslmaths "$OUTDIR/${human}" -sub 1024 "$OUTDIR/${human}"
      # fslmaths "$OUTDIR/${human}" -thr 1024 "$OUTDIR/${human}"

      brainvol=`fslstats "$OUTDIR/${fpmask}" -V | awk '{ print $2/1000 }'`;
      echo "Brain volume is $brainvol"
      if [ $result ]
      then
        echo "Deleting ${zeroed} for cleanup"
        rm "${zeroed}"
      fi 

