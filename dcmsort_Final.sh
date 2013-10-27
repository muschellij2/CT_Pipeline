#!/bin/bash

# Author: Chuck Theobald, March 2012
# This script depends upon DCMTK's dcmdump program.
# It will copy a directory full of DICOM files into the
# given output directory with a sub-directory named for
# the ProtocolName group,element.  Alternative sortings
# may be had by judicious editing and modification of this
# script.

function usage {
  echo "Usage: $0 -D <DICOM directory>"
  echo "          -o Output directory"
  echo "          -m Move instead of copy"
  echo "          -x Exclude Localizers"
  echo "          -h This page"
  echo "          -s Series Date/Time instead of Series Date/Time"
}

cmd="cp";
exclude='';
while getopts "hD:o:mxs" flag
do
  case "$flag" in
    D)
      DICOMDIR=$OPTARG
      ;;
    o)
      OUTDIR=$OPTARG
      ;;
    m)
      cmd="mv"
      ;;
    x)
      exclude='true'
      ;;   
    s)
      study='true'
      ;;            
    h|?)
      usage
      exit 2
      ;;
  esac
done


# DICOMDIR=$(printf '%q' "$DICOMDIR")
# OUTDIR=$(printf '%q' "$OUTDIR")

if [ -z "${DICOMDIR}" ]; then
  echo "DICOMDIR is required"
  usage
  exit 2
fi

if [ -z "${OUTDIR}" ]; then
  echo "OUTDIR is required"
  usage
  exit 3
else
  mkdir -p "${OUTDIR}"
fi

# Select each file in given DICOM directory.
for d in "${DICOMDIR}"/*.dcm
do
  ## checking for directories
  if [[ -d $d ]]; then
      echo "$d is a directory"
  elif [[ -f $d ]]; then
      # echo "$PASSED is a file"

    echo "${d}"    

    dcmdump "${d}" > tmp.txt
    
    PID=$(cat tmp.txt | grep -e '^(0010,0020)' | sed -e 's/^(0010,0020).*\[\(.*\)\].*/\1/')
    StudyDate=$(cat tmp.txt | grep -e '^(0008,0020)' | sed -e 's/^(0008,0020).*\[\(.*\)\].*/\1/')
    StudyTime=$(cat tmp.txt | grep -e '^(0008,0030)' | sed -e 's/^(0008,0030).*\[\(.*\)\].*/\1/')
    NUMBER=$(echo $StudyTime | awk '{ print $0 / 100 }')
    NUMBER=$(echo $NUMBER | awk -F. '{ print $1 }')
    NUMBER=$(printf "%04.0f" $NUMBER)

    # ProtocolName=$(cat tmp.txt | grep -e '^(0018,1030)' | sed -e 's/^(0018,1030).*\[\(.*\)\].*/\1/')
    # PP=$(echo $ProtocolName | awk '{ gsub(" ","_"); print }')
    # PP=$(echo $PP | awk '{ gsub("/","_"); print }')
    # PP=$(echo $PP | awk '{ gsub("/","_"); print }')
    StudyDesc=$(cat tmp.txt | grep -e '^(0008,1030)' | sed -e 's/^(0008,1030).*\[\(.*\)\].*/\1/')
    StudyDesc=$(echo $StudyDesc | awk '{ gsub(" ","_"); print }')
    StudyDesc=$(echo $StudyDesc | awk '{ gsub("\\(0008,1030\\)_LO_\\(no_value_available\\)","unnamed"); print }')
    StudyDesc=$(echo $StudyDesc | awk '{ gsub("/",""); print }')
    StudyDesc=$(echo $StudyDesc | awk '{ gsub(/[:\47]*/,""); print }')
    

    SeriesDesc=$(cat tmp.txt | grep -e '^(0008,103e)' | sed -e 's/^(0008,103e).*\[\(.*\)\].*/\1/')
    SeriesDesc=$(echo $SeriesDesc | awk '{ gsub(" ","_"); print }')
    SeriesDesc=$(echo $SeriesDesc | awk '{ gsub("\\(0008,103e\\)_LO_\\(no_value_available\\)","unnamed"); print }')
    SeriesDesc=$(echo $SeriesDesc | awk '{ gsub("\\(no_value_available\\)","unnamed"); print }')
    SeriesDesc=$(echo $SeriesDesc | awk '{ gsub("_#_0,_0_SeriesDescription",""); print }')
    SeriesDesc=$(echo $SeriesDesc | awk '{ gsub("/",""); print }')
    SeriesDesc=$(echo $SeriesDesc | awk '{ gsub(/[:\47]*/,""); print }')

    SeriesDate=$(cat tmp.txt | grep -e '^(0008,0020)' | sed -e 's/^(0008,0020).*\[\(.*\)\].*/\1/')
    SeriesTime=$(cat tmp.txt | grep -e '^(0008,0031)' | sed -e 's/^(0008,0031).*\[\(.*\)\].*/\1/')
    SENUMBER=$(echo $SeriesTime | awk '{ print $0 / 100 }')
    SENUMBER=$(echo $SENUMBER | awk -F. '{ print $1 }')
    SENUMBER=$(printf "%04.0f" $SENUMBER)

    SID=$(cat tmp.txt | grep -e '^(0020,000d)' | sed -e 's/^(0020,000d).*\[\(.*\)\].*/\1/')
  #  echo $SeriesDesc
    Modality=$(cat tmp.txt | grep -e '^(0008,0060)' | sed -e 's/^(0008,0060).*\[\(.*\)\].*/\1/')
  # (0020,000d)
    UUID=$(cat tmp.txt | grep -e '^(0008,0018)' | sed -e 's/^(0008,0018).*\[\(.*\)\].*/\1/')
    # UUID=$(cat 3008.txt | grep -e '^(0008,0018)' | sed -e 's/^(0008,0018).*\[\(.*\)\].*/\1/')
    UUID=$(echo $UUID | sed -e 's/\(.*\)\.\(.*\)$/\1/')
    UUID=$(echo $UUID | sed -e 's/\(.*\)\.\(.*\)$/\1/')
    
    SNUM=$(cat tmp.txt | grep -e '^(0020,0011)' | sed -e 's/^(0020,0011).*\[\(.*\)\].*/\1/')
    
  # 
    # PNAME="${PID}_${StudyDate}_${NUMBER}_${PP}"
    # PNAME="${PID}_${StudyDate}_${NUMBER}_${StudyDesc}_${SeriesDesc}"
    # PNAME="${PID}_${StudyDate}_${NUMBER}_${Modality}_${SeriesDesc}"
    # PNAME="${PID}_${StudyDate}_${NUMBER}_${Modality}_${SID}_${SeriesDesc}"
    # PNAME="${PID}_${StudyDate}_${NUMBER}_${Modality}_${UUID}_${SeriesDesc}"
    # PNAME="${PID}_${StudyDate}_${NUMBER}_${Modality}_${SID}_${SNUM}_${SeriesDesc}"
    # PNAME="${PID}_${StudyDate}_${NUMBER}_${Modality}_${SID}_${SNUM}_${StudyDesc}_${SeriesDesc}"
    DATER="${SeriesDate}_${SENUMBER}"
    if [[ ! -z "${study}" ]]; then
        DATER="${StudyDate}_${NUMBER}"
    fi

    PNAME="${PID}_${DATER}_${Modality}_${SNUM}_${StudyDesc}_${SeriesDesc}"
    PNAME=$(echo $PNAME | awk '{ gsub(" ","_"); print }')

    if [[ ! -z "${exclude}" ]]; then
        # echo "is it excluded?"
        ITYPE=$(cat tmp.txt | grep -e '^(0008,0008)' | sed -e 's/^(0008,0008).*\[\(.*\)\].*/\1/')
        ITYPE=$(echo $ITYPE | tr '[:lower:]' '[:upper:]')
        change=$(echo $ITYPE | grep -o LOCALIZER)
        # echo "$change"
        # echo "$ITYPE"
        if [[ -n "${change}" ]]; then
          # echo "Localizer"
          PNAME="LOCALIZER_${PNAME}"
        fi
    fi


      # echo "Exclude is $exclude"
      # echo "$PNAME"
      # exit 1

    DESTDIR="${OUTDIR}/${PNAME}"
    # $(printf '%q' "$DICOMDIR")
    if [ ! -d "${DESTDIR}" ]; then
      mkdir "${DESTDIR}"
    fi
    $cmd "${d}" "${DESTDIR}"/
    rm tmp.txt
    
  fi    
done
# ls ${OUTDIR}

