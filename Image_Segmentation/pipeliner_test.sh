#!/bin/bash 

rootdir="/dexter/disk2/smart/stroke_ct/ident"
basedir="${rootdir}/Test_5"
progdir="${rootdir}/programs"
nfolds=75;
#"\\d\\d\\d-(\\d|)\\d\\d\\d"
nid=`find $basedir -mindepth 1 -maxdepth 1 -type d -regex '.*[0-9][0-9][0-9]-[0-9]?[0-9][0-9][0-9]$' | wc -l`
# find $basedir -regex '{[0-9]}'

cd $progdir

### Make NIfTI Files
qsub -cwd -N test1 -l mem_free=2G,h_vmem=8G -r no \
test.sh 

### Brain extraction and mask generation
qsub -cwd -hold_jid test1 -N test2 \
test.sh

echo 0 > check.txt
while read line           
do           
    if [ "$line" -ne "0" ]
	then
	exit 1
	fi
done < check.txt 
echo 0 >> check.txt
