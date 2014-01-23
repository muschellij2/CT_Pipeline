#!/bin/bash 
R --no-save < test.R
# exstat=$?
if [ "$?" -ne "0" ]
then
	exit 100
fi
