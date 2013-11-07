#!/bin/bash 
exit_code=100
R --no-save < make_nifti.R
exit $exit_code