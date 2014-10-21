#!/bin/bash 
#$ -t 1-112
# maxvmem= 20G
R --no-save < Predict_From_Model.R
