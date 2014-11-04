#!/bin/bash 
#$ -t 1-112
### Maxvmem = 22G
R --no-save < Create_Predictors.R
