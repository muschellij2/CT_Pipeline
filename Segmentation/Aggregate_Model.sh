#!/bin/bash 
#$ -t 1,4
R --no-save < Aggregate_Model.R
### Maxvmem = 31G
