#!/bin/bash 
#$ -t 1-5
R --no-save < Aggregate_Model.R
### Maxvmem = 31G
