#!/bin/bash 
#$ -t 1-9
R --no-save < Aggregate_Model.R
### Maxvmem = 31G
