#!/bin/bash
#PBS -q hotel
#PBS -N 16S_dada_err2
#PBS -l nodes=1:ppn=8
#PBS -l walltime=40:00:00
#PBS -o /projects/ps-shurinlab/users/cbwall/YoZoop_microbiome_21.22/output/outs/Rscript_out
#PBS -e /projects/ps-shurinlab/users/cbwall/YoZoop_microbiome_21.22/output/errs/Rscript_err
#PBS -V
#PBS -M cbwall@ucsd.edu
#PBS -m ae
#PBS -A shurin-group

export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
module load R/4.1.2
Rscript --no-save /projects/ps-shurinlab/users/cbwall/scripts/16S_dada_err_2.R

