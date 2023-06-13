#!/bin/bash

############################# BEGIN PBS DIRECTIVES #############################

### mail on job exit
#PBS -m e

### email address to send job updates
#PBS -M jeremy.schraw@bcm.edu

### keep job output and error files in your home directory
### this overrides -e and -o options to specify error/output file locations
#PBS -k oe

### request 4 CPUs on 1 node
#PBS -l nodes=1:ppn=4

### request 48gb virtual memory TOTAL
#PBS -l vmem=48gb

### request 1 day walltime
#PBS -l walltime=24:00:00

############################## END PBS DIRECTIVES ##############################

### enable 'module' commands
source /etc/profile.d/modules.sh

### JOB PREP ###

module load R/4.0.5

workingdir="/mount/genepi2/Old_genepi2/Jeremy/GOBACK/R_scripts/GOBACK-R-Scripts/"

## JOBS ##

## SINGLE VARIANT
#################

Rscript $workingdir/data_cleaning_clean_and_append_SC_data.R
