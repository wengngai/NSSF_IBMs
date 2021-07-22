##!/bin/bash

PATH=/share/apps/R/R-3.6.2/bin/:/share/apps/gcc-8.3/bin:$PATH
LD_LIBRARY_PATH=/share/apps/gcc-8.3/lib64:$LD_LIBRARY_PATH

###############################
# SGE settings Here
# Basically, if a line starts with "#$", then you can enter any
# qsub command line flags .  See the qsub(1) man page.
# Email address, and options, to whom reports should be sent.
# The -M is for the address, the -m is for when to send email.
# a=abort b=begin e=end s=suspend n=none
#$ -M haoran.lai@canterbury.ac.nz
#$ -m ae

# Redirect the STDOUT and STDERR files to the ~/jobs directory
#$ -o /home/hrl23/jobs/
#$ -e /home/hrl23/jobs/

# This script, ladies and gentlemen, is in bash
#$ -S /bin/bash

# Do some validation checking, and bail on errors
##$ -w e

# Operate in the current directory
#$ -cwd

# Jobs variable
#$ -t 1

# Number of cores
##$ -pe mpi_16_tasks_per_node 16

# End SGE Settings
###############################
#you can put your scripts here

# /share/apps/R/R-3.6.2/bin/Rscript code/Prunus_IBM_v1.R $SGE_TASK_ID
/share/apps/R/R-3.6.2/bin/Rscript code/four_spp_IBM.R $SGE_TASK_ID


# run: qsub -S /bin/bash -pe smp 24 code/run_sim.sh