#!/bin/bash 
#PBS -r n

### Output files
#PBS -e $HOME/pbs/log/job_${PBS_JOBID}.err
#PBS -o $HOME/pbs/log/job_${PBS_JOBID}.log

### Choose queue
#PBS -q medium

### Mail to user
# - - _BS -m ae
# - - _BS -M doug.schouten@triumf.ca

set -f

#--------------------------------------------------------
# Pass arguments here
#--------------------------------------------------------
export WORKDIR=/tmp/${PBS_JOBID}/

#--------------------------------------------------------
# Operational Settings here
#--------------------------------------------------------

#---------------------------------------------
# Do the execution
#---------------------------------------------
mkdir -p ${WORKDIR}
mkdir -p ${WORKDIR}/log
mkdir -p ${WORKDIR}/outputs

cd $WORKDIR

$PBS_PRE_EXEC | tee -a log/pre_exec.log

$PBS_EXEC $PBS_MACRO | tee -a log/exec.log

#---------------------------------------------
# Copy and Clean up
#---------------------------------------------
## mkdir -p $PBS_OUTPUTAREA/$PBS_JOBID
## cp -rf $WORKDIR/outputs/* $PBS_OUTPUTAREA/$PBS_JOBID

cd; rm -rf ${WORKDIR}
#--------------------------------------------------------

set +f
