#!/bin/sh
#PBS -q medium
#PBS -e $HOME/code/gridgen/log/job_${PBS_JOBID}.err
#PBS -o $HOME/code/gridgen/log/job_${PBS_JOBID}.log

#
# configuration parameters:
#
# LHEF_FILE - path to input LHEF file 
# GRID_FILE - path to output grid file
#
# SOURCE_LOCATION - path to where all scripts and configuration files are stored 
#
# EVNT_JOBOPTS - path to job options file to configure parton shower
# EVNT_MAXEVENTS - maximum number of events to read from LHEF file
# SLIM_COMMAND - script for parsing D3PD truth table to finally prepare grid
# RND_SEED - random number seed
#

export SOURCE_LOCATION

export LHEF_FILE
export GRID_FILE

export RND_SEED
export EVNT_JOBOPTS
export EVNT_MAXEVENTS
export SLIM_COMMAND

quit() {
    echo ${@}
    exit -1
}

[ "x${LHEF_FILE}" == "x" ] && quit "no input file defined"

echo "SOURCE_LOCATION=${SOURCE_LOCATION}"
echo "LHEF_FILE=${LHEF_FILE}"
echo "GRID_FILE=${GRID_FILE}"
echo "EVNT_JOBOPTS=${EVNT_JOBOPTS}"
echo "SLIM_COMMAND=${SLIM_COMMAND}"
echo "RND_SEED=${RND_SEED}"

echo
echo " ==> Step #0 ------------------------- prepare work directory ------------------------- "
echo

export WORKDIR=$( mktemp -d /tmp/tmp.XXXXXX )
cd ${WORKDIR}

echo "Working in: ${WORKDIR}"

cp -rvf ${SOURCE_LOCATION}/*.* .

echo
echo " ==> Step #1 ------------------------------ parton shower ----------------------------- "
echo

./evgen.sh 

export EVNT_FILE="$( ls | grep .pool.root )"

echo
echo " ==> Step #2 ------------------------------- D3PD ntuple ------------------------------ "
echo

./d3pd.sh

export D3PD_FILE="$( ls | grep .d3pd.root )"

echo
echo " ==> Step #3 ------------------------------- grid ntuple ------------------------------ "
echo

ROOTSETUP=$( /home/dschouten/bashlib/myroot.sh -v 5.32.02 )
source ${ROOTSETUP}

cp -rv ${SOURCE_LOCATION}/parse/*.py .
python ${SLIM_COMMAND} -f ${GRID_FILE} ${D3PD_FILE}

echo
echo " ==> DONE! "
echo

cd /tmp
rm -rf ${WORKDIR}
