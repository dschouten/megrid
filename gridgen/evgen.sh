#!/bin/sh

mkdir evgen/ && cd evgen

DSNAME=$( basename ${LHEF_FILE//.tar.gz/} )
DSFOLDER=$( dirname ${LHEF_FILE} )

# configure Athena environment

source /atlas/ATLASLocalRootBase/user/atlasLocalSetup.sh
source $AtlasSetup/scripts/asetup.sh 17.2.8.3,AtlasProduction,here 

export JOBOPTSEARCHPATH=/cvmfs/atlas.cern.ch/repo/sw/Generators/MC12JobOptions/latest/common:$JOBOPTSEARCHPATH

# create the input files

cp ${LHEF_FILE} .

# create the job options for Athena

cp ${SOURCE_LOCATION}/${EVNT_JOBOPTS} joboptions.py

Generate_trf.py --omitvalidation=testEventMinMax ecmEnergy=8000 firstEvent=0 maxEvents=${EVNT_MAXEVENTS} randomSeed=${RND_SEED} \
    jobConfig=joboptions.py runNumber=1 outputEVNTFile=$PWD/${DSNAME}.pool.root \
    inputGeneratorFile=$( basename ${LHEF_FILE} ) &> evgen.log

echo "random: ${RND_SEED}" >> evgen.log

mv ${DSNAME}.pool.root ${WORKDIR}/.
tar czf ${WORKDIR}/evgen.log.tgz evgen.log 
