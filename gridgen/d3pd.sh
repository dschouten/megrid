#!/bin/sh

mkdir d3pd/ && cd d3pd

# get the input file

ln -s ${WORKDIR}/${EVNT_FILE} .

# configure Athena environment

source /atlas/ATLASLocalRootBase/user/atlasLocalSetup.sh
source $AtlasSetup/scripts/asetup.sh 17.2.2,AtlasOffline,here 

# create joboptions files

sed -e "s:__INPUT__:${EVNT_FILE}:g" -e "s:__OUTPUT__:${EVNT_FILE//.pool.root/.d3pd.root}:g" \
    ${SOURCE_LOCATION}/jobopts/TruthD3PDfromEVGEN_topOptions.py > joboptions.py

cp ${SOURCE_LOCATION}/jobopts/TruthD3PDfromEVGEN_prodOptions.py .

athena.py -b -l INFO joboptions.py &> d3pd.log

tar czf ../d3pd.log.tgz d3pd.log

mv -v ${EVNT_FILE//.pool.root/.d3pd.root} ${WORKDIR}/.
