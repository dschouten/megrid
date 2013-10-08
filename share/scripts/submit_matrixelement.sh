#!/bin/sh

source $( $HOME/bashlib/myroot.sh )

export HEPPY=$HOME/code/higgs

numEvents() 
{
    ARGS=( ${@} )
    root -l -q -b "nevents.C( \"${ARGS[0]}\", \"${ARGS[1]}\" )" 2>&1 | tail -1;
}

# ------------------------------------------------------------------------------------------------

export MYMAXNJOBS=1000 ## max number of sub-jobs per dataset
export MYNICE=10 ## niceness level for batch jobs

export MYQSUB="-q long" ## qsub directives

export MYOPTS="" ## an environment variable the running job may find useful

export MYOUTPUTDIR=/data/ds03_1/dschouten/outputs/ ## choose a different location for different runs

#########################################################################

## for running over MVA ntuples
export MYSCRIPT="$HEPPY/mymeanalysis.git/megrid/share/calculateME.py"
export MYTREE="HWWTree"
export MEEXEC=""

#########################################################################

MELIST="XYZ" ## for cfg_XYZ config files 

ALLINPUTS="ntuple.root"

export MYNUMEVENTS_DEF=1000 # default number of events to process in each sub-job
export MYBLOCKSIZE_DEF=1000 # default block of events to assign to each sub-job 
                            # (allows to process < 100% of input events randomly)

export MYNUMEVENTS_MIN=100
export MYBLOCKSIZE_MIN=100

# ------------------------------------------------------------------------------------------------

[ -f config_matrixelement.sh ] && source config_matrixelement.sh

# ------------------------------------------------------------------------------------------------

mkdir -p ${MYOUTPUTDIR}
if [ $? -ne 0 ]; then
    exit 1
fi

for INPUT in ${ALLINPUTS//,/ }; do 
    export MYNUMEVENTS=${MYNUMEVENTS_DEF}
    export MYBLOCKSIZE=${MYBLOCKSIZE_DEF}
    
    export MYINPUT=${INPUT}
    
    NEVNTS=$( numEvents ${INPUT} "${MYTREE}" )

    if [ ${NEVNTS} -lt $(( ${MYBLOCKSIZE} * ${MYMAXNJOBS} )) ]; then
	export MYNUMEVENTS=${MYNUMEVENTS_MIN}
	export MYBLOCKSIZE=${MYBLOCKSIZE_MIN}
    fi
    
    iME=0
    for ME in ${MELIST//,/ }; do
	if [ ${iME} -eq 0 ]; then
	    export MYAPPENDFLAG="-a"
	    iME=1
	else
	    export MYAPPENDFLAG=""
	fi
	export MYME=${ME}
	export MYMEEXEC=${MEEXEC}
	./matrixelement.sh
    done
    # break
done
