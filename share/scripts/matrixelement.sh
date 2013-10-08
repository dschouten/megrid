#!/bin/sh

## site-specific: path to folder containing ME code 
export HEPPY=$HOME/code/higgs

set -e
set -f

if [ "x${MYINPUT}" == "x" ]; then
    echo "ERROR variable MYINPUT not defined"
    exit 1;
fi

if [ "x${MYME}" == "x" ]; then
    echo "ERROR variable MYME not defined"
    exit 1;
fi

if [ "x${MYNUMEVENTS}" == "x" ]; then
    export MYNUMEVENTS=500
fi

if [ "x${MYBLOCKSIZE}" == "x" ]; then
    export MYBLOCKSIZE=2000
fi

if [ "x${MYMAXNJOBS}" == "x" ]; then
    export MYMAXNJOBS=100
fi

if [ "x${MYNICE}" == "x" ]; then
    export MYNICE=0
fi

if [ "x${MYSCRIPT}" == "x" ]; then
    export MYSCRIPT="$HEPPY/mymeanalysis.git/megrid/share/calculateME.py"
fi

numJobsRequired() 
{
    ARGS=( ${@} )
    root -l -q -b "njobs.C( \"${ARGS[0]}\", \"${ARGS[1]}\", ${ARGS[2]} )" 2>&1 | tail -1;
}

## "-o log/ -k oe" "-l nodes=lhcxx.phys.sfu.ca"

submitJobs() {
    nFAIL=0	
    MYRANDOMSEL=${MYNUMEVENTS}
    if [ ${MYNUMEVENTS} -eq $VAR_NEVENTS_PERJOB ]; then
	MYRANDOMSEL=-1
    fi
    for iJOB in $( seq $VAR_NJOBS ); do
	iBEGIN=$(( ($iJOB - 1) * $VAR_NEVENTS_PERJOB ))
	iEND=$(( $iJOB * $VAR_NEVENTS_PERJOB ))
	PBS_PRE_EXEC="source $VAR_LOC/mymeanalysis.git/megrid/share/scripts/loadenv.sh"
	PBS_EXEC="python"
	PBS_MACRO="$VAR_SCRIPT -v ${MYAPPENDFLAG} -t 120 --begin=${iBEGIN} --end=${iEND}" 
        PBS_MACRO="$PBS_MACRO --random=${MYRANDOMSEL} --seed=${iJOB}" 
	PBS_MACRO="$PBS_MACRO --input=${VAR_IN} --output=${VAR_OUT}_$( printf %03d ${iJOB} ).root" 
	PBS_MACRO="$PBS_MACRO --cfg=configs/cfg_${VAR_ME}.py ${MYMEEXEC}"
	export PBS_PRE_EXEC
	export PBS_EXEC
	export PBS_MACRO
	export PBS_JOBDESCR=$VAR_JOBDESCR
	echo $PBS_PRE_EXEC
	echo $PBS_EXEC $PBS_MACRO
	echo "qsub -V $MYQSUB $VAR_LOC/mymeanalysis.git/megrid/share/scripts/job.sh"
	echo -n "job id: "
	qsub -V $MYQSUB $VAR_LOC/mymeanalysis.git/megrid/share/scripts/job.sh && sleep 0.1
	RET=$?
	while [ $RET -ne 0 ]; do
	    nFAIL=$(( ${nFAIL} + 1 ))
	    echo
	    echo "ERROR job not submitted!"
	    echo 
	    echo -n -e "\t retry: "
	    qsub -V $MYQSUB $VAR_LOC/mymeanalysis.git/megrid/share/scripts/job.sh && sleep 0.1
	    RET=$?
	done
    done
    echo 
    echo "${nFAIL} jobs failed to submit"
    echo
}

export VAR_LOC="$HEPPY"
export VAR_SCRIPT=$MYSCRIPT
export VAR_NICE=$MYNICE
export VAR_OPTS=$MYOPTS

export VAR_NJOBS=$( numJobsRequired $MYINPUT "${MYTREE}" $MYBLOCKSIZE )
if [ $VAR_NJOBS -gt $MYMAXNJOBS ]; then
    export VAR_NJOBS=$MYMAXNJOBS
fi

export VAR_ID="$( basename $MYINPUT | sed -e 's:.root::g' )"

# ------------------------------------------------------------------------------------------------

export VAR_IN="${MYINPUT}"
export VAR_ME="${MYME}"
export VAR_NEVENTS_PERJOB=$MYBLOCKSIZE
export VAR_JOBDESCR="${MYME}_${VAR_ID}"  
export VAR_OUT="${MYOUTPUTDIR}/ME_${MYME}_${VAR_ID}"

echo
echo "------------------------------------------------------------------------------------------------"
echo "submitting ${VAR_NJOBS} jobs for ${VAR_IN}"
echo "------------------------------------------------------------------------------------------------"
echo

submitJobs

set +f
set +e
