#!/bin/sh

export SOURCE_LOCATION=$PWD

# ==========================================================================================

# # -------- Wt Powheg + Pythia
# export EVNT_JOBOPTS="jobopts/Wt_PY8.py"
# export EVNT_MAXEVENTS=30000
# export SLIM_COMMAND="wt.py -v --ps=PYTHIA8 --dojets --leptons=e,mu"
# export SEARCH_TAG="*/*.110141.st_Wtchan_dilepton_DR_8TeV.*/*.tar.gz"
# export OUTPUT_FOLDER="/global/dschouten/grids/outputs/Wt/"

# -------- tT Powheg + Pythia
# export EVNT_JOBOPTS="jobopts/tT_PY8.py"
# export EVNT_MAXEVENTS=30000
# export SLIM_COMMAND="ttbar.py -v --ps=PYTHIA8 --dojets --leptons=e,mu"
# export SEARCH_TAG="*/*.110004.ttbar_dilepton_8TeV.*/*.tar.gz"
# export OUTPUT_FOLDER="/global/dschouten/grids/outputs/tT/"

# -------- ggF Powheg + Pythia
export EVNT_JOBOPTS="jobopts/ggF_PY8.py"
export EVNT_MAXEVENTS=30000
export SLIM_COMMAND="ww.py -v --ps=PYTHIA8 --dojets --leptons=e,mu"
export SEARCH_TAG="*/*.116703.ggH_SM_M125_8TeV.*/*.tar.gz"
export OUTPUT_FOLDER="/global/dschouten/grids/outputs/ggF/"

# ==========================================================================================

mkdir -p ${OUTPUT_FOLDER}

IFILE=0
for FILE in $( find /global/dschouten/grids/inputs -wholename "${SEARCH_TAG}" | sort ); do
    IFILE=$(( ${IFILE} + 1 ))

    if [ ${IFILE} -le 1 ]; then
    	continue
    fi

    export LHEF_FILE=${FILE}
    export GRID_FILE=${OUTPUT_FOLDER}/$( basename ${FILE//.tar.gz/.grid.root} )
    export RND_SEED=${IFILE}

    echo "${LHEF_FILE} -> ${GRID_FILE}"

    if [ -f ${GRID_FILE} ]; then
	continue;
    fi
    
    qsub -V gridgen.sh
    while [ $? -ne 0 ]; do
    	qsub -V gridgen.sh
    done

    # if [ ${IFILE} -ge 5 ]; then
    # 	break
    # fi
done

