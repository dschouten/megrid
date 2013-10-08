
export MYMAXNJOBS=9 ## max number of sub-jobs per dataset
export MYNICE=10 ## niceness level for batch jobs

export MYQSUB="-q medium" ## qsub directives

export MYOPTS="" ## an environment variable the running job may find useful

export MYOUTPUTDIR=/data/ds03_1/dschouten/outputs/fastme/testing ## choose a different location for different runs

#########################################################################

## for running over MVA ntuples
export MYSCRIPT="$HEPPY/mymeanalysis.git/megrid/share/calculateME.py"
export MYTREE="HWWTree"
export MEEXEC=""

#########################################################################

MELIST="WW,ggF125" 

ALLINPUTS="/global/dschouten/inputs/00-02-07/mva/OS/ww_0j.root
/global/dschouten/inputs/00-02-07/mva/OS/ggf125_0j.root"

# ALLINPUTS="/global/dschouten/inputs/00-02-17/vbfbdt/Nominal/vbf125_slim.root 
# /global/dschouten/inputs/00-02-17/vbfbdt/Nominal/ttbar_slim.root" 

export MYNUMEVENTS_DEF=250  # default number of events to process in each sub-job
export MYBLOCKSIZE_DEF=2000 # default block of events to assign to each sub-job 
                            # (allows to process < 100% of input events randomly)

export MYNUMEVENTS_MIN=250
export MYBLOCKSIZE_MIN=2000
