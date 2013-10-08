tupleFileOutput = '__OUTPUT__'

###################################################################
# Define the input file here.
#

from AthenaCommon.AthenaCommonFlags import athenaCommonFlags
athenaCommonFlags.FilesInput= ["__INPUT__"]

###################################################################
# Define the output file here.
#
from TruthD3PDMaker.TruthD3PDMakerFlags import TruthD3PDFlags
if not globals().get('tupleFileOutput'):
    tupleFileOutput = '__OUTPUT__'
TruthD3PDFlags.TruthD3PDOutputFileName = tupleFileOutput

athenaCommonFlags.EvtMax = -1

from RecExConfig.RecFlags import rec

rec.doHist = False
rec.doMonitoring = False
rec.doCBNT = False
rec.doWriteESD = False
rec.doWriteAOD = False
rec.doWriteTAG = False

rec.doPerfMon.set_Value_and_Lock( False )
rec.doDetailedPerfMon = False
rec.doSemiDetailedPerfMon = False

from PerfMonComps.PerfMonFlags import jobproperties
jobproperties.PerfMonFlags.doMonitoring.set_Value_and_Lock( False )

from D3PDMakerConfig.D3PDMakerFlags import D3PDMakerFlags
D3PDMakerFlags.TruthSGKey = 'GEN_EVENT,GEN_AOD,TruthEvent'
D3PDMakerFlags.DoTrigger = False

from AthenaCommon.GlobalFlags import globalflags
globalflags.DetGeo.set_Value_and_Lock('atlas')
globalflags.ConditionsTag.set_Value_and_Lock('OFLCOND-SDR-BS7T-04-13')
globalflags.DetDescrVersion.set_Value_and_Lock("ATLAS-GEO-16-00-00")

rec.doDPD=True
rec.readAOD.set_Value_and_Lock(True)
rec.AutoConfiguration.set_Value_and_Lock(['ProjectName','BeamType','RealOrSim','DoTruth','InputType'])

rec.doFileMetaData.set_Value_and_Lock(False)

from JetRec.JetRecFlags import jetFlags
jetFlags.noStandardConfig.set_Value_and_Lock(True)
jetFlags.evgenJets.set_Value_and_Lock(True)

rec.DPDMakerScripts.append("TruthD3PDfromEVGEN_prodOptions.py")

include("RecExCommon/RecExCommon_topOptions.py")
