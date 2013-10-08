#=================================================================================================
#
# from MC12.1610NN.PowhegPythia8_AU2CT10_VBFHXYZ_WW2lep_EF_15_5.py
#
evgenConfig.description = "POWHEG+PYTHIA8 VBF H->WW->lvlv with AU2 CT10"
evgenConfig.keywords = ["SMhiggs", "VBF", "W","leptonic"]
evgenConfig.inputfilecheck = "VBFH_SM_M125"

include("MC12JobOptions/PowhegPythia8_AU2_CT10_Common.py")
# ... Photos
include ( "MC12JobOptions/Pythia8_Photos.py" )

topAlg.Pythia8.Commands += [
                            '25:onMode = off',#decay of Higgs
                            '25:onIfMatch = 24 24',
                            '24:onMode = off',#decay of W
                            '24:mMin = 2.0',
                            '24:onIfMatch = 11 12',
                            '24:onIfMatch = 13 14',
                            '24:onIfMatch = 15 16'
                            ]

# ... Filter
from GeneratorFilters.GeneratorFiltersConf import MultiLeptonFilter
topAlg += MultiLeptonFilter("Multi1TLeptonFilter")
topAlg += MultiLeptonFilter("Multi2LLeptonFilter")

Multi1TLeptonFilter = topAlg.Multi1TLeptonFilter
Multi1TLeptonFilter.Ptcut = 15000.
Multi1TLeptonFilter.Etacut = 5.0
Multi1TLeptonFilter.NLeptons = 1

Multi2LLeptonFilter = topAlg.Multi2LLeptonFilter
Multi2LLeptonFilter.Ptcut = 5000.
Multi2LLeptonFilter.Etacut = 5.0
Multi2LLeptonFilter.NLeptons = 2

StreamEVGEN.RequireAlgs  = [ "Multi1TLeptonFilter" ]
StreamEVGEN.RequireAlgs += [ "Multi2LLeptonFilter" ]

#=================================================================================================
