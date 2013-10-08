#=================================================================================================
#
# from MC12.110004.PowhegJimmy_AUET2CT10_ttbar_dilepton.py
#

evgenConfig.description = "Powheg+fHerwig/Jimmy ttbar production with a two lepton filter and AUET2 CT10 tune"
evgenConfig.keywords = ["top", "ttbar", "dileptonic"]
evgenConfig.contact  = ["doug.schouten@cern.ch"]
evgenConfig.inputfilecheck = 'ttbar' 

include("MC12JobOptions/PowhegJimmy_AUET2_CT10_Common.py")
evgenConfig.generators += [ "Lhef"]

include("MC12JobOptions/Jimmy_Tauola.py")
include("MC12JobOptions/Jimmy_Photos.py")

include("MC12JobOptions/TTbarWToLeptonFilter.py")
topAlg.TTbarWToLeptonFilter.Ptcut = 1000.0 # 1 GeV
topAlg.TTbarWToLeptonFilter.NumLeptons = 2

#=================================================================================================
