#=================================================================================================
#
# from MC12.110202.PowhegPythia8_AU2CT10_ttbar_LeptonFilter.py
#

evgenConfig.description = "POWHEG+Pythia8 ttbar production with a 1 MeV lepton filter and AUET2B CT10 tune"
evgenConfig.keywords = ["top" ,"ttbar", "leptonic", "powheg", "pythia8" ]
evgenConfig.inputfilecheck = 'ttbar' 

include("MC12JobOptions/PowhegPythia8_AU2_CT10_Common.py")

include("MC12JobOptions/Pythia8_Tauola.py")
include("MC12JobOptions/Pythia8_Photos.py")

include("MC12JobOptions/TTbarWToLeptonFilter.py")
topAlg.TTbarWToLeptonFilter.Ptcut = 1000.0 # 1 GeV
topAlg.TTbarWToLeptonFilter.NumLeptons = 2

#=================================================================================================
