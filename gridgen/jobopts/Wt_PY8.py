#=================================================================================================
#
# from MC12.110141.PowhegPythia_P2011C_st_Wtchan_dilepton_DR.py
#

evgenConfig.description = "POWHEG+Pythia6 Wt dilepton production with PythiaPerugia2011C tune"
evgenConfig.keywords = ["top", "Wt", "dilepton", "DR"]
evgenConfig.contact  = ["doug.schouten@cern.ch"]
evgenConfig.inputfilecheck = 'st_Wtchan_dilepton_DR'

include("MC12JobOptions/PowhegPythia_Perugia2011C_Common.py")

include("MC12JobOptions/Pythia_Tauola.py")
include("MC12JobOptions/Pythia_Photos.py")

#=================================================================================================
