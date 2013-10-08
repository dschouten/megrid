# @author doug schouten <doug dot schouten at triumf dot ca>

from math import sqrt

#
# this snippet is responsible for reading in the
# BDT regression estimate for recoil px, py
#

class BDTConfig:
    def __init__( self, weights, varstring, index=-1 ):
        self.weights   = weights
        self.varstring = varstring
        self.index = index

if 'tmvaInitialized' not in vars():
    configs = {
        'ggf_0j' : BDTConfig( 'recoilweights_ggf_0j.xml', 'dilep_pt:dilep_pxpy:trackmet_pxpy:calomet_pxpy:numprimvertices:mu' ),
        'ww_0j'  : BDTConfig( 'recoilweights_ww_0j.xml',  'dilep_pt:dilep_pxpy:trackmet_pxpy:calomet_pxpy:numprimvertices:mu' ),
        'ggf_Nj' : BDTConfig( 'recoilweights_ggf_Nj.xml', 'dilep_pt:dilep_pxpy:totjet_pxpy:totjet_scalarsumpt:trackmet_pxpy:calomet_pxpy:numjets:numprimvertices:mu' ),
        'ww_Nj'  : BDTConfig( 'recoilweights_ww_Nj.xml' , 'dilep_pt:dilep_pxpy:totjet_pxpy:totjet_scalarsumpt:trackmet_pxpy:calomet_pxpy:numjets:numprimvertices:mu' ) }
    
    import ROOT

    ROOT.gSystem.Load( 'libTMVA.so' )
    ROOT.gSystem.Load( 'libTMVAWrapper.so' )
    index = 0
    for (k, v) in configs.items():
        if os.path.exists( v.weights ):
            log.info( 'loading BDTG estimate from [%s], assigned to index #%d'%( v.weights, index ) )
            ROOT.TMVATools.loadReader( 'BDT::BDTG', v.weights, v.varstring, index )
            v.index = index
            index += 1
        else:
            continue
    tmvaInitialized = True

from math import cos
from math import sin

lepPx0 = itree.lepPt0 * cos( itree.lepPhi0 )
lepPy0 = itree.lepPt0 * sin( itree.lepPhi0 )

lepPx1 = itree.lepPt1 * cos( itree.lepPhi1 )
lepPy1 = itree.lepPt1 * sin( itree.lepPhi1 )

sumJetPx = 0.
sumJetPy = 0.
sumJetPT = 0.

numJets = itree.m_jet_n

for ijet in xrange( numJets ):
    try:
        sumJetPx += itree.m_jet_pt[ijet] * cos( itree.m_jet_phi[ijet] )
        sumJetPy += itree.m_jet_pt[ijet] * sin( itree.m_jet_phi[ijet] )
        sumJetPT += itree.m_jet_pt[ijet]
    except AttributeError:
        if ijet <= 2:
            sumJetPx += getattr( itree, 'jetPt%d'%( ijet ) ) * cos( getattr( itree, 'jetPhi%d'%( ijet ) ) )
            sumJetPy += getattr( itree, 'jetPt%d'%( ijet ) ) * sin( getattr( itree, 'jetPhi%d'%( ijet ) ) )
            sumJetPT += getattr( itree, 'jetPt%d'%( ijet ) ) 
        else:
            continue

if numJets == 0:
    args_x = [ itree.Ptll/GeV            , lepPx0/GeV + lepPx1/GeV,
               itree.MET_x_TrackHWW/GeV  , itree.MET_x/GeV,
               itree.Nvxp                , itree.averageIntPerXing ]
    
    args_y = [ itree.Ptll/GeV            , lepPy0/GeV + lepPy1/GeV,
               itree.MET_y_TrackHWW/GeV  , itree.MET_y/GeV,
               itree.Nvxp                , itree.averageIntPerXing ]
else:
    args_x = [ itree.Ptll/GeV           , lepPx0/GeV + lepPx1/GeV,
               sumJetPx/GeV             , sumJetPT/GeV,
               itree.MET_x_TrackHWW/GeV , itree.MET_x/GeV,
               numJets                  , itree.Nvxp,
               itree.averageIntPerXing ]
    
    args_y = [ itree.Ptll/GeV           , lepPy0/GeV + lepPy1/GeV,
               sumJetPy/GeV             , sumJetPT/GeV,
               itree.MET_y_TrackHWW/GeV , itree.MET_y/GeV,
               numJets                  , itree.Nvxp,
               itree.averageIntPerXing ]

key = None

if mehndl.get_name().lower().startswith( 'hww' ) or \
   mehndl.get_name().lower().startswith( 'gww' ) or \
   mehndl.get_name().lower().startswith( 'ggf' ):
    key = 'ggf'
elif mehndl.get_name().lower().startswith( 'ww' ):
    key = 'ww'

if key == None:
    raise RuntimeError

if numJets == 0: key += '_0j'
else           : key += '_Nj'

reader_index = configs[key].index

px = ROOT.TMVATools.evaluateRegressionMVA( reader_index, *args_x )
py = ROOT.TMVATools.evaluateRegressionMVA( reader_index, *args_y )

measured_recoil.SetPxPyPzE( -px, -py, 0, sqrt(px**2 + py**2) ) ## note the sign !! recoil = -higgs pT

log.info( 'recoil from R-BDT %s estimate (pT, phi): %0.2f, %0.2f'%( key,
                                                                    measured_recoil.Pt(),
                                                                    measured_recoil.Phi() ) )
