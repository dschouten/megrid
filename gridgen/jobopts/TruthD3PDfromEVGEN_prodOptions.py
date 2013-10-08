# Build truth missing Et
include( "ParticleBuilderOptions/MissingEtTruth_jobOptions.py" )
METPartTruth.TruthCollectionName="GEN_EVENT"
topSequence.METAlg+=METPartTruth

###################################################################
# Make the D3PD.
#

#--------------------------------------------------------------------------
# Configuration
#--------------------------------------------------------------------------
# compatibility with jets
from RecExConfig.RecConfFlags import jobproperties
jobproperties.RecConfFlags.AllowBackNavigation = True

from JetRec.JetRecFlags import jobproperties as jobpropjet
jobpropjet.JetRecFlags.inputFileType = "GEN"

# Build list of particles stored in D3PD
from TruthD3PDAnalysis.truthParticleConfig           import truthParticleConfig
truthParticleConfig (topSequence)

# Build list of particles used to jet building
from TruthD3PDMaker.TruthJetFilterConfig             import TruthJetFilterConfig
TruthJetFilterConfig (topSequence, writePartons = False, writeHadrons = True,
                      excludeWZdecays = True)

# Build truth particle (hadron-level) jets
# flags for AOD
from ParticleBuilderOptions.AODFlags import AODFlags
AODFlags.MissingEtTruth = True
AODFlags.TruthParticleJet = True
AODFlags.McEventKey="GEN_EVENT"

# Reconstruct standard ATLAS truth jets
from JetRec.JetGetters                               import *
antikt4truthAlgStd = make_StandardJetGetter('AntiKt',0.4,'Truth',disable=False).jetAlgorithmHandle()
antikt6truthAlgStd = make_StandardJetGetter('AntiKt',0.6,'Truth',disable=False).jetAlgorithmHandle()

#--------------------------------------------------------------------------
# Make the D3PD
#--------------------------------------------------------------------------
from TruthD3PDMaker.TruthD3PDMakerFlags import TruthD3PDFlags
from OutputStreamAthenaPool.MultipleStreamManager import MSMgr
alg = MSMgr.NewRootStream( "StreamNTUP_TRUTH", TruthD3PDFlags.TruthD3PDOutputFileName(), "truth" )

from EventCommonD3PDMaker.EventInfoD3PDObject        import EventInfoD3PDObject
alg += EventInfoD3PDObject (10,'EventInfo')

from TruthD3PDMaker.GenEventD3PDObject               import GenEventD3PDObject
alg += GenEventD3PDObject (1, 'GenEvent')

from TrackD3PDMaker.TruthVertexD3PDObject            import TruthVertexD3PDObject
from D3PDMakerConfig.D3PDMakerFlags                  import D3PDMakerFlags
alg += TruthVertexD3PDObject(1, 'TruthVertex', sgkey=D3PDMakerFlags.TruthSGKey(),prefix='mcVx',label='mcVx_')
alg += TruthVertexD3PDObject(1, 'TruthVertexFake', sgkey=D3PDMakerFlags.TruthSGKey(),prefix='vx',label='vx_')

from JetD3PDMaker.JetD3PDObject                      import JetD3PDObject

alg += JetD3PDObject (2, 'Jet_Truth', sgkey='AntiKt4TruthJets', prefix='jet_akt4_', include=['TrueFlavorComponents'],
                      exclude=['L1Kinematics','L2Kinematics','EFKinematics','El02Match','Mu02Match'])
alg += JetD3PDObject (2, 'Jet_NotTruth', sgkey='AntiKt4TruthJets', prefix='jet_antiKt4Truth_',
                      include=['TrueFlavorComponents'],
                      exclude=['L1Kinematics','L2Kinematics','EFKinematics','El02Match','Mu02Match'])

alg += JetD3PDObject (2, 'Jet_Truth', sgkey='AntiKt6TruthJets', prefix='jet_akt6_', include=['TrueFlavorComponents'],
                      exclude=['L1Kinematics','L2Kinematics','EFKinematics','El02Match','Mu02Match'])
alg += JetD3PDObject (2, 'Jet_NotTruth', sgkey='AntiKt6TruthJets', prefix='jet_antiKt6Truth_',
                      include=['TrueFlavorComponents'],
                      exclude=['L1Kinematics','L2Kinematics','EFKinematics','El02Match','Mu02Match'])

from SUSYD3PDMaker.SUSYTruthParticleD3PDObject       import SUSYTruthParticleD3PDObject
alg += SUSYTruthParticleD3PDObject (1, 'SUSYTruthParticle')

from TruthD3PDMaker.Atlfast1MissingETD3PDObject      import TruthMETD3PDObject
alg += TruthMETD3PDObject (level=10)

#from MissingETD3PDMaker.MissingETD3PDMakerFlags      import MissingETD3PDMakerFlags
#MissingETD3PDMakerFlags.doMissingETRegions=False
#from MissingETD3PDMaker.MissingETD3PDObject          import RefFinalMETD3PDObject
#alg += RefFinalMETD3PDObject( 10 , name='MET_Simplified20' , lebel=1, sgkey='MET_Truth_NonInt' , prefix='MET_Simplified20_' , label='MET_Simplified20_' , exclude=['MET_Regions'] , include =  ['MET','MET_Phi','MET_Et','MET_SumEt'] )

from TruthD3PDMaker.TruthParticleFakerObject         import *
alg += TruthParticleFakerElD3PDObject (level=10, sgkey='GEN_EVENT', prefix="el_", label='el_')
alg += TruthParticleFakerElD3PDObject (level=10, sgkey='GEN_EVENT', prefix="el_truth_", label='el_t_')
alg += TruthParticleFakerMuD3PDObject (level=10, sgkey='GEN_EVENT', prefix="mu_", label='mu_')
alg += TruthParticleFakerMuD3PDObject (level=10, sgkey='GEN_EVENT', prefix="mu_truth_", label='mu_t_')
alg += TruthParticleFakerPhD3PDObject (level=10, sgkey='GEN_EVENT', prefix="ph_", label='ph_')
alg += TruthParticleFakerPhD3PDObject (level=10, sgkey='GEN_EVENT', prefix="ph_truth_", label='ph_t_')
alg += TruthParticleFakerTauD3PDObject (level=10, sgkey='GEN_EVENT', prefix="tau_", label='tau_')
alg += TruthParticleFakerTauD3PDObject (level=10, sgkey='GEN_EVENT', prefix="trueTau_", label='tau_t_')

