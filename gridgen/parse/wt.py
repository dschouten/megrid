#!/bin/env python
#
#
# wt.py - create grid from events with tt->W(lv)b,W(lv)b final states
#
#
# USAGE: wt.py --ps= --leptons= [-f] output.root input.root [... input.root]
#
# --ps       : identify parton shower program used (Pythia6/Pythia8/Herwig/Herwig++)
# --leptons  : list of allowed lepton decays (e,mu,tau)
# -f,--force : force overwrite if output file exists
#

import os
import sys
import math
import random
import logging
import string

from optparse import OptionParser
from math import sqrt
from array import array

#-------------------------------------------------------------------------------------------------------#

def usage():
    print 'ww.py [-v,--verbose] [-w,--warn] [-f,--force] --leptons= --ps= output.root input.root [input.root ... input.root]'
    sys.exit( -1 )

parser = OptionParser()
parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='verbose output' )
parser.add_option('-w', '--warn', dest='warn', action='store_true', default=False, help='halt on warn (and delete the output file)' )
parser.add_option('-f', '--force', dest='force', action='store_true', default=False, help='overwrite output file if exists already' )
parser.add_option('-t', '--tree', dest='tree', action='store', type='string', default='truth', help='name of tree in input file(s)' )
parser.add_option('-p', '--ps', dest='ps', action='store', type='string', default='herwig', help='parton shower program used' )
parser.add_option('-l', '--leptons', dest='leptons', action='store', type='string', default='e,mu,tau', help='allowed lepton decays' )
parser.add_option('-j', '--dojets', dest='dojets', action='store_true', default=False, help='flag to save jet momenta' )

(options, args) = parser.parse_args()

if len( args ) < 2:
    usage()

#-------------------------------------------------------------------------------------------------------#
    
import ROOT
ROOT.gROOT.SetBatch( True )

from truthparser import *

#-------------------------------------------------------------------------------------------------------#

def warn( msg ):
    global options
    logging.warn( msg )
    if options.warn:
        sys.exit( -1 )
    pass

def error( msg ):
    logging.error( msg )
    sys.exit( -1 )
    pass

lepdecays = [ PDG[p] for p in map( string.lower, options.leptons.split( ',' ) ) ]

logging.basicConfig( level='INFO' )

if options.verbose:
    logging.basicConfig( level='DEBUG' )

parser = None
if options.ps.lower() == 'pythia6':
    parser = parse_wt_pythia
elif options.ps.lower() == 'pythia8':
    parser = parse_wt_pythia
elif options.ps.lower() == 'herwig':
    parser = parse_wt_herwig

if parser == None:
    error( 'invalid PS or not implemented' )

matchb    = True
minjetpt  = 15.0
maxjeteta = 5.0

#-------------------------------------------------------------------------------------------------------#

ofilen, ifilenlist = args[0], args[1:]

if os.path.exists( ofilen ) and not options.force:
    error( 'output file [%s] exists'%( ofilen ) )
else:
    warn( 'overwriting existing file [%s]'%( ofilen ) )
        
ohndl = ROOT.TFile.Open( ofilen, 'recreate' )

if ohndl == None or ohndl.IsZombie():
    error( 'could not create output file [%s]'%( ofilen ) )

if not os.path.exists( 'libvec_cc.so' ):
    hndl = open( 'libvec.cc', 'w' )
    hndl.write( '#include "TLorentzVector.h"\n' )
    hndl.write( '#ifdef __CINT__\n' )
    hndl.write( '#pragma link C++ class vector<vector<int> >+;\n' )
    hndl.write( '#pragma link C++ class vector<TLorentzVector>+;\n' )
    hndl.write( '#endif\n' )
    hndl.close()
ROOT.gROOT.ProcessLine( '.L libvec.cc+' )

itree = ROOT.TChain( options.tree, options.tree )
for ifilen in ifilenlist:
    itree.Add( ifilen )

if itree.GetEntries() == 0:
    error( 'no input data to parse' )

itree.SetBranchStatus( '*', 0 )
itree.SetBranchStatus( 'mc_*', 1 )
itree.SetBranchStatus( 'mcevt_*', 1 )

if options.dojets:
    itree.SetBranchStatus( 'jet_antiKt4Truth_n',   1 )
    itree.SetBranchStatus( 'jet_antiKt4Truth_E',   1 )
    itree.SetBranchStatus( 'jet_antiKt4Truth_pt',  1 )
    itree.SetBranchStatus( 'jet_antiKt4Truth_eta', 1 )
    itree.SetBranchStatus( 'jet_antiKt4Truth_phi', 1 )

if options.verbose:
    itree.Show()

ohndl.cd()

otree = ROOT.TTree( 'grid', '' )

br_lp_px = array( 'd', [0] )
br_lm_px = array( 'd', [0] )
br_lp_py = array( 'd', [0] )
br_lm_py = array( 'd', [0] )
br_lp_pz = array( 'd', [0] )
br_lm_pz = array( 'd', [0] )
br_lp_id = array( 'i', [0] )
br_lm_id = array( 'i', [0] )

br_vl_px = array( 'd', [0] )
br_vr_px = array( 'd', [0] )
br_vl_py = array( 'd', [0] )
br_vr_py = array( 'd', [0] )
br_vl_pz = array( 'd', [0] )
br_vr_pz = array( 'd', [0] )

br_b_px = array( 'd', [0] )
br_b_py = array( 'd', [0] )
br_b_pz = array( 'd', [0] )
br_b_m = array( 'd', [0] )
br_b_id = array( 'i', [0] )

br_qsqr = array( 'd', [0] )
br_x1 = array( 'd', [0] )
br_x2 = array( 'd', [0] )
br_fl1 = array( 'd', [0] )
br_fl2 = array( 'd', [0] )

otree.Branch( 'lp_px', br_lp_px, 'lp_px/D' )
otree.Branch( 'lm_px', br_lm_px, 'lm_px/D' )
otree.Branch( 'lp_py', br_lp_py, 'lp_py/D' )
otree.Branch( 'lm_py', br_lm_py, 'lm_py/D' )
otree.Branch( 'lp_pz', br_lp_pz, 'lp_pz/D' )
otree.Branch( 'lm_pz', br_lm_pz, 'lm_pz/D' )
otree.Branch( 'lp_id', br_lp_id, 'lp_id/I' )
otree.Branch( 'lm_id', br_lm_id, 'lm_id/I' )

otree.Branch( 'vl_px', br_vl_px, 'vl_px/D' )
otree.Branch( 'vr_px', br_vr_px, 'vr_px/D' )
otree.Branch( 'vl_py', br_vl_py, 'vl_py/D' )
otree.Branch( 'vr_py', br_vr_py, 'vr_py/D' )
otree.Branch( 'vl_pz', br_vl_pz, 'vl_pz/D' )
otree.Branch( 'vr_pz', br_vr_pz, 'vr_pz/D' )

otree.Branch( 'b_px', br_b_px, 'b_px/D' )
otree.Branch( 'b_py', br_b_py, 'b_py/D' )
otree.Branch( 'b_pz', br_b_pz, 'b_pz/D' )
otree.Branch( 'b_m', br_b_m, 'b_m/D' )
otree.Branch( 'b_id', br_b_id, 'b_id/I' )

otree.Branch( 'x1', br_x1, 'x1/D' )
otree.Branch( 'x2', br_x2, 'x2/D' )
otree.Branch( 'fl1', br_fl1, 'fl1/D' )
otree.Branch( 'fl2', br_fl2, 'fl2/D' )
otree.Branch( 'qsqr', br_qsqr, 'qsqr/D' )

br_njets = array( 'i', [0] )
br_jet_id = ROOT.std.vector('int')()
br_jet_px = ROOT.std.vector('double')()
br_jet_py = ROOT.std.vector('double')()
br_jet_pz = ROOT.std.vector('double')()

if options.dojets:    
    otree.Branch( 'm_jet_n', br_njets, 'm_jet_n/I' )
    otree.Branch( 'jet_px', br_jet_px )
    otree.Branch( 'jet_py', br_jet_py )
    otree.Branch( 'jet_pz', br_jet_pz )
    if matchb:
        otree.Branch( 'jet_id', br_jet_id )

#-------------------------------------------------------------------------------------------------------#

nevents = itree.GetEntries()

#
# use TTree::Draw() to find all the particle indices for tops
# to speed up the processing (avoid long Python loops)
#

indexmap = { }
    
nvals = itree.Draw( 'Entry$:Iteration$', 'abs(mc_pdgId)==%d || abs(mc_pdgId)==%d'%(PDG['t'],PDG['W']), 'goff' )
itree.GetPlayer().SetEstimate( nvals )
nvals = itree.Draw( 'Entry$:Iteration$', 'abs(mc_pdgId)==%d || abs(mc_pdgId)==%d'%(PDG['t'],PDG['W']), 'goff' )

evt_indices = [ int(itree.GetVal( 0 )[i]) for i in xrange( nvals ) ]
par_indices = [ int(itree.GetVal( 1 )[i]) for i in xrange( nvals ) ]

for iind, evt in enumerate( evt_indices ):
    if evt not in indexmap:
        indexmap[evt] = []
    indexmap[evt].append( par_indices[iind] )

print 'will run over [%d] events'%( nevents ) 
sys.stdout.flush()

for ievt in xrange( nevents ):
    itree.GetEntry( ievt )

    if ievt > 0 and ievt % ( nevents / 20 ) == 0:
        print 'processed %d entries'%( ievt )
        sys.stdout.flush()

    ( tlv_b, tlv_lm, tlv_vr, tlv_lp, tlv_vl,
     pdg_lm, pdg_lp, pdg_b ) = parser( itree, lepdecays = lepdecays, version=options.ps.upper(), indices = indexmap[ievt] )

    if None in [ tlv_b, tlv_lm, tlv_vr, tlv_lp, tlv_vl, pdg_lm, pdg_lp, pdg_b ]:
        print 'no t->W,b or W->l,v signature found in event [%d]'%( ievt ) 
        continue
    
    br_lp_px[0] = tlv_lp.Px()
    br_lm_px[0] = tlv_lm.Px()
    br_lp_py[0] = tlv_lp.Py()
    br_lm_py[0] = tlv_lm.Py()
    br_lp_pz[0] = tlv_lp.Pz()
    br_lm_pz[0] = tlv_lm.Pz()
    br_lp_id[0] = pdg_lp
    br_lm_id[0] = pdg_lm
    
    br_vl_px[0] = tlv_vl.Px()
    br_vr_px[0] = tlv_vr.Px()
    br_vl_py[0] = tlv_vl.Py()
    br_vr_py[0] = tlv_vr.Py()
    br_vl_pz[0] = tlv_vl.Pz()
    br_vr_pz[0] = tlv_vr.Pz()
    
    br_b_px[0] = tlv_b.Px()
    br_b_py[0] = tlv_b.Py()
    br_b_pz[0] = tlv_b.Pz()
    br_b_m[0] = tlv_b.M()
    br_b_id[0] = pdg_b    
    
    br_x1[0] = itree.mcevt_pdf_x1[0]
    br_x2[0] = itree.mcevt_pdf_x2[0]
    br_fl1[0] = itree.mcevt_pdf_id1[0]
    br_fl2[0] = itree.mcevt_pdf_id2[0]
    br_qsqr[0] = itree.mcevt_pdf_scale[0]

    if options.dojets:
        br_njets[0] = 0
        br_jet_px.clear()
        br_jet_py.clear()
        br_jet_pz.clear()
        if matchb:
            br_jet_id.clear()
        if itree.jet_antiKt4Truth_n > 0:
            for ij in xrange( itree.jet_antiKt4Truth_n ):
                j = ROOT.TLorentzVector()
                j.SetPtEtaPhiE( itree.jet_antiKt4Truth_pt[ij],
                                itree.jet_antiKt4Truth_eta[ij],
                                itree.jet_antiKt4Truth_phi[ij],
                                itree.jet_antiKt4Truth_E[ij] )
                if j.Pt() > minjetpt and abs(j.Eta()) < maxjeteta and j.DeltaR( tlv_lm ) > 0.4 and j.DeltaR( tlv_lp ) > 0.4:
                    br_njets[0] = br_njets[0]+1 
                    br_jet_px.push_back( j.Px() )
                    br_jet_py.push_back( j.Py() )
                    br_jet_pz.push_back( j.Pz() )
                    if matchb:
                        if j.DeltaR( tlv_b ) < 0.4:
                            br_jet_id.push_back( 5 )
                        else:
                            br_jet_id.push_back( 0 )

    otree.Fill()

otree.Write()
ohndl.Close()

print '==> Done.'
sys.stdout.flush()
