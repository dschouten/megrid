#!/usr/bin/python

#
# @file calculateME.py
#
# Input:
#  cfg=   : python module(s) that configure the calculation
#  input= : ntuple file with the events to loop over
#  ouput= : the file to store the calculations in
#  

import sys
import os
import random
import math
import array
import getopt
import glob
import signal
import time
import pickle
import string

from commands import getstatusoutput as shell
from glob import glob
from array import array

#_________________________________________________________________

append      = False
update      = False
verbose     = False
ifilen      = None
ofilen      = None
begevent    = 0
endevent    = 0
randomsel   = 0
randomseed  = 0
cfgmodules  = []
cfgexec     = ''
itreen      = 'HWWTree'

def usage():
    print 'USAGE: [-v] [-a] [-u] --cfg=cfg.py [--exec="mehndl.setX()"] --input=ww.root --output=me_ww.root'
    print '       [--begin=0] [--end=0] [--random=0] [--seed=0]'
    sys.exit( -1 )

try:
    opts, args = getopt.getopt( sys.argv[1:], 'vat:c:x:i:o:b:e:r:s:q:',
                                ['timeout=',
                                 'cfg=',
                                 'exec=',
                                 'input=',
                                 'output=',
                                 'begin=',
                                 'end=',
                                 'random=',
                                 'seed=',
                                 'tree='] )
except getopt.GetoptError, error:
    print str( error )
    usage()

for o, a in opts:
    if o in ( '-v' ):
        verbose = True
    if o in ( '-a' ):
        append = True
    if o in ( '-u' ):
        update = True
    if o in ( '-o', '--output' ):
        ofilen = a
    if o in ( '-i', '--input' ):
        ifilen = sorted( a.split( ',' ) )
    if o in ( '-c', '--cfg' ):
        cfgmodules += [ a ]
    if o in ( '-x', '--exec' ):
        cfgexec = a
    if o in ( '-b', '--begin' ):
        begevent = int( a )
    if o in ( '-e', '--end' ):
        endevent = int( a )
    if o in ( '-r', '--random' ):
        randomsel = int( a )
    if o in ( '-s', '--seed' ):
        randomseed = int( a )
    if o in ( '-q', '--tree' ):
        itreen = a

if len(cfgmodules)==0 or ofilen==None or ifilen==None:
    usage()    

#_________________________________________________________________

from ROOT import gROOT
gROOT.SetBatch( True )

import ROOT
import matrixelement
import msg

#_________________________________________________________________
class Electron( ROOT.TLorentzVector ):
    def __init__( self, tlv, charge ):
        ROOT.TLorentzVector.__init__( self, tlv.X(), tlv.Y(), tlv.Z(), tlv.E() )
        self.charge = charge
        self.type = 'electron'
        
#_________________________________________________________________
class Muon( ROOT.TLorentzVector ):
    def __init__( self, tlv, charge ):
        ROOT.TLorentzVector.__init__( self, tlv.X(), tlv.Y(), tlv.Z(), tlv.E() )
        self.charge = charge
        self.type = 'muon'

#_________________________________________________________________
def fbranch():
    return array('d', [0])

def ibranch():
    return array('i', [0])

def farrbranch():
    return ROOT.std.vector('double')()

def normalize_string( s ):
    for char in ',./<>?;\':"[]\{}\!@#$%^&*()+-=':
        s = s.replace( char, '' )
    return s

def buffer_result( r, t=None ):
    rfile = open( 'result_buffer.pickle', 'w' )
    pickle.dump( r, rfile )
    rfile.close()
    if t != None:
        t.Fill()

def id_generator( size=6, chars=string.ascii_uppercase + string.digits ):
    return ''.join( random.choice(chars) for x in range(size) )

#_________________________________________________________________

#
# define some useful objects / shortcuts / globals
#

pdg = ROOT.TDatabasePDG( )

lepMasses = { 13 : pdg.GetParticle( 'mu-' ).Mass( ),
              11 : pdg.GetParticle( 'e-' ).Mass( ) }
lepTypes  = { 13 : Muon,
              11 : Electron }

tlv = ROOT.TLorentzVector

GeV = 1000.

log = None
if verbose: log = msg.msglog( 'calculateME', 'debug', useColor = False )
else:       log = msg.msglog( 'calculateME', 'info',  useColor = False )

ROOT.gErrorIgnoreLevel = ROOT.kError

timer = ROOT.TStopwatch( )

doBuffering = False

me_value  = fbranch()
me_error  = fbranch()
me_neval  = fbranch()
me_time   = fbranch()
me_entry  = ibranch()
me_evtnum = ibranch()
me_runnum = ibranch()
me_extra  = farrbranch()

#
# quick sanity checks
#

if update and len( ifilen )>1:
    log.fatal( 'cannot run in update mode with more than one input file' )

#
# the following variables steer the calculation
# and can be set in the configuration modules
#

selection            = None   # selection cuts to apply to events before calculating the ME
boostCalculator      = None   # name of python file containing code to estimate system recoil for each event
    
log.info( 'start initializing ME grid(s)' )
timer.Start()
for module in cfgmodules:
    if os.path.exists( module ):
        execfile( module )
    else:
        print 'ERROR: configuration [%s] does not exist'%( module )
        sys.exit( -1 )

if len( cfgexec ) > 0:
    print
    print
    print cfgexec
    print
    print
    exec( cfgexec )
timer.Stop()
title = normalize_string( mehndl.get_name() )

log.info( 'finished initializing ME grid(s) in [%d] seconds'%( timer.CpuTime() ) )

stat, host = shell( 'hostname' )
stat, time = shell( 'date' )
log.info( 'running on [%s] @ %s'%( host, time ) )
log.info( 'will perform [%s] ME calculation on file(s) %s'%( title, str(ifilen)  ) )
    
# ------------------------------------------------------------------------
#
# merge and load the input file(s)
#
# ------------------------------------------------------------------------

hadded = False
if len( ifilen ) > 1:
    ifilen = '%s.root'%( id_generator() )
    stat, out = shell( ' '.join( ['hadd -f %s.root'%( )] + ifilen ) )
    log.debug( out )
    hadded = True
    log.info( 'merged files into single input [%s]'%( ifilen ) )
else:
    ifilen = ifilen.pop()

mode = 'read'
if update:
    mode = 'update'
    log.warning( 'updating original file' )    
ifile = ROOT.TFile.Open( ifilen, mode )
if ifile == None or ifile.IsZombie():
    log.fatal( 'file [%s] not found or cannot be opened'%( ifilen ) )
    
itree = ifile.Get( itreen )
itree.SetBranchStatus( '*', True )

cyclemap = {}
keyobjs  = []
keynames = []

for item in ifile.GetListOfKeys():
    keyobjs += [item]

for (key, item) in zip( map( lambda k: k.GetName(), keyobjs ), keyobjs ):
    keynames += key
    if key not in cyclemap:
        cyclemap[key] = []
    cyclemap[key] += [ item.GetCycle() ]
        
for item in keyobjs:
    if item.ReadObj().Class() == 'TTree' and \
           item.GetCycle() == max( cyclemap[item.GetName()] ):
        itree.AddFriend( item.GetName() )

if update:
    itree.Branch( 'P_%s'%(title), me_value, 'P_%s/D'%(title) )

# ------------------------------------------------------------------------
#
# load the output file
#
# ------------------------------------------------------------------------

ofile = ROOT.TFile.Open( ofilen, 'recreate' )
if ofile == None or ofile.IsZombie():
    log.fatal( 'file [%s] can not be accessed'%( ofilen ) )

otreename = title
otreedesc = title

spectatortrees = []

if append:
    keynames = []       
    for item in keyobjs:
        if item.GetCycle() != max( cyclemap[item.GetName()] ):
            continue
        item = item.ReadObj()
        if item.Class().InheritsFrom( ROOT.TTree.Class() ):
            tree = ifile.Get( item.GetName() )
            tree.SetBranchStatus( '*', True )
            ofile.cd( )
            clone = tree.CloneTree( 0 )
            spectatortrees.append( (tree,clone) )
        if item.Class().InheritsFrom( ROOT.TH1.Class() ):
            ofile.cd( )
            buf = item.Clone( item.GetName() )
            buf.Write( )
    otreename = title
    inc = 0
    while otreename in keynames: 
        otreename = title + '_v%02d'%( inc )
        inc += 1

otree = ROOT.TTree( otreename, otreedesc )

me_value_br  = otree.Branch( 'me'     , me_value , 'me/D' )
me_error_br  = otree.Branch( 'error'  , me_error , 'error/D' )
me_neval_br  = otree.Branch( 'neval'  , me_neval , 'neval/D' )
me_time_br   = otree.Branch( 'time'   , me_time  , 'time/D' )
me_entry_br  = otree.Branch( 'entry'  , me_entry , 'entry/I' )
me_evtnum_br = otree.Branch( 'evtnum' , me_evtnum, 'evtnum/I' )
me_runnum_br = otree.Branch( 'runnum' , me_runnum, 'runnum/I' )
me_extra_br  = otree.Branch( 'extra'  , me_extra )

# ------------------------------------------------------------------------
#
# loop over the events
#
# ------------------------------------------------------------------------

a = min( begevent, itree.GetEntries() )
b = min( endevent, itree.GetEntries() )
if b <= 0:
    b = itree.GetEntries()

entries = range( a, b )

elist = None
if selection != None and len(selection) > 0:
    itree.Draw( ">>entrylist", selection )
    elist = ROOT.gROOT.FindObject( "entrylist" )
    if elist != None:
        log.info( 'using [%d] selected events of [%d] total'%( elist.GetN(), itree.GetEntries() ) )
    else:
        log.warning( 'failed to apply selection [%s]'%( selection ) )

random.seed( randomseed ) 

if randomsel > 0:
    entries = [ elist.GetEntry(index) for index in xrange( elist.GetN() ) ]
    entries = filter( lambda ientry: a<ientry and ientry<=b, entries )
    entries = random.sample( entries, min( len(entries), randomsel ) )

nprocessed = 0

for ientry in entries:
    measured_lp       = tlv() 
    measured_lm       = tlv() 
    measured_jets     = []
    measured_met      = tlv() 
    measured_recoil   = tlv()

    me_extra.clear()

    # ------------------------------------------------------------------------
    # initialize the results table 
    # ------------------------------------------------------------------------
    
    result = { 'me' : -1, 'ievent' : -1, 'irun' : -1, 'ientry' : ientry, 'time' : -1 }
    if doBuffering: buffer_result( result )

    # ------------------------------------------------------------------------
    # selection applied here; use TEntryList with proper TCuts definition
    # ------------------------------------------------------------------------
       
    if ( elist != None and (not elist.Contains( ientry )) ):
        continue
    
    nb = itree.GetEntry( ientry )
    if nb <= 0:
        log.error( 'reading event #%d'%( ientry ) )
        sys.exit( -1 )
    
    log.debug( 'reading event #%d'%( ientry ) )

    ievent = getattr( itree, 'EventNumber' ) 
    irun   = getattr( itree, 'RunNumber' )

    log.info( "event number: ", ievent )
    
    if hasattr( itree, 'mc_channel_number' ):
        irun = getattr( itree, 'mc_channel_number' )
        
    result = { 'me' : -1, 'ievent' : ievent, 'irun' : irun, 'ientry' : ientry, 'time' : -1 }

    me_value[0]  = result['me']
    me_entry[0]  = result['ientry']
    me_time[0]   = result['time']
    me_evtnum[0] = result['ievent']
    me_runnum[0] = result['irun']
     
    # ------------------------------------------------------------------------
    # buffer the results table (if integration times out or fails, sensible
    # entries will be put in the output buffer)
    # ------------------------------------------------------------------------
    
    if doBuffering: buffer_result( result )
    
    # ------------------------------------------------------------------------
    # get lepton flavours and charges
    # ------------------------------------------------------------------------

    lepID0 = itree.lepID0 
    lepID1 = itree.lepID1

    random.seed( ) 

    if lepID0 * lepID1 > 0: ## deal with same-sign charged leptons by randomizing the charge assignment
        s = -2 * (random.uniform(0,1) >= 0.5) + 1
        lepID0 = ( s ) * abs(lepID0)
        lepID1 = (-s ) * abs(lepID1)
    
    if lepID0 > 0: ## lepID contains charge & PDG code information
        measured_lp.SetPtEtaPhiM( itree.lepPt0 / GeV, itree.lepEta0, itree.lepPhi0, lepMasses[abs(lepID0)] )
        measured_lm.SetPtEtaPhiM( itree.lepPt1 / GeV, itree.lepEta1, itree.lepPhi1, lepMasses[abs(lepID1)] )
        measured_lp = lepTypes[abs(lepID0)]( measured_lp,  1)
        measured_lm = lepTypes[abs(lepID1)]( measured_lm, -1)
    else: 
        measured_lm.SetPtEtaPhiM( itree.lepPt0 / GeV, itree.lepEta0, itree.lepPhi0, lepMasses[abs(lepID0)] )
        measured_lp.SetPtEtaPhiM( itree.lepPt1 / GeV, itree.lepEta1, itree.lepPhi1, lepMasses[abs(lepID1)] )
        measured_lm = lepTypes[abs(lepID0)]( measured_lm, -1)
        measured_lp = lepTypes[abs(lepID1)]( measured_lp,  1)

    log.info( '\t lp (pT, eta, phi): %0.2f, %0.2f, %0.2f'%( measured_lp.Pt(), measured_lp.Eta(), measured_lp.Phi() ) )
    log.info( '\t lm (pT, eta, phi): %0.2f, %0.2f, %0.2f'%( measured_lm.Pt(), measured_lm.Eta(), measured_lm.Phi() ) )
    log.info( '\t ll (pT): %0.2f'%( (measured_lp + measured_lm).Pt() ) )

    # ------------------------------------------------------------------------
    # get the jets for the ME
    # ------------------------------------------------------------------------
       
    for ijet in xrange( itree.m_jet_n ):
        measured_jets.append( tlv() )
        thejet = measured_jets[-1]
        try:
            thejet.SetPtEtaPhiM( itree.m_jet_pt[ijet] / GeV, itree.m_jet_eta[ijet], itree.m_jet_phi[ijet], 0. )
        except AttributeError:
            if ijet < 2:
                thejet.SetPtEtaPhiM( getattr( itree, 'jetPt%d'%( ijet ) ) / GeV,
                                     getattr( itree, 'jetEta%d'%( ijet ) ),
                                     getattr( itree, 'jetPhi%d'%( ijet ) ), 0 )
            else:
                continue
        log.info( '\t jet #%d (pT, eta, phi): %0.2f, %0.2f, %0.2f'%( ijet+1, thejet.Pt(), thejet.Eta(), thejet.Phi() ) )
    
    # ------------------------------------------------------------------------
    # get the MET & recoil estimate
    # ------------------------------------------------------------------------
    
    measured_met.SetPxPyPzE( itree.MET_x / GeV, itree.MET_y / GeV, 0., itree.MET / GeV )
    log.info( '\t MET (pT, phi): %0.2f, %0.2f'%( measured_met.Pt(), measured_met.Phi() ) )

    if boostCalculator != None:
        execfile( '%s'%( boostCalculator ) ) # << put the boost snippet in another file 
        log.info( '\t recoil (pT, phi): %0.2f, %0.2f'%( measured_recoil.Pt(), measured_recoil.Phi() ) )       
        me_extra.push_back( measured_recoil.Px() )
        me_extra.push_back( measured_recoil.Py() )
    
    buffer_result( result )
    
    # ------------------------------------------------------------------------
    # start integrating the ME calculation
    # ------------------------------------------------------------------------
    
    timer.Start( )

    result['me'] = mehndl( [measured_lp, measured_lm], measured_jets, measured_met, measured_recoil )
    buffer_result( result )
    
    timer.Stop( )
    
    result['time'] = timer.CpuTime()
    
    log.info( 'ME result:', result )
    
    me_value[0]  = result['me']
    me_entry[0]  = result['ientry']
    me_time[0]   = result['time']
    me_evtnum[0] = result['ievent']
    me_runnum[0] = result['irun']

    log.info( 'completed event #%d in %0.2f (s)'%( nprocessed + 1, me_time[0] ) )

    sys.stdout.flush()

    if append:
        for pair in spectatortrees:
            orig, clone = pair
            orig.GetEntry( ientry )
            if orig.Class().GetName() == 'TNtuple':
                clone.Fill( array( 'f', ( getattr( orig, br.GetName() ) for br in orig.GetListOfBranches() ) ) )
            else:
                clone.Fill( )
    
    otree.Fill( )
    if update:
        itree.Fill( )
    
    buffer_result( result )

    nprocessed += 1
    
if append:
    for pair in spectatortrees:
        orig, clone = pair
        clone.Write( )

otree.Write( )
ofile.Close( )

if update:
    itree.Write( itree.GetName(), ROOT.TObject.kOverwrite )
ifile.Close( )

if hadded:
    shell( 'rm -rf %s'%( ifilen ) )
