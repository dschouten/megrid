#!/usr/bin/env python
#
# @file rotate.py
# @author doug schouten <doug dot schouten at triumf dot ca>
#
# objective: reduce the dimensionality of an ME grid by removing redundant
# azimuthal degree of freedom (although can in principle rotate vectors by arbitrary angle)
#
# note: independent method as used in agf::transformations library
#
# USAGE: reduce [-v,--verbose] -i,--input=/path/to/file:tree -o,--output=/path/to/file:tree
#        -c,--coordinates="x|y|z" -k,--keys="x|y|z" -r,--rotation="atan(x/y)"
#
# -c,--coordinates : expression (any valid TTreeFormula) that identifies the grid coordinates,
#                    separated by '|'
# -k,--keys : simple names that identify both (a) the mapping of coordinates to particle
#             momenta, and (b) the names of branches in the output file
#             NOTE that names are expected to be of the format "px" or "py" where 'p' is a
#                  particle name identifier and 'x', 'y' are the particle momenta
# -r,--rotation : expression that defines the angle with which to rotate about the z-axis in
#                 the counter-clockwise direction
# 

import ROOT

import getopt
import megrid
import sys
import msg
import re
import warnings

from array import array 
from commands import getstatusoutput as shell

warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )

def usage():
    print 'USAGE: reduce [-v,--verbose] -i,--input=/path/to/file:tree -o,--output=/path/to/file:tree'
    print '       -c,--coordinates="x:y:z" -k,--keys="x:y:z" -r,--rotation="atan(x/y)"'
    sys.exit( -1 )

# 
# normalizes a string
#
def normalize( s ):
    charmap = { '+' : '_plus_', '-' : '_minus_', '*' : '_times_', '/' : '_div_' }
    for char in '~!@#$%^&()=;?><.,|][}{':
        charmap[char] = ''
    for k,v in charmap.items():
        s = s.replace( k, v )
    return s
    
try:
    opts, args = getopt.getopt( sys.argv[1:], 'vi:o:c:k:r:',
                                ['verbose',
                                 'input=',
                                 'output=',
                                 'coords=',
                                 'keys=',
                                 'rotation='] )
except getopt.GetoptError, error:
    print str( error )
    usage()

verbose     = False
input       = None
output      = None
coord_exprs = None
coord_keys  = None
rot_expr    = None

for o, a in opts:
    if o in ( '-v', '--verbose' ):
        verbose = True
    elif o in ( '-i', '--input' ):
        input = a
    elif o in ( '-o', '--output' ):
        output = a
    elif o in ( '-c', '--coords' ):
        coord_exprs = a.split( '|' )
    elif o in ( '-k', '--keys' ):
        coord_keys = a.split( '|' )
    elif o in ( '-r', '--rotation' ):
        rot_expr = a
    else:
        usage()

if input==None or output==None or coord_exprs==None:
    usage()

log = msg.msglog( 'reduce', 'debug', useColor = True )
log.setPrintOnly( False )

# -------=======------- -------=======------- -------=======-------

if len( coord_keys ) != len( coord_exprs ):
    usage()

coord_pairs = sorted( zip( coord_keys, coord_exprs ) )
coord_keys, coord_exprs = zip( *coord_pairs )

particles = list( set( map( lambda item: item[0:-1], coord_keys ) ) )

ifilen, itreen = input.split(':')
ofilen, otreen = output.split(':')

ihndl = ROOT.TFile.Open( ifilen, 'read' )
if ihndl == None or ihndl.IsZombie():
    log.error( 'input file [%s] not accessible'%( ifilen ) )

if ihndl.Get( itreen ) == None:
    log.error( 'input file does not contain [%s]'%( itreen ) )

itree = ihndl.Get( itreen )

coord_hndls = { }
for index, item in enumerate( coord_keys ):
    key = item.lower()
    coord_hndls[key] = { 'expression' : ROOT.TTreeFormula( 'c_%d'%(index), coord_exprs[index], itree ),
                         'value' : array('d',[0]) }
    coord_hndls[key]['expression'].SetQuickLoad( True )
    
rot_hndl = ROOT.TTreeFormula( 'angle', rot_expr, itree )
rot_hndl.SetQuickLoad( True )

rot_val  = 0.

ohndl = ROOT.TFile.Open( ofilen, 'recreate' )
otree = ROOT.TTree( otreen, otreen )

for index, item in enumerate( coord_keys ):
    key = item.lower()
    otree.Branch( item, coord_hndls[key]['value'], item + '/D' )

nentries = itree.GetEntries()

for ientry in xrange( nentries ):
    itree.GetEntry( ientry )

    for key in coord_keys:
        coord_hndls[key]['value'][0] = coord_hndls[key]['expression'].EvalInstance(0)
        
    rot_val = rot_hndl.EvalInstance(0)

    for part in particles:
        x = coord_hndls['%sx'%(part)]['value'][0]
        y = coord_hndls['%sy'%(part)]['value'][0]

        tlv = ROOT.TLorentzVector()
        tlv.SetXYZM( x, y, 0, 0 )
        tlv.RotateZ( rot_val )
        
        coord_hndls['%sx'%(part)]['value'][0] = tlv.X()
        coord_hndls['%sy'%(part)]['value'][0] = tlv.Y()
    
    otree.Fill()

    if nentries < 10 or ientry % (nentries / 10) == 0:
        log.debug( 'processed %d grid points'%( ientry ) )
    
ihndl.Close()

otree.Write()
ohndl.Close()
