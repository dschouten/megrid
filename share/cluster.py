#!/usr/bin/env python
#
#
# @file cluster.py
# @author doug schouten <doug dot schouten at triumf dot ca>
#
# loads an ME grid and performs clustering on the grid points
# to reduce the size of the grid, thererby reducing memory footprint
# and speeding up the characterization of test points
#
# the clustering is performed using Euclidian metric (only option for now),
# so given a large enough grid, then clustering from N -> M (or some
# reasonable number of clusters) should not reduce the precision of the
# grid too much
#
# note that clustering is very expensive, up to O(n^{3}) ... !
# however, due to the fact that the grids are i.i.d. random points
# in phase space, it should be possible to (within statistic uncertainty)
# perform clustering on subsets of the grids and combine results
# such that the global clustering is ~ equivalent to the sum of clustered subsets
#

import ROOT

import getopt
import megrid
import sys
import msg

from commands import getstatusoutput as shell

from megrid import megrid
from megrid import transform

def usage():
    print 'USAGE: cluster [-v,--verbose] --input=/path/to/file:tree --output=/path/to/file:tree'
    print '                              --nclusters=N --coords=A:B:C --min=N --max=N [--algorithm=KMEANS]'
    sys.exit( -1 )

def normalize( s ):
    charmap = { '+' : '_plus_', '-' : '_minus_', '*' : '_times_', '/' : '_div_' }
    for char in '~!@#$%^&()=;?><.,|][}{':
        charmap[char] = ''
    for k,v in charmap.items():
        s = s.replace( k, v )
    return s

try:
    opts, args = getopt.getopt( sys.argv[1:], 'vi:o:n:a:c:',
                                ['verbose',
                                 'input=',
                                 'output=',
                                 'nclusters=',
                                 'min=',
                                 'max=',
                                 'algorithm=',
                                 'coords='] )
except getopt.GetoptError, error:
    print str( error )
    usage()

verbose   = False
input     = None
output    = None
nclusters = None
coords    = None
alg       = None
method    = None
imin      = -1
imax      = -1

for o, a in opts:
    if o in ( '-v', '--verbose' ):
        verbose = True
    elif o in ( '-i', '--input' ):
        input = a
    elif o in ( '-o', '--output' ):
        output = a
    elif o in ( '-n', '--nclusters' ):
        nclusters = int( a )
    elif o in ( '--min' ):
        imin = int( a )
    elif o in ( '--max' ):
        imax = int( a )
    elif o in ( '-c', '--coords' ):
        coords = a
    elif o in ( '-a', '--algorithm' ):
        if a.lower() == 'kmeans':
            alg = megrid.KMEANS
            method = 'm'
        if a.lower() in [ 'tree', 'hierarchical' ]:
            alg = megrid.HIERARCHICAL
            method = 'a'
    else:
        usage()

if input==None or output==None or nclusters==None or coords==None or alg==None:
    usage()

log = msg.msglog( 'cluster', 'debug', useColor = True )

# -------=======------- -------=======------- -------=======-------

ifilen, itreen = input.split(':')
ofilen, otreen = output.split(':')

ihndl = ROOT.TFile.Open( ifilen, 'read' )
if ihndl == None or ihndl.IsZombie():
    log.error( 'input file [%s] not accessible'%( ifilen ) )

if ihndl.Get( itreen ) == None:
    log.error( 'input file does not contain [%s]'%( itreen ) )

if imax < 0:
    imax = ihndl.Get( itreen ).GetEntries()
if imin < 0:
    imin = 0

ihndl.Close()

grid = megrid( 'grid', ifilen, itreen, coords, '(Entry$ >= %d)&&(Entry$ < %d)'%( imin, imax ), abs(imax-imin) )
grid.cluster( alg, nclusters, method, 'e' )
grid.store( ofilen, otreen, normalize(coords) )

