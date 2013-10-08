
#
# script to combine the VBF grids
# optionally, perform clustering on the grids 
#

import os
import sys
import ROOT
import megrid
import msg
import locale

from glob import glob
from commands import getstatusoutput as shell

log = msg.msglog( 'debug' )
log.setPrintOnly( False )

# --------------------------------------------------------------

input_path     = '/global/dschouten/grids/outputs/VBF-ii'
output_path    = '/global/dschouten/grids/outputs/merged/VBF'

do_cluster     = False
size_threshold = 5000000
cluster_brick  = 50000

coordinates    = [ 'lp_px','lp_py','lp_pz',
                   'lm_px','lm_py','lm_pz',
                   'jet_px[0]','jet_py[0]','jet_pz[0]',
                   'jet_px[1]','jet_py[1]','jet_pz[1]' ]

rotation       = None # 'atan(lm_px/lm_py)' ## use runtime transformation instead  

# --------------------------------------------------------------

def get_key( item ):
    kmap = { 'jet_px[0]':'jet1_px',
             'jet_py[0]':'jet1_py',
             'jet_pz[0]':'jet1_pz',
             'jet_px[1]':'jet2_px',
             'jet_py[1]':'jet2_py',
             'jet_pz[1]':'jet2_pz' }
    if item in kmap.keys():
        return kmap[item]
    return item

exists = os.path.exists

coord_defs = ''
coord_keys = ''
for item in coordinates:
    coord_defs += '%s|'%(item)
    coord_keys += '%s|'%(get_key(item))
coord_defs = coord_defs[0:-1]
coord_keys = coord_keys[0:-1]

#
# step 1: rebuild the inclusive grid
#

if not exists( '/tmp/VBF_inclusive.root' ):
    for iset in xrange( 0, 20 ):
        if len( glob( '%s/*.mc12_v2._%03d*.grid.root'%( input_path, iset ) ) ) <= 1:
            continue
        stat, out = shell( 'hadd -f /tmp/VBF_inclusive_%d.root %s/*.mc12_v2._%03d*.grid.root'%( iset, input_path, iset ) )
        if stat != 0:
            log.fatal( out )
        log.info( '/tmp/VBF_inclusive_%d.root'%( iset ) )
    stat, out = shell( 'hadd -f /tmp/VBF_inclusive.root /tmp/VBF_inclusive_*.root' )
    if stat != 0:
        log.fatal( out )
    stat, out = shell( 'split.py --input=/tmp/VBF_inclusive.root --tree=grid --output="m_jet_n>=2:/tmp/VBF_inclusive_filtered.root"' )
    if stat != 0:
        log.fatal( out )
    shell( 'mv /tmp/VBF_inclusive_filtered.root /tmp/VBF_inclusive.root' )

#
# step 2: cluster the individual samples
#
#
# @TODO complete the clustering steps ...

#
# step 3: reduce the grid dimensionality
#

if rotation != None:
    log.info( 'apply rotation to the grid' )
    
    stat, out = shell( 'rotate.py -i /tmp/VBF_inclusive.root:grid' +
                       ' -o /tmp/VBF_inclusive_rotated.root:grid' + 
                       ' --coords="%s" --keys="%s" --rotation="%s"'%( coord_defs, coord_keys, rotation ) )
    if stat != 0:
        log.fatal( out )
    log.info( out )
    stat, out = shell( 'mv /tmp/VBF_inclusive_rotated.root /tmp/VBF_inclusive.root' )

#
# step 4: store in output location
#

stat, out = shell( 'cp /tmp/VBF_inclusive.root %s/.'%( output_path ) )
stat, out = shell( 'rm -rf /tmp/VBF_inclusive.root /tmp/VBF_inclusive_*.root' )
