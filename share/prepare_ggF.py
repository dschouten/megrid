
#
# script to combine the ggF grids
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

input_path     = '/global/dschouten/grids/outputs/ggF-ii'
output_path    = '/global/dschouten/grids/outputs/merged/ggF'

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

coord_defs_0j = ''
coord_keys_0j = ''
coord_defs_1j = ''
coord_keys_1j = ''
coord_defs_2j = ''
coord_keys_2j = ''
for item in coordinates:
    if 'jet' not in item:
        coord_defs_0j += '%s|'%(item)
        coord_keys_0j += '%s|'%(get_key(item))
    if 'jet' in item and item.endswith('[0]'):
        coord_defs_1j += '%s|'%(item)
        coord_keys_1j += '%s|'%(get_key(item))
    if 'jet' in item and (item.endswith('[0]') or
                          item.endswith('[1]')):
        coord_defs_2j += '%s|'%(item)
        coord_keys_2j += '%s|'%(get_key(item))
coord_defs_0j = coord_defs_0j[0:-1]
coord_keys_0j = coord_keys_0j[0:-1]
coord_defs_1j = coord_defs_1j[0:-1]
coord_keys_1j = coord_keys_1j[0:-1]
coord_defs_2j = coord_defs_2j[0:-1]
coord_keys_2j = coord_keys_2j[0:-1]

#
# step 1: rebuild the inclusive grid
#

if not exists( '/tmp/ggF_inclusive.root' ):
    for iset in xrange( 0, 20 ):
        if len( glob( '%s/*.mc12_v3._%03d*.grid.root'%( input_path, iset ) ) ) <= 1:
            continue
        stat, out = shell( 'hadd -f /tmp/ggF_inclusive_%d.root %s/*.mc12_v3._%03d*.grid.root'%( iset, input_path, iset ) )
        if stat != 0:
            log.fatal( out )
        log.info( '/tmp/ggF_inclusive_%d.root'%( iset ) )
    stat, out = shell( 'hadd -f /tmp/ggF_inclusive.root /tmp/ggF_inclusive_*.root' )
    if stat != 0:
        log.fatal( out )
    stat, out = shell( 'split.py --input=/tmp/ggF_inclusive.root --tree=grid --output="m_jet_n==0:/tmp/ggF_inclusive_0j.root"' )
    if stat != 0:
        log.fatal( out )
    stat, out = shell( 'split.py --input=/tmp/ggF_inclusive.root --tree=grid --output="m_jet_n==1:/tmp/ggF_inclusive_1j.root"' )
    if stat != 0:
        log.fatal( out )
    stat, out = shell( 'split.py --input=/tmp/ggF_inclusive.root --tree=grid --output="m_jet_n>=2:/tmp/ggF_inclusive_2j.root"' )
    if stat != 0:
        log.fatal( out )
        
#
# step 2: cluster the individual samples
#
#
# @TODO complete the clustering steps ...


#
# step 3: reduce the grid dimensionality
#

if rotation != None:
    for njets in [ '0j', '1j', '2j' ]:
        defs = vars()['coord_defs_%s'%(njets)]
        keys = vars()['coord_keys_%s'%(njets)]
        
        log.info( 'apply rotation to the grid' )
        
        stat, out = shell( 'rotate.py -i /tmp/ggF_inclusive_%s.root:grid'%(njets) +
                           ' -o /tmp/ggF_%s_rotated.root:grid'%(njets) + 
                           ' --coords="%s" --keys="%s" --rotation="%s"'%( defs, keys, rotation ) )
        if stat != 0:
            log.fatal( out )
        log.info( out )
        stat, out = shell( 'mv /tmp/ggF_%s_rotated.root /tmp/ggF_%s.root'%(njets,njets) )

#
# step 4: store in output location
#

for njets in [ '0j', '1j', '2j' ]:
    stat, out = shell( 'cp /tmp/ggF_%s.root %s/.'%(njets, output_path) )
    stat, out = shell( 'rm -rf /tmp/ggF_%s.root'%(njets) )
stat, out = shell( 'rm -rf /tmp/ggF_inclusive_*.root' )
