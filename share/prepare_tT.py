
#
# script to combine the tT grids
# and then split by jet flavour topology:
#    - both leading jets are b's
#    - leading is b, subleading is light flavour/gluon
#    - leading is light flavour/gluon, subleading is b
#
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

input_path     = '/global/dschouten/grids/outputs/tT'
output_path    = '/global/dschouten/grids/outputs/merged/tT/pythia'

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

if not exists( '/tmp/tT_inclusive.root' ):
    for iset in xrange( 0, 20 ):
        if len( glob( '%s/*.mc12_v*._%03d*.grid.root'%( input_path, iset ) ) ) <= 1:
            continue
        stat, out = shell( 'hadd -f /tmp/tT_inclusive_%d.root %s/*.mc12_v*._%03d*.grid.root'%( iset, input_path, iset ) )
        if stat != 0:
            log.fatal( out )
        log.info( '/tmp/tT_inclusive_%d.root'%( iset ) )
    stat, out = shell( 'hadd -f /tmp/tT_inclusive.root /tmp/tT_inclusive_*.root' )
    if stat != 0:
        log.fatal( out )

#
# step 2: split the inclusive sample by jet flavour topology
#

if not exists( '/tmp/tT_bb.root' ) or \
   not exists( '/tmp/tT_bg.root' ) or \
   not exists( '/tmp/tT_gb.root' ) or \
   not exists( '/tmp/tT_gg.root' ):
    log.info( 'split by jet flavour topology' )
    stat, out = shell( 'split.py --input=/tmp/tT_inclusive.root --tree=grid' +
                       ' --output="m_jet_n>=2 && jet_id[0]==5 && jet_id[1]==5:/tmp/tT_bb.root"' +
                       ' --output="m_jet_n>=2 && jet_id[0]==5 && jet_id[1]!=5:/tmp/tT_bg.root"' +
                       ' --output="m_jet_n>=2 && jet_id[0]!=5 && jet_id[1]==5:/tmp/tT_gb.root"' +
                       ' --output="m_jet_n>=2 && jet_id[0]!=5 && jet_id[1]!=5:/tmp/tT_gg.root"' )
    if stat != 0:
        log.fatal( out )
    log.info( out )
    
#
# step 3: cluster the individual samples
#
#
# @TODO complete the clustering steps ...

#
# step 4: reduce the grid dimensionality
#

if rotation != None:
    for jetcomb in [ 'bb', 'bg', 'gb', 'gg' ]:
        stat, out = shell( 'rotate.py -i /tmp/tT_%s.root:grid'%(jetcomb) +
                           ' -o /tmp/tT_%s_rotated.root:grid'%(jetcomb) + 
                           ' --coords="%s" --keys="%s" --rotation="%s"'%( coord_defs, coord_keys, rotation ) )
        if stat != 0:
            log.fatal( out )
        log.info( out )
        stat, out = shell( 'mv /tmp/tT_%s_rotated.root /tmp/tT_%s.root'%(jetcomb, jetcomb) )

#
# step 5: store in output location
#

for jetcomb in [ 'bb', 'bg', 'gb', 'gg' ]:
    stat, out = shell( 'cp /tmp/tT_%s.root %s/.'%(jetcomb, output_path) )
    stat, out = shell( 'rm -rf /tmp/tT_%s.root'%(jetcomb) )
stat, out = shell( 'rm -rf /tmp/tT_inclusive_%d.root' )
stat, out = shell( 'rm -rf /tmp/tT_inclusive.root' )
