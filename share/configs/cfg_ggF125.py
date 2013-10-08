
from megrid import megrid

from megrid import lvlv
from megrid import transform

coords  = 'lp_px|lp_py|lp_pz|'
coords += 'lm_px|lm_py|lm_pz|'
coords += '(wp_px+wm_px)|(wp_py+wm_py)'

mehndl = lvlv( 'ggF125', '/global/dschouten/grids/outputs/merged/ggF/ggF_inclusive.root', 'grid', coords, max_size=-1 )

mehndl.set_expressions( [ 'lp.X()', 'lp.Y()', 'lp.Z()',
                          'lm.X()', 'lm.Y()', 'lm.Z()',
                          'recoil.X()', 'recoil.Y()' ] )
                        
ndim = mehndl.get_ndim()

mehndl.set_wc( 8 )
mehndl.set_maxk( 10000 ) 
mehndl.set_option_pdf( megrid.AGF )

# scale and normalize the grid axes
scale_tfm = transform.scale_fixed( array('f', [1.0/GeV for idim in xrange( ndim )]), ndim )
mehndl.apply_transformation( scale_tfm, False ) # convert grid to units of GeV

# norm_tfm = transform.norm_gauss( ndim  )
# mehndl.apply_transformation( norm_tfm )

# rotate away one of the azimuthal degrees of freedom
# rotate_tfm = transform.rotate_phi( 'X0:Y0:Z0:X1:Y1:Z1:X2:Y2' )
# mehndl.apply_transformation( rotate_tfm )

# configure recoil estimate
boostCalculator = 'configs/tmvaRecoilEstimator.py'

# apply selection cuts
selection = 'Mll/1000>10 && abs(lepID0)!=abs(lepID1) && Ptll/1000>10 && METRel/1000>25 && Mll/1000<100'
selection = 'METRel/1000>25 && Ptll/1000>20 && Mll/1000<80 && DPhill<2.8 && DPhillMET>%04g'%( math.pi / 2  )
