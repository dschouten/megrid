
from megrid import megrid

from megrid import lvlvjj
from megrid import transform

coords  = 'lp_px|lp_py|lp_pz|'
coords += 'lm_px|lm_py|lm_pz|'
coords += 'jet_px[0]|jet_py[0]|jet_pz[0]|'
coords += 'jet_px[1]|jet_py[1]|jet_pz[1]'

bg_grid = lvlvjj( 'Wt_bg', '/global/dschouten/grids/outputs/merged/Wt/Wt_bg.root', 'grid', coords )
gb_grid = lvlvjj( 'Wt_gb', '/global/dschouten/grids/outputs/merged/Wt/Wt_gb.root', 'grid', coords )

mehndl = lvlvjj( 'Wt', ndim=len(coords.split('|')) )
mehndl += bg_grid
mehndl += gb_grid

del bg_grid
del gb_grid

mehndl.set_expressions( [ 'lp.X()', 'lp.Y()', 'lp.Z()',
                          'lm.X()', 'lm.Y()', 'lm.Z()',
                          'leadj.X()', 'leadj.Y()', 'leadj.Z()',
                          'subleadj.X()', 'subleadj.Y()', 'subleadj.Z()' ] )
                        
ndim = mehndl.get_ndim()

mehndl.set_wc( 50 )
mehndl.set_maxk( 2500 )
mehndl.set_option_pdf( megrid.AGF )

# scale and normalize the grid axes
scale_tfm = transform.scale_fixed( array('f', [1.0/GeV for idim in xrange( ndim )]), ndim )
mehndl.apply_transformation( scale_tfm, False ) # convert grid to units of GeV

norm_tfm = transform.norm_gauss( ndim  )
mehndl.apply_transformation( norm_tfm )

# rotate away one of the azimuthal degrees of freedom
rotate_tfm = transform.rotate_phi( 'X0:Y0:Z0:X1:Y1:Z1:X2:Y2:Z2:X3:Y3:Z3' )
mehndl.apply_transformation( rotate_tfm )

# no recoil calculation
boostCalculator = None

# apply selection cuts
selection = 'nJets_Pt25_MV1_75==0 && BDT_75_allBkgs > -0.5'
