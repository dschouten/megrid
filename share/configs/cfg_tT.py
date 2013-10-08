
from megrid import megrid

from megrid import lvlvjj
from megrid import transform

coords  = 'lp_px|lp_py|lp_pz|'
coords += 'lm_px|lm_py|lm_pz|'
coords += 'jet_px[0]|jet_py[0]|jet_pz[0]|'
coords += 'jet_px[1]|jet_py[1]|jet_pz[1]'

import tempfile
f,n = tempfile.mkstemp( suffix='.C', prefix='helper', dir='/tmp', text=True )
f = file( n, 'write' )
f.write( '#include <iostream>\n'                                      )
f.write( '#include <cmath>\n'                                         )
f.write( '\n'                                                         )
f.write( 'using namespace std;\n'                                     )
f.write( '\n'                                                         )
f.write( 'double eta(double x, double y, double z)\n'                 )
f.write( '{\n'                                                        )
f.write( '    double theta=atan(sqrt(pow(x,2)+pow(y,2))/fabs(z));\n'  )
f.write( '    //std::cout << log(tan(theta/2)) << std::endl;\n'       )
f.write( '    return -log(tan(theta/2.));\n'                          )
f.write( '}\n'                                                        )
f.close()
ROOT.gROOT.ProcessLine( '.L %s+'%(n) )
del f,n

wtfun_bg = '(1-0.75)*(abs(eta(jet_px[0],jet_px[0],jet_pz[0]))<2.5)+(abs(eta(jet_px[0],jet_px[0],jet_pz[0]))>2.5)'
wtfun_gb = '(1-0.75)*(abs(eta(jet_px[1],jet_px[1],jet_pz[1]))<2.5)+(abs(eta(jet_px[1],jet_px[1],jet_pz[1]))>2.5)'

bb_grid = lvlvjj( 'tT_bb', '/global/dschouten/grids/outputs/merged/tT/pythia/tT_bg.root', 'grid', coords, '(%s)*(%s)'%(wtfun_bg,wtfun_gb), max_size=-1 )
bg_grid = lvlvjj( 'tT_bg', '/global/dschouten/grids/outputs/merged/tT/pythia/tT_bg.root', 'grid', coords, wtfun_bg, max_size=-1 )
gb_grid = lvlvjj( 'tT_gb', '/global/dschouten/grids/outputs/merged/tT/pythia/tT_gb.root', 'grid', coords, wtfun_gb, max_size=-1 )
gg_grid = lvlvjj( 'tT_gb', '/global/dschouten/grids/outputs/merged/tT/pythia/tT_gg.root', 'grid', coords, max_size=-1 )

mehndl = lvlvjj( 'tT', ndim=len(coords.split('|')) )
mehndl += bg_grid
mehndl += gb_grid
mehndl += gg_grid
mehndl += bb_grid

del bg_grid
del gb_grid
del gg_grid
del bb_grid

mehndl.set_expressions( [ 'lp.X()', 'lp.Y()', 'lp.Z()',
                          'lm.X()', 'lm.Y()', 'lm.Z()',
                          'leadj.X()', 'leadj.Y()', 'leadj.Z()',
                          'subleadj.X()', 'subleadj.Y()', 'subleadj.Z()' ] )
                        
ndim = mehndl.get_ndim()

mehndl.set_wc( 20 )
mehndl.set_maxk( 5000 ) 
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
selection = 'nJets_Pt25_MV1_75==0 && BDT_75_allBkgs > -0.75'
