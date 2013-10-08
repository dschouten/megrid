import ROOT

from ROOT import gSystem
from ROOT import gROOT

tlv = ROOT.TLorentzVector

gSystem.Load( 'libTMVA' )
gSystem.Load( 'libmegrid' )
gSystem.Load( 'libmegrid_dict' )

transform = ROOT.transform
megrid    = ROOT.megrid

class megrid( ROOT.megrid ):
    __valid_objects = [ 'leadj', 'subleadj', 'lp', 'lm', 'recoil', 'met' ]
    __valid_attrs   = [ 'X()', 'Y()', 'Z()', 'Px()', 'Py()', 'Pz()' ]
    
    def __init__( self, name, gridfile_pattern=None, tree_name=None, coord_definition=None, weights_expression='', max_size=-1, ndim=-1 ):
        if( gridfile_pattern != None ):
            if not ( tree_name != None and \
                     coord_definition != None ):
                raise RuntimeError
            ROOT.megrid.__init__( self, name, gridfile_pattern, tree_name, coord_definition, weights_expression, max_size )
        else:
            if ndim == 0:
                raise RuntimeError
            ROOT.megrid.__init__( self, name, ndim )
        self._exprs = []
        pass

    def __call__( self, leptons, jets, met, recoil ):
        raise NotImplementedError

    def set_expressions( self, exprs ):
        for expr in exprs:
            if sum( [ expr.count( obj ) for obj in self.__valid_objects ] ) == 0 or \
               sum( [ expr.endswith( attr ) for attr in self.__valid_attrs ] ) == 0:
                raise RuntimeError
        self._exprs = exprs 
        pass

class lvlv( megrid ):
    def __init__( self, name, gridfile_pattern=None, tree_name=None, coord_definition=None, weights_expression='', max_size=-1, ndim=-1 ):
        megrid.__init__( self, name, gridfile_pattern, tree_name, coord_definition, weights_expression, max_size, ndim )
        pass

    def _parse_leptons( self, leptons ):
        try:
            lp, lm = leptons
        except:
            raise RuntimeError
        if lp.charge / lm.charge != -1:
            raise RuntimeError
        if lp.charge < 0:
            lp, lm = lm, lp
        return lp, lm
    
    def __call__( self, leptons, jets, met, recoil ):
        lp, lm = self._parse_leptons( leptons )
        args = [ eval(expr) for expr in self._exprs ]
        return self.pdf( *args )

class lvlvj( lvlv ):
    def __init__( self, name, gridfile_pattern=None, tree_name=None, coord_definition=None, weights_expression='', max_size=-1, ndim=-1 ):
        lvlv.__init__( self, name, gridfile_pattern, tree_name, coord_definition, weights_expression, max_size, ndim )
        pass

    def __call__( self, leptons, jets, met, recoil ):
        lp, lm = self._parse_leptons( leptons )
        leadj = jets[0]
        args = [ eval(expr) for expr in self._exprs ]
        return self.pdf( *args )

class lvlvjj( lvlv ):
    def __init__( self, name, gridfile_pattern=None, tree_name=None, coord_definition=None, weights_expression='', max_size=-1, ndim=-1 ):
        lvlv.__init__( self, name, gridfile_pattern, tree_name, coord_definition, weights_expression, max_size, ndim )
        pass

    def __call__( self, leptons, jets, met, recoil ):
        lp, lm = self._parse_leptons( leptons )
        leadj = jets[0]
        subleadj = jets[1]
        args = [ eval(expr) for expr in self._exprs ]
        return self.pdf( *args )
