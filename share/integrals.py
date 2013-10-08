from megrid import megrid
from megrid import transform

from pycuba import Vegas as int_vegas
from pycuba import Miser as int_miser

def global_int_function( ndim, x, ncomp, f, userdata, weight=None, iter=None ):
    

class integrand( megrid ):
    def __init__(  self, name, gridfile_pattern=None, tree_name=None, coord_definition=None, weights_expression='', max_size=-1, ndim=-1 ):
        megrid.__init__( self, name, gridfile_pattern, tree_name, coord_definition, weights_expression, max_size, ndim )
        self._int_definition = {}
        self._int_parameters = []
        pass
    def add_integration( self, ivar, limits=None, transfer_function=None ):
        if limits==None and transfer_function==None:
            raise ValueError
        if limits==None:
            self._int_definition[ivar] = { 'limits' : ( transfer_function.get_xlo,
                                                        transfer_function.get_xhi ),
                                           'weights' : transfer_function }
        else:
            a,b = limits
            self._int_definition[ivar] = { 'limits' : (a,b)
                                           'weights' : transfer_function }
    def evaluate( self, leptons, jets, met, recoil ):
        raise NotImplemented

    def __call__( ndim, x, ncomp, f, userdata, weight=None, iter):
        
# @TODO
