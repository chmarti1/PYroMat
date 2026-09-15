import pyromat as pm
import numpy as np



class Cycle:
    """Cycle prototype class
    c = Cycle(substance)
    
DO NOT USE THIS CLASS.

This class only exists to provide common methods to all cycle classes.
"""
    def __init__(self, substance):
        # The PYroMat substance
        self.subst = substance
        

class BraytonC(Cycle):
    pass
    
class RankineC(Cycle):
    pass
    
class OttoC(Cycle):
    pass
    
class AtkinsonC(Cycle):
    pass
    
class DieselC(Cycle):
    pass
    
class RefrigerationC(Cycle):
    pass
