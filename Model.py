from numpy import *
from numpy.linalg import *

class Model:

    def __init__(self):
        
        # Diffusivity of the goop inside and outside the axons.
        self.D_parallel = 0.0007         # mm^2/s
        self.D_perpendicular = 0.0007    # mm^2/s
        
        # Unweighted signal
        self.S_unweighted = 290          # unitless
        
        # Axon direction
        self.axon_direction = array([0.0, 0.0, 1.0])
        
        self.k = 1.88
        self.squiggly_theta = 1.28e-4
    
    '''
    def __init__(self, S_unweighted, axon_direction, k = 1.88, squiggly_theta = 1.28e-4,
                                     D_parallel = 1e-3, D_perpendicular = 1e-3):
    
        # Diffusivity of the goop inside and outside the axons.
        self.D_parallel = D_parallel         # mm^2/s
        self.D_perpendicular = D_perpendicular    # mm^2/s
        
        # Unweighted signal
        self.S_unweighted = S_unweighted          # unitless
        
        # Axon direction
        self.axon_direction = axon_direction
        
        self.k = k
        self.squiggly_theta = squiggly_theta
        '''