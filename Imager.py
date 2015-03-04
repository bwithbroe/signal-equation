from numpy import *
from numpy.linalg import *

class Imager:

    def __init__(self):
    
        # Typical values for scan constants
        # Panagiotaki et al., Neuroimage 2012
        self.little_delta = .003 # seconds
        self.big_delta = .02     # seconds   (time from excitation to refocusing; usually 1/2 T_E)
        self.tau = .02           # seconds   (1/2 T_E)
        self.gamma = 41.065e6    # Hz/T
        self.G = 0.00125         # T/mm
        # self.G = 0.0018
        
        # Gradient direction in an axon-aligned coordinate frame
        # (axon goes along the Z direction)
        self.g = array([0.0, 0.0, 1.0])  # unitless / normalized
        # self.g_axon = float(g[2])
        # self.g_xy = norm(g[0:2])
        
        # b/"diffusion weighting factor".  Typical values are 250 -- 6000 s/mm^2.
        self.b = self.gamma**2.0 * self.G**2.0 * self.little_delta**2.0 * (self.big_delta - self.little_delta / 3.0)
        print self.b
    
    '''
    def __init__(self, little_delta, big_delta, tau, gamma, G, g):
    
        # Typical values for scan constants
        # Panagiotaki et al., Neuroimage 2012
        self.little_delta = little_delta # seconds
        self.big_delta = big_delta     # seconds
        self.tau = tau           # seconds
        self.gamma = gamma    # Hz/T
        #G = 0.00125         # T/mm
        self.G = G
        
        # Gradient direction in an axon-aligned coordinate frame
        # (axon goes along the Z direction)
        self.g = g  # unitless / normalized
        # self.g_axon = float(g[2])
        # self.g_xy = norm(g[0:2])
        
        # b/"diffusion weighting factor".  Typical values are 250 -- 6000 s/mm^2.
        self.b = gamma**2.0 * G**2.0 * little_delta**2.0 * (big_delta - little_delta / 3.0)
        '''     