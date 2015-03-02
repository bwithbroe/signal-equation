from numpy import *
from numpy.linalg import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D

# Typical values for scan constants
# Panagiotaki et al., Neuroimage 2012
little_delta = .003 # seconds
big_delta = .02     # seconds
tau = .02           # seconds
gamma = 41.065e6    # Hz/T
#G = 0.00125         # T/mm
G = 0.0018

# Gradient direction in an axon-aligned coordinate frame
# (axon goes along the Z direction)
g = matrix('.5; .1; .86023')  # unitless / normalized
g_z = float(g[2])
g_xy = norm(g[0:2])

# Diffusivity of the goop inside and outside the axons.
D_parallel = 0.0007         # mm^2/s
D_perpendicular = 0.0007    # mm^2/s

# Unweighted signal
S_unweighted = 290          # unitless

# ================= DERIVED VALUES =================

# b/"diffusion weighting factor".  Typical values are 250 -- 6000 s/mm^2.
b = gamma**2.0 * G**2.0 * little_delta**2.0 * (big_delta - little_delta / 3.0)

# Components of the q-vector in the axon's coordinate frame.
q_parallel = (gamma * G * g_z * little_delta) / (2.0 * pi)
q_perpendicular = (gamma * G * g_xy * little_delta) / (2.0 * pi)

print "Effective b-value:", b
print "S-T exponent:     ", -b * D_parallel
print "Expected signal:  ", (S_unweighted * exp(-b * D_parallel))


def S_cyl(radius):
    # Assaf et al., MRM 2004
    log_S_restr_plel = -4.0 * pi**2.0 * q_parallel**2.0 * (big_delta - (little_delta/3.0)) * D_parallel
    # BUG --- this should reduce to the previous equation as radius -> inf
    log_S_restr_perp = (-4.0 * pi**2.0 * radius**4.0 * q_perpendicular**2.0 * D_perpendicular * tau *
                        (7.0/96.0) * (2.0-(99.0/112.0)*(radius**2.0/(D_perpendicular*tau))))
    # print "rest_plel exponent:", log_S_restr_plel
    # print "rest_perp exponent:", log_S_restr_perp
    return exp(log_S_restr_plel + log_S_restr_perp)

def Gamma(radius, k, sqiggly_theta):
    return ( (radius**(k-1.0) * exp(-(radius/squiggly_theta))) /
             (squiggly_theta**k * float(math.factorial(k))) )

S_in = 0.0
Gamma_sum = 0.0
k = 2
squiggly_theta = 1.28e-4

S_out = S_unweighted * exp(-b * D_parallel)
charmed_S = S_unweighted * S_cyl(1e-3)
print "Predicted signal with axon radius 1e-3 mm", charmed_S


width = 5.0e-6
for i in range(0, 100):
    middle_radius = width/2 + i * width
    #print "S_Cyl at", middle_radius, " is:", S_unweighted * S_cyl(middle_radius) * 5.0e-5
    S_in_middle = S_cyl(middle_radius) * Gamma(middle_radius, k, squiggly_theta)
    S_in += S_in_middle * width
    Gamma_sum += Gamma(middle_radius, k, squiggly_theta) * width

print ""
print "Total signal over range of radii:", S_in + S_out
print "Gamma sum", Gamma_sum
