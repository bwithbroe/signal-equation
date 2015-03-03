from numpy import *
from numpy.linalg import *
from scipy.special import gamma as gamma_func
from Model import *
from Imager import *

print "\n==========================="
imager = Imager()
model = Model()

def S_cyl(radius):
    # Assaf et al., MRM 2004
    log_S_restr_plel = -4.0 * pi**2.0 * q_parallel**2.0 * (imager.big_delta - (imager.little_delta/3.0)) * model.D_parallel
    # BUG --- this should reduce to the previous equation as radius -> inf
    log_S_restr_perp = (-((4.0 * pi**2.0 * radius**4.0 * q_perpendicular**2.0) / (model.D_perpendicular * imager.tau)) *
                        (7.0/96.0) * (2.0-(99.0/112.0)*(radius**2.0/(model.D_perpendicular * imager.tau))))
    # print "rest_plel exponent:", log_S_restr_plel
    # print "rest_perp exponent:", log_S_restr_perp
    return exp(log_S_restr_plel + log_S_restr_perp)

def Gamma(radius, k, squiggly_theta):
    return ( (radius**(k-1.0) * exp(-(radius/squiggly_theta))) /
             (squiggly_theta**k * gamma_func(k)) )

# Projections of gradient direction along and perpendicular to the axon direction
g_plel = dot(imager.g, model.axon_direction)
g_perp = sqrt(g_plel**2 + norm(imager.g)**2)

# Components of the q-vector in the axon's coordinate frame.
q_parallel = (imager.gamma * imager.G * g_plel * imager.little_delta) / (2.0 * pi)
q_perpendicular = (imager.gamma * imager.G * g_perp * imager.little_delta) / (2.0 * pi)

# Printing some simplified variables
print "Effective b-value:", imager.b
print "S-T exponent:     ", -imager.b * model.D_parallel
print "Expected signal:  ", (model.S_unweighted * exp(-imager.b * model.D_parallel))

S_out = model.S_unweighted * exp(-imager.b * model.D_parallel)
charmed_S = model.S_unweighted * S_cyl(.001)
print "Predicted signal with axon radius 1e-3 mm", charmed_S

# Integrating over the distribution of radii
S_in = 0.0
Gamma_sum = 0.0
width = 5.0e-6
for i in range(0, 300):
    middle_radius = width/2 + i * width
    #print "S_Cyl at", middle_radius, " is:", S_unweighted * S_cyl(middle_radius) * 5.0e-5
    S_in_middle = S_cyl(middle_radius) * Gamma(middle_radius, model.k, model.squiggly_theta)
    S_in += S_in_middle * width
    Gamma_sum += Gamma(middle_radius, model.k, model.squiggly_theta) * width

print ""
print "Total signal over range of radii:", S_in + S_out
print "Gamma sum", Gamma_sum

print ""