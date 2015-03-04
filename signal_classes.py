from numpy import *
from numpy.linalg import *
from scipy.special import gamma as gamma_func
from scipy.integrate import quad
from matplotlib.pyplot import *
from Model import *
from Imager import *
from random import *


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
    print q_parallel
    return exp(log_S_restr_plel + log_S_restr_perp)

def Gamma(radius, k, squiggly_theta):
    return ( (radius**(k-1.0) * exp(-(radius/squiggly_theta))) /
             (squiggly_theta**k * gamma_func(k)) )

# Algorithm from top answer to this question:
# http://stackoverflow.com/questions/6283080/random-unit-vector-in-multi-dimensional-space
def genVectors(num_vectors):
    vector_list = []
    vector_count = 0
    while vector_count < num_vectors:
        random_nums_array = array([0.0, 0.0, 0.0])
        for i in range(3):
            random_nums_array[i] = uniform(-1.0,1.0)
        
        r_squared = sum(random_nums_array**2)
        normalized_vector = array([0.0,0.0,0.0])
        if 0 < r_squared <=1.0:
            for i in range(3):
                normalized_vector[i] = random_nums_array[i]/sqrt(r_squared)
            vector_list.append(normalized_vector)
            vector_count += 1
    
    return vector_list

# Projections of gradient direction along and perpendicular to the axon direction
g_plel = dot(imager.g, model.axon_direction)
g_perp = sqrt(g_plel**2 + norm(imager.g)**2)

# Components of the q-vector in the axon's coordinate frame.
q_parallel = (imager.gamma * imager.G * g_plel * imager.little_delta) / (2.0 * pi)
q_perpendicular = (imager.gamma * imager.G * g_perp * imager.little_delta) / (2.0 * pi)

# Printing some simplified variables
print "\n==========================="
print "Effective b-value:", imager.b
print "S-T exponent:     ", -imager.b * model.D_parallel
print "Expected signal:  ", (model.S_unweighted * exp(-imager.b * model.D_parallel))

S_out = exp(-imager.b * model.D_parallel)
charmed_S = model.S_unweighted * S_cyl(1e-3)
print
print "Predicted signal with axon radius 1e-3 mm", charmed_S
print "Calculated S_out:", S_out

# Integrating over the distribution of radii
def S_in_instance(radius):
    return S_cyl(radius)*Gamma(radius, model.k, model.squiggly_theta)

S_in = quad(S_in_instance, 0, .008)[0]
S = model.S_unweighted*(model.f*S_in + (1-model.f)*S_out)

print "Total S_in over range of radii:", S_in
print "Combined Signal:", S

num_vectors = 3
vector_list = genVectors(num_vectors)
orientation_list = []
axon_1_signal_list = []

# Recomputes the variables with the new axon directions and gradient directions
def recomputeVariables(gradient_direction):
    g_plel = dot(gradient_direction, model.axon_direction)
    g_perp = sqrt(g_plel**2 + norm(gradient_direction)**2)
    
    q_parallel = (imager.gamma * imager.G * g_plel * imager.little_delta) / (2.0 * pi)
    q_perpendicular = (imager.gamma * imager.G * g_perp * imager.little_delta) / (2.0 * pi)
    
    print q_parallel

    S_in = quad(S_in_instance, 0, .008)[0]
    S = model.S_unweighted*(model.f*S_in + (1-model.f)*S_out)
    return S

# Calculating signal at each gradient
for vector in vector_list:
    axon_1_signal_list.append(recomputeVariables(vector))

orientation_list.append(axon_1_signal_list)

# print orientation_list[0]


print ""

# Graph of the signal equation and Gamma over the range of integration
'''
radii = arange(0, 0.0045, 0.00001)
prob = Gamma(radii, model.k, model.squiggly_theta)
signal = S_cyl(radii)

figure()
subplot(211)
plot(radii, prob)

subplot(212)
plot(radii, signal)
show()
'''

