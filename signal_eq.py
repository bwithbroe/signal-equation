from numpy import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D

little_delta = .003 # seconds
big_delta = .02     # seconds
tau = .02           # seconds
gamma = 41.065e6    # Hz/T

g = matrix('.5; .1; .86023')
g_z = g[2]
g_xy = norm(g[0:2])

#
b = gamma**2.0 * g_z**2.0 * little_delta**2.0 * (big_delta - little_delta / 3.0)
print b

#
q_parallel = (gamma * g_z * little_delta) / (2.0 * pi)
q_perpendicular = (gamma * g_xy * little_delta) / (2.0 * pi)

#
D_parallel = 1.0 * 10.0**-3
D_perpendicular = 1.0 * 10.0**-3

def S_cyl(radius):
    print exp(-4.0 * pi**2.0 * q_parallel**2.0 * (big_delta - (little_delta/3.0)) * D_parallel)

    return (exp(-4.0 * pi**2.0 * q_parallel**2.0 * (big_delta - (little_delta/3.0)) * D_parallel) *
            exp(-4.0 * pi**2.0 * radius**4.0 * q_perpendicular**2.0 * D_perpendicular * tau * (7.0/96.0) *
                (2.0-(99.0/112.0)*(radius**2.0/(D_perpendicular*tau)))))

def Gamma(radius, k, sqiggly_theta):
    return ( (radius**(k-1.0) * exp(-(radius/squiggly_theta))) / (squiggly_theta**k * float(math.factorial(k))) )

S_in = 0.0

k = 2
squiggly_theta = 1.28


# for i in range(0, 100):
#     radius = i * 5.0 * 10.0**-4
#     S_in += S_cyl(radius) * Gamma(radius, k, squiggly_theta)
