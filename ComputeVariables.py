from numpy import *
from numpy.linalg import *
from scipy.special import gamma as gamma_func
from scipy.integrate import quad
from matplotlib.pyplot import *
from Model import *
from Imager import *
from random import *

class ComputeVariables:

    def __init__(self, imager, model):
        self.imager = imager
        
        # Projections of gradient direction along and perpendicular to the axon direction
        self.g_plel = dot(imager.g, model.axon_direction)
        self.g_perp = sqrt(self.g_plel**2 + norm(imager.g)**2)
        
        # Components of the q-vector in the axon's coordinate frame.
        self.q_parallel = (imager.gamma * imager.G * self.g_plel * imager.little_delta) / (2.0 * pi)
        self.q_perpendicular = (imager.gamma * imager.G * self.g_perp * imager.little_delta) / (2.0 * pi)
        
        self.S = self.calcS(model)
    
    def updateVariables(self, model):
        self.g_plel = dot(imager.g, model.axon_direction)
        self.g_perp = sqrt(self.g_plel**2 + norm(imager.g)**2)
    
        self.q_parallel = (imager.gamma * imager.G * self.g_plel * imager.little_delta) / (2.0 * pi)
        self.q_perpendicular = (imager.gamma * imager.G * self.g_perp * imager.little_delta) / (2.0 * pi)

        self.S = self.calcS(model)
    
    def S_cyl(self, radius, model):
        # Assaf et al., MRM 2004
        log_S_restr_plel = -4.0 * pi**2.0 * self.q_parallel**2.0 * (imager.big_delta - (imager.little_delta/3.0)) * model.D_parallel
        # BUG --- this should reduce to the previous equation as radius -> inf
        log_S_restr_perp = (-((4.0 * pi**2.0 * radius**4.0 * self.q_perpendicular**2.0) / (model.D_perpendicular * imager.tau)) *
                            (7.0/96.0) * (2.0-(99.0/112.0)*(radius**2.0/(model.D_perpendicular * imager.tau))))
        # print "rest_plel exponent:", log_S_restr_plel
        # print "rest_perp exponent:", log_S_restr_perp
        return exp(log_S_restr_plel + log_S_restr_perp)
        
    def S_restr_plel(self, model):
        return exp( -4.0 * pi**2.0 * self.q_parallel**2.0 * (imager.big_delta - (imager.little_delta/3.0)) * model.D_parallel )
        
    def S_restr_perp(self, radius, model):
        return exp((-((4.0 * pi**2.0 * radius**4.0 * self.q_perpendicular**2.0) / (model.D_perpendicular * imager.tau)) *
                            (7.0/96.0) * (2.0-(99.0/112.0)*(radius**2.0/(model.D_perpendicular * imager.tau)))))
                            
    def ST(self, model):
        return exp(-imager.b * model.D_parallel)
                            
    def Gamma(self, radius, k, squiggly_theta):
        return ( (radius**(k-1.0) * exp(-(radius/squiggly_theta))) /
                 (squiggly_theta**k * gamma_func(k)) )
                 
    def calcSout(self, model):
        return exp(-imager.b * model.D_parallel)
        
    def calcSin(self, model):
        return quad(self.S_in_instance, 0, .008, args=(model))[0]
    
    def calcS(self, model):
        S_in = self.calcSin(model)
        S_out = self.calcSout(model)
        return model.S_unweighted*(model.f*S_in + (1-model.f)*S_out)
    
    def S_in_instance(self, radius, model):
        return self.S_cyl(radius, model)*self.Gamma(radius, model.k, model.squiggly_theta)
        
    def genRandomVectors(self, num_vectors):
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
        
    def genVectorsWithPhi(self, phi, num_vectors):
        vector_list = []
        z = cos(phi)
        for i in range(num_vectors):
            theta = 2*pi * i / (num_vectors - 1)
            x = cos(theta) * sin(phi)
            y = sin(theta) * sin(phi)
            new_vector = array([x, y, z])
            vector_list.append(new_vector)
        return vector_list
    
    def genVectorsWithTheta(self, theta, num_vector):
        vector_list = []
        for i in range(num_vectors):
            phi = (pi/2) * i / (num_vectors - 1)
            z = cos(phi)
            x = cos(theta) * sin(phi)
            y = sin(theta) * sin(phi)
            new_vector = array([x, y, z])
            vector_list.append(new_vector)
        return vector_list
        
    def compareSignals(self, orig_axon_signal_list):
        found = False
        num_checked = 0
        while not found:
            found = True
            
            new_model = Model()
            new_model.S_unweighted = uniform(100,500)
            new_model.axon_direction = self.genVectors(1)[0]
            
            new_axon_signal_list = []
            for vector in vector_list:
                self.imager.g = vector
                self.updateVariables(new_model)
                new_axon_signal_list.append(self.calcS(new_model))
            for i in range(len(orig_axon_signal_list)):
                # print "orig:", orig_axon_signal_list[i], " new:", new_axon_signal_list[i]
                if abs(orig_axon_signal_list[i] - new_axon_signal_list[i]) > 10:
                    found = False
                    break
            num_checked += 1
            print "nope", num_checked

        print "ta-da!"
        return new_model

    
    def plotSignal(self, model):
        radii = arange(0, 0.005, 0.00001)
        signal_plel = empty(500)
        signal_plel.fill(self.S_restr_plel(model))
        signal_perp = self.S_restr_perp(radii, model)
        signal_ST = empty(500)
        signal_ST.fill(self.ST(model))
        
        figure()
        
        plot(radii, signal_plel)
        plot(radii, signal_perp)
        plot(radii, signal_ST)
        
        show()
        
    
        
if __name__ == '__main__':
    imager = Imager()
    model = Model()
    
    c = ComputeVariables(imager, model)
    num_vectors = 10
    vector_list = c.genVectorsWithTheta(0, num_vectors)
    
    axon_1_signal_list = []
    for vector in vector_list:
        c.imager.g = vector
        c.updateVariables(model)
        axon_1_signal_list.append(c.calcS(model))
        print vector
        c.plotSignal(model)
    #print axon_1_signal_list
    
    
    '''
    y = array(axon_1_signal_list)
    x = arange(0, pi/2, pi/2 / num_vectors)
    
    figure()
    plot(x, y)
    show()
    '''
    
    
    '''
    print
    new_model = c.compareSignals(axon_1_signal_list)
    print model.axon_direction, "vs.", new_model.axon_direction
    print model.S_unweighted, "vs.", new_model.S_unweighted
    '''