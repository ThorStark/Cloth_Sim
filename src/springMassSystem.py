"""
Mass-spring system:

Date of creation: 10. Feb 2016
Author: Thor Staerk Stenvang
"""
import numpy as np
from enum import Enum

class explicit_method(Enum):
    forward_euler = 0
    runge_kutta_2   = 1
    runge_kutta_4   = 2
    
class Spring:
    """ Class for Springs connecting vertices in the spring system """

    def __init__(self, L0, ks, kd, I, J):
        self.l0_ = L0   # Rest length
        self.ks_ = ks   # Spring constant
        self.kd_ = kd   # Spring damping coeff
        self.indI_ = I  # index to first particle connected to spring
        self.indJ_ = J  # index to second ------||------    

class Cloth:
    
    def __init__(self, dimX, dimY, m):
        """
            Initialize
        """
        self.nElements = (dimX*dimY) #Number of elements/Particles
        self.X = np.zeros((self.nElements, 3)) #Positions
        self.V = np.zeros((self.nElements, 3)) #Velocities
        self.F = np.zeros((self.nElements, 3)) #Forces
        self.mass = m
        self.sGravity = np.array([0.0,0.0,-9.8]) #Standard gravity acceleration
        
        #Create Particles in uniform mesh
        self.particles = []
        for i in range(0,dimY):
            for j in range(0,dimX):
                pos = np.asarray([i,j,1.0])
                self.X[i*dimY+j] = pos
        print("Particles initialized")
        
        #Create springs
        self.springs = []
        for i, xi in enumerate(self.X):
            for j, xj in enumerate(self.X):
                taxicab_dist = np.linalg.norm((xi-xj),ord=1)
                euclid_dist = np.linalg.norm((xi-xj))
                if((taxicab_dist <= 2) and i<j):
                    #Create springs!
                    if(taxicab_dist == 1 ):
                        #create stretch spring
                        spring = Spring(1.0,5.0e+4,0.8,i,j)
                        self.springs.append(spring)
                    if(euclid_dist < 2 and taxicab_dist == 2):
                        #create shear springs
                        spring = Spring(euclid_dist,0.5e+4,0.2,i,j)
                        self.springs.append(spring)
                    if(euclid_dist == 2.0 and taxicab_dist == 2.0):
                        #create Bend springs
                        spring = Spring(euclid_dist,0.05e+7,0.15,i,j)
                        self.springs.append(spring)
        print("Springs initialized")
                
    def force(self,X,V):
        """
            Calculates the force on all particels in the system
        """
        #Gravity
        for i in range(0,self.nElements):
            #Clear forces
            self.F[i] = [0.0, 0.0, 0.0]
            #Add gravitational force
            self.F[i] += self.sGravity
        #Add spring forces
        for s in self.springs:
            """
            Calculates the spring force acting on particle i and j
            """
            xi = X[s.indI_] #Pos of particle i
            xj = X[s.indJ_] #Pos of particle j
            vi = V[s.indI_] #velocity of particle i
            vj = V[s.indJ_] #velocity of particle j
            deltaX = xj-xi
            norm2 = np.linalg.norm(deltaX)
            spring_force = s.ks_ * deltaX/norm2 *(norm2-s.l0_)     # Spring force
            spring_damp_force = -s.kd_ * np.dot((vj-vi),deltaX)/norm2  #Damping on string
            #Add forces
            self.F[s.indI_] += spring_force+spring_damp_force
            self.F[s.indJ_] -= spring_force+spring_damp_force
        #Check if constrained
        for c in self.constrIdx:
            self.F[c] = [0.0,0.0,0.0]
    def simUpdateExplicit(self,stepT,method):
        """
            simulation update: Uses explicit Euler
        """
        #Compute forces
        self.force(self.X,self.V)
        #Choose method
        if method == explicit_method.forward_euler:
            self.forwardEuler(stepT)
        elif method == explicit_method.runge_kutta_2:
            self.RK2(stepT)
        elif method == explicit_method.runge_kutta_4:
            print "Not implemented"

        

    def forwardEuler(self, stepT):
        oldV = self.V
        self.V += stepT*self.F/self.mass
        self.X += stepT*oldV
            

    def RK2(self, stepT):      
        a1 = self.V
        a2 = self.F/self.mass
        
        b1 = self.V + stepT/2*a2
        X2 = self.X+stepT/2*a1
        V2 = self.V+stepT/2*a2
        self.force(X2,V2)
        b2 = self.F/self.mass

        self.X += stepT*b1
        self.V += stepT*b2

    def RK4(self, stepT):
        for c in self.constrIdx:
            self.F[c] = [0.0,0.0,0.0]
        
        
    def constrain(self,constrIdx):
        self.constrIdx = constrIdx
        a1 = self.V
        a2 = self.F/self.mass

#Test
        """
c = Cloth(3,3,1)
for i in range(0,10):
    c.simUpdateExplicit(0.002,explicit_method.runge_kutta_2)
    print c.F
    print ""
    """


