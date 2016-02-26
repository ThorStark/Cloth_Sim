"""
Mass-spring system:

Date of creation: 10. Feb 2016
Author: Thor Staerk Stenvang
"""
import numpy as np
from enum import Enum

CONST_GRAVITY       = -9.8      # Gravitational acceleration (m/s^2)
CONST_KS_STRETCH    = 5.0e+4    # Spring constant (N/m)
CONST_KS_SHEAR      = 0.5e+4    # --||--
CONST_KS_BEND       = 0.05e+7   # --||--
CONST_KD_STRETCH    = 0.8       # Damping coefficient (N*s/m)
CONST_KD_SHEAR      = 0.2       # --||--
CONST_KD_BEND       = 0.2       # --||--

class explicit_method(Enum):
    forward_euler = 0
    runge_kutta_2   = 1
    runge_kutta_4   = 2
class spring_type(Enum):
    tension = 0
    torsion = 1
    
class Spring:
    """ Class for Springs connecting vertices in the spring system """
    def __init__(self, L0, ks, kd, I, J, K, t):
        self.l0_ = L0   # Rest length
        self.ks_ = ks   # Spring constant
        self.kd_ = kd   # Spring damping coeff
        self.indI_ = I  # index to first particle connected to spring
        self.indJ_ = J  # index to second ------||------
        self.indK_ = K  # index two third particle (used for torsion spring)
        self.type_  = t  # Type of spring (tension or torsion)

class Cloth:
    
    def __init__(self, dimX, dimY, m):
        """
            Initialize
        """
        self.nElements = (dimX*dimY) #Number of elements/Particles
        self.X = np.zeros((self.nElements, 3)) #Positions
        self.V = np.zeros((self.nElements, 3)) #Velocities
        self.F = np.zeros((self.nElements, 3)) #Forces
        self.mass = m # Mass
        self.M = np.diag(np.zeros(self.nElements)+self.mass) #diagonal mass matrix
        self.Minv = np.linalg.inv(self.M) #Invers Mass matrix
        self.sGravity = np.array([0.0,0.0,CONST_GRAVITY]) #Standard gravity acceleration
        self.constrIdx = np.array([])
        
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
                        spring = Spring(euclid_dist,CONST_KS_STRETCH,CONST_KD_STRETCH,i,j,-1,spring_type.tension)
                        self.springs.append(spring)
                    if(euclid_dist < 2 and taxicab_dist == 2):
                        #create shear springs
                        spring = Spring(euclid_dist,CONST_KS_SHEAR,CONST_KD_SHEAR,i,j,-1,spring_type.tension)
                        self.springs.append(spring)
                    if(euclid_dist == 2.0 and taxicab_dist == 2.0):
                        #create Bend springs
                        spring = Spring(euclid_dist,CONST_KS_BEND,CONST_KD_BEND,i,j,-1,spring_type.tension)
                        self.springs.append(spring)
                        """
                        #TODO. Finish new implementation of bend torsion spring.
                        #create torsion bend springs.
                        for k, xk in enumerate(self.X):
                            #find particle between the two
                            taxicab_dist1 = np.linalg.norm((xi-xk),ord=1)
                            taxicab_dist2 = np.linalg.norm((xj-xk),ord=1)
                            if(taxicab_dist1 == 1 and taxicab_dist2 == 1):
                                spring = Spring(0,0.05e+7,0.15,i,j,k,spring_type.torsion)
                                self.springs.append(spring)
                        """
                    
                        
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
            if(s.type_ == spring_type.tension):
                xi = X[s.indI_] #Pos of particle i
                xj = X[s.indJ_] #Pos of particle j
                vi = V[s.indI_] #velocity of particle i
                vj = V[s.indJ_] #velocity of particle j
                deltaX = xj-xi
                norm2 = np.linalg.norm(deltaX)
                spring_force = s.ks_ * deltaX/norm2 *(norm2-s.l0_)     # Spring force
                spring_damp_force = -s.kd_ * np.dot((vj-vi),deltaX/norm2)  #Damping on string
                #Add forces
                self.F[s.indI_] += spring_force+spring_damp_force
                self.F[s.indJ_] -= spring_force+spring_damp_force
            elif(s.type_ == spring_type.torsion):
                xi = X[s.indI_] #Pos of particle i
                xj = X[s.indJ_] #Pos of particle j
                xk = X[s.indK_] #Pos of particle k
                #Find angle between xij and xik
                xij = xi-xj
                xik = xi-xk
                xij_u = xij/np.linalg.norm(xij) #unit vector
                xik_u = xik/np.linalg.norm(xik) #unit vector
                angle = np.arccos(np.dot(xij_u,xik_u))
                if np.isnan(angle):
                    if(xij_u == xik_u).all():
                        angle = 0.0
                    else:
                        angle = np.pi
                #Calculate force
                #spring_torsion += 
                #print angle

        #Check if constrained
        for c in self.constrIdx:
            self.F[c] = [0.0,0.0,0.0]

    def forceDerivatives(self,X,V):
        for s in self.springs:
            xi = X[s.indI_] #Pos of particle i
            xj = X[s.indJ_] #Pos of particle j
            vi = V[s.indI_] #velocity of particle i
            vj = V[s.indJ_] #velocity of particle j
            deltaX = xj-xi
            deltaXtransposed = np.transpose(deltaX)
            norm2 = np.linalg.norm(deltaX)
            I = np.identity(3)
            #Position jacobian
            
        

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
            self.RK4(stepT)

    def forwardEuler(self, stepT):
        oldV = self.V
        self.V += stepT*self.F/self.mass
        self.X += stepT*oldV  

    def RK2(self, stepT):
        #Step 1
        a1 = self.V
        a2 = self.F/self.mass

        #Step 2
        b1 = self.V + stepT/2*a2
        X2 = self.X+stepT/2*a1
        V2 = self.V+stepT/2*a2
        self.force(X2,V2)
        b2 = self.F/self.mass

        #Update Pos and Vel
        self.X += stepT*b1
        self.V += stepT*b2

    def RK4(self, stepT):
        #Step 1
        a1 = self.V
        a2 = self.F/self.mass

        #Step 2
        b1 = self.V + stepT/2*a2
        Xtmp = self.X+stepT/2*a1
        Vtmp = self.V+stepT/2*a2
        self.force(Xtmp,Vtmp)
        b2 = self.F/self.mass

        #Step 3
        c1 = self.V + stepT/2*b2
        Xtmp = self.X+stepT/2*b1
        Vtmp = self.V+stepT/2*b2
        self.force(Xtmp,Vtmp)
        c2 = self.F/self.mass

        #Step 4
        d1 = self.V + stepT*c2
        Xtmp = self.X+stepT*c1
        Vtmp = self.V+stepT*c2
        self.force(Xtmp,Vtmp)
        d2 = self.F/self.mass

        #Update Pos and Vel
        self.X += stepT/6*(a1+2*b1+2*c1+d1)
        self.V += stepT/6*(a2+2*b2+2*c2+d2)

    def ImplictEuler(self, stepT):
        """
            (I - h*M^-1 * df/dv - h2*M^-1*df_dx)deltaV = h*M^-1(f0+h*df/dx*v0)
            A*deltaV = b
            A = (I - h*M^-1 * df/dv - h2*M^-1*df_dx)
            b = h*M^-1(f0+h*df/dx*v0)
        """
        
        
        
    def constrain(self,constrIdx):
        self.constrIdx = constrIdx #Index to constrains

"""
c = Cloth(3,3,0.2)
for i in range(0,200):
    c.simUpdateExplicit(0.005,explicit_method.forward_euler)
    print ""
print len(c.springs)
"""

