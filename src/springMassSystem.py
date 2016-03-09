"""
Mass-spring system:

Date of creation: 10. Feb 2016
Author: Thor Staerk Stenvang
"""
import numpy as np
from enum import Enum

CONST_GRAVITY       = -9.8      # Gravitational acceleration (m/s^2)
CONST_KD_DRAG       = 0.02       # Air drag damper coefficient (N*s/m)
CONST_KS_STRETCH    = 2.0e+3    # Spring constant (N/m)
CONST_KS_SHEAR      = 1.0e+3    # --||--
CONST_KS_BEND       = 1.0e+3    # --||--
CONST_KD_STRETCH    = 0.001      # Damping coefficient
CONST_KD_SHEAR      = 0.001      # --||--
CONST_KD_BEND       = 0.001      # --||--

class explicit_method(Enum):
    fe    = 0 # forward euler
    rk2   = 1 # Runga kutta 2
    rk4   = 2 # runga kutta 4
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
        self.type_  = t # Type of spring (tension or torsion)

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
        self.M = np.diag(np.zeros(3*self.nElements)+self.mass) #diagonal mass matrix
        self.sGravity = np.array([0.0,0.0,CONST_GRAVITY]) #Standard gravity acceleration
        self.Jx = np.zeros((3*self.nElements, 3*self.nElements)) #Force Jacobian. J(x), used by implicit methods
        self.Jv = np.zeros((3*self.nElements, 3*self.nElements)) #Force Jacobian. J(v), used by implicit methods
        self.constrIdx = np.array([])
        
        #Create Particles in uniform mesh
        self.particles = []
        for i in range(0,dimY):
            for j in range(0,dimX):
                pos = np.asarray([i,j,1.0])
                self.X[i*dimX+j] = pos
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
        #Add damping force (Air drag)
        self.F -= CONST_KD_DRAG*self.V
            
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
                spring_damp_force = s.kd_ * np.dot((vi-vj),deltaX/norm2)*deltaX/norm2  #Damping on string
                #print([np.linalg.norm(spring_damp_force), np.linalg.norm(spring_force)]) #test magnitude of force
                #print([np.dot(spring_damp_force,spring_force),(np.linalg.norm(spring_damp_force)*(np.linalg.norm(spring_force)))]) #Test direction of forces
                """
                if(np.linalg.norm(spring_damp_force) > np.linalg.norm(spring_force)):
                    print("K ",spring_force)
                    print("d ",spring_damp_force)
                    print("s ",s.indI_,s.indJ_)
                    print ""
                """
                
                #Add forces
                self.F[s.indI_] += spring_force +spring_damp_force
                self.F[s.indJ_] -= spring_force +spring_damp_force
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

    def forceDerivatives(self,X,V):
        self.Jx = np.zeros((3*self.nElements, 3*self.nElements))
        self.Jv = np.zeros((3*self.nElements, 3*self.nElements))
        for s in self.springs:
            xi = X[s.indI_] #Pos of particle i
            xj = X[s.indJ_] #Pos of particle j
            deltaX = xj-xi
            deltaXT = np.transpose(deltaX)
            norm2 = np.linalg.norm(deltaX)
            I = np.identity(3)
            
            #Position jacobian
            jx_sub = s.ks_*(-I + s.l0_/norm2*(I - np.outer(deltaX,deltaX)/(norm2*norm2)))
            """
            Insert into jacobian
            Jx = | Jx_ii  Jx_ij | = | j_sub  -j_sub |
                 | Jx_ji  Jx_jj |   | -j_sub  j_sub |
            """
            self.Jx[3*s.indI_:3*s.indI_+3, 3*s.indI_:3*s.indI_+3] += jx_sub
            self.Jx[3*s.indI_:3*s.indI_+3, 3*s.indJ_:3*s.indJ_+3] += -jx_sub
            self.Jx[3*s.indJ_:3*s.indJ_+3, 3*s.indI_:3*s.indI_+3] += -jx_sub
            self.Jx[3*s.indJ_:3*s.indJ_+3, 3*s.indJ_:3*s.indJ_+3] += jx_sub

            #Velocity jacobian
            """
            Todo
            """
            
            
              

    def simUpdateExplicit(self,stepT,method):
        """
            simulation update: Uses explicit Euler
        """
        #Compute forces
        self.force(self.X,self.V)
        #Choose method
        if method == explicit_method.fe:
            self.forwardEuler(stepT)
        elif method == explicit_method.rk2:
            self.RK2(stepT)
        elif method == explicit_method.rk4:
            self.RK4(stepT)

    def forwardEuler(self, stepT):
        oldV = self.V
        self.V += stepT*self.F/self.mass
        #Check if constrained
        for c in self.constrIdx:
            self.V[c] = [0.0,0.0,0.0]
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
        deltaX = stepT*b1
        deltaV = stepT*b2
        #Check if constrained
        for c in self.constrIdx:
            deltaV[c] = [0.0,0.0,0.0]
            deltaX[c] = [0.0,0.0,0.0]
        self.X += deltaX
        self.V += deltaV

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
        deltaX = stepT/6*(a1+2*b1+2*c1+d1)
        deltaV = stepT/6*(a2+2*b2+2*c2+d2)
        #Check if constrained
        for c in self.constrIdx:
            deltaV[c] = [0.0,0.0,0.0]
            deltaX[c] = [0.0,0.0,0.0]
        self.X += deltaX
        self.V += deltaV

    def ImplictEuler(self, stepT):
        """
            (M-h*df/dv-h^2*df/dx)deltaV=h(f0+h*df/dx*v0)
            A*deltaV = b
            A = (M-h*df/dv-h^2*df/dx)
            b = h(f0+h*df/dx*v0)
        """
        #Calculate forces and jacobians
        self.force(self.X,self.V)
        self.forceDerivatives(self.X,self.V)
        #setup linear equation
        A = self.M - stepT*self.Jv - stepT*stepT*self.Jx
        b = stepT*(self.F.flatten() + stepT * np.dot(self.Jx,self.V.flatten()))
        
        #Solve equation
        deltaV = np.linalg.solve(A,b)

        #reshape
        deltaV = deltaV.reshape((-1,3))

        #Check if constrained
        for c in self.constrIdx:
            deltaV[c] = [0.0,0.0,0.0]
        

        #Update pos and vel
        self.X += stepT*(self.V + deltaV)
        self.V += deltaV
        
        
    def constrain(self,constrIdx):
        self.constrIdx = constrIdx #Index to constrains

"""
c = Cloth(2,2,0.2)
constr = np.array([0,2])
c.constrain(constr)
for i in range(0,5):
    print""
    c.simUpdateExplicit(0.00014,explicit_method.fe)
"""







