# -*- coding: cp1252 -*-
"""
Mass-spring system:

Date of creation: 10. Feb 2016
Author: Thor Staerk Stenvang
"""
import numpy as np
from enum import Enum

CONST_GRAVITY       = -9.8      # Gravitational acceleration (m/s^2)
CONST_KD_DRAG       = 0.001     # Air drag damper coefficient (N*s/m)
CONST_KS_STRETCH    = 2.0e+3    # Spring constant (N/m)
CONST_KS_SHEAR      = 1.0e+3    # --||--
CONST_KS_BEND       = 1.0e+3    # --||--
CONST_KD_STRETCH    = 0.001     # Damping coefficient
CONST_KD_SHEAR      = 0.001     # --||--
CONST_KD_BEND       = 0.00      # --||--

#Constants for torsions spring
CONST_KS_TOR = 1.0e+4
CONST_KD_TOR = 0.5

class explicit_method(Enum):
    fe    = 0 # forward euler
    rk2   = 1 # Runga kutta 2
    rk4   = 2 # runga kutta 4
class spring_type(Enum):
    tension = 0
    torsion = 1
    
class Spring:
    """ Class for Springs connecting vertices in the spring system """
    def __init__(self, L0, ks, kd, I, J):
        self.l0_ = L0   # Rest length
        self.ks_ = ks   # Spring constant
        self.kd_ = kd   # Spring damping coeff
        self.indI_ = I  # index to first particle connected to spring
        self.indJ_ = J  # index to second ------||------

class TorSpring:
    """ class for Torsion Spring """
    def __init__(self,a,ks,kd,o,i,j,theta):
        self.a_ = a #Rest angel
        self.ks_ = ks #Spring constants
        self.kd_ = kd #Spring damping
        self.indI_ = i  # index to first particle connected to spring
        self.indJ_ = j  # index to second ------||------
        self.indO_ = o  # index to point of revolution
        self.angle_ = theta #Used for calculating damping
        
class Cloth:
    
    def __init__(self, dimX, dimY, rho, s):
        """
            Initialize
        """
        self.s = s
        self.dX = dimX # Mesh dimension
        self.dY = dimY # Mesh dimension
        self.nElements = (dimX*dimY) #Number of elements/Particles
        self.X = np.zeros((self.nElements, 3)) #Positions
        self.V = np.zeros((self.nElements, 3)) #Velocities
        self.F = np.zeros((self.nElements, 3)) #Forces
        self.mass = rho*s*s # Mass = density * spacing^2
        self.M = np.diag(np.zeros(3*self.nElements)+self.mass) #diagonal mass matrix
        self.sGravity = np.array([0.0,0.0,CONST_GRAVITY]) #Standard gravity acceleration
        self.Jx = np.zeros((3*self.nElements, 3*self.nElements)) #Force Jacobian. J(x), used by implicit methods
        self.Jv = np.zeros((3*self.nElements, 3*self.nElements)) #Force Jacobian. J(v), used by implicit methods
        self.constrIdx = np.array([])
        
        self.__UniMeshParticleCreator(dimX,dimY,s)
        print("Particles initialized")
        self.springs = []
        self.torSprings = []
        self.__NewSpringCreator(self.X,s)
        print("Springs initialized")

    def __UniMeshParticleCreator(self,dimY, dimX, pSpacing):
        #Create Particles in uniform mesh
        self.particles = []
        for i in range(0,dimY):
            for j in range(0,dimX):
                pos = np.asarray([i*pSpacing,j*pSpacing,1.0])
                self.X[i*dimX+j] = pos
        
    def __OldSpringCreator(self,X,pSpacing):
        #Create springs
        self.springs = []
        for i, xi in enumerate(self.X):
            for j, xj in enumerate(self.X):
                taxicab_dist = np.linalg.norm((xi-xj),ord=1)
                euclid_dist = np.linalg.norm((xi-xj))
                if((taxicab_dist <= (2*pSpacing) ) and i<j):
                    #Create springs!
                    if(taxicab_dist == pSpacing ):
                        #create stretch spring
                        spring = Spring(euclid_dist,CONST_KS_STRETCH,CONST_KD_STRETCH,i,j)
                        self.springs.append(spring)
                    if(euclid_dist < (2*pSpacing) and taxicab_dist == (2*pSpacing)):
                        #create shear springs
                        spring = Spring(euclid_dist,CONST_KS_SHEAR,CONST_KD_SHEAR,i,j)
                        self.springs.append(spring)
                    if(euclid_dist == (2*pSpacing) and taxicab_dist == (2*pSpacing)):
                        #create Bend springs
                        spring = Spring(euclid_dist,CONST_KS_BEND,CONST_KD_BEND,i,j)
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
    def __NewSpringCreator(self,X,pSpacing):
        #Create springs
        for i, xi in enumerate(self.X):
            for j, xj in enumerate(self.X):
                taxicab_dist = np.linalg.norm((xi-xj),ord=1)
                euclid_dist = np.linalg.norm((xi-xj))
                if((taxicab_dist <= (2*pSpacing) ) and i<j):
                    #Create springs!
                    if(taxicab_dist == pSpacing ):
                        #create stretch spring
                        spring = Spring(euclid_dist,CONST_KS_STRETCH,CONST_KD_STRETCH,i,j)
                        self.springs.append(spring)
                    if(euclid_dist < (2*pSpacing) and taxicab_dist == (2*pSpacing)):
                        #create shear springs
                        spring = Spring(euclid_dist,CONST_KS_SHEAR,CONST_KD_SHEAR,i,j)
                        self.springs.append(spring)
        #Create Torsion springs :)
        """ ONLY WORKS FOR GRID MESH STRUCTURE!!!
            
        """
        #Convert array to grid
        X3d = self.X.reshape((self.dY,self.dX,3))
        for y in range(0,self.dY):
            for x in range(0,self.dX):
                o  = ((y+0)*self.dY+x)
                xp = ((y+0)*self.dY+(x+1))
                xm = ((y+0)*self.dY+(x-1))
                yp = ((y+1)*self.dY+(x+0))
                ym = ((y-1)*self.dY+(x+0))
                #print([o,xp,xm,yp,ym]) #used for debugging
                
                if((y == 0 or y == (self.dY-1)) and (x == 0 or x == (self.dX-1))):
                    pass
                else:
                    if((y == 0 or y == (self.dY-1))):
                        tspring = TorSpring(0.0,CONST_KS_TOR,CONST_KD_TOR,o,xm,xp,0.0)
                        self.torSprings.append(tspring)
                    elif((x == 0 or x == (self.dX-1))):
                        tspring = TorSpring(0.0,CONST_KS_TOR,CONST_KD_TOR,o,ym,yp,0.0)
                        self.torSprings.append(tspring)
                    else:
                        tspring = TorSpring(0.0,CONST_KS_TOR,CONST_KD_TOR,o,xm,xp,0.0)
                        self.torSprings.append(tspring)
                        tspring = TorSpring(0.0,CONST_KS_TOR,CONST_KD_TOR,o,ym,yp,0.0)
                        self.torSprings.append(tspring)
    
        #print([X3d[1,0],self.X[self.dY]])
        
                    
    def Energy(self):
        """
            Calculates and returns the energy of the system
            Types of energy:
            
            Energy stored in a spring:
                Tension:        U_el1 = 1/2*k*x^2
                Torsion:        U_el2 = 1/2*k*theta^2
            Energy of gravity:  U_g = m*g*h
            Kinetic energy:     U_k =1/2*m*v^2
        """
        #Calculate energy for particles:
        #Gravitational energy:
        self.U_g = self.mass * -CONST_GRAVITY * self.X[:,2]
        #Kinetic energy
        vnorm = np.linalg.norm(self.V,axis=1) #Row wise norm
        self.U_k = 0.5 * self.mass * vnorm*vnorm

        #Calculate spring energy
        self.U_el1 = np.zeros(len(self.springs))
        self.U_el2 = np.zeros(len(self.torSprings))
        #Tension
        for i, s in enumerate(self.springs):
            xi = self.X[s.indI_]
            xj = self.X[s.indJ_]
            deltaX = xj-xi
            Xnorm = np.linalg.norm(deltaX)
            self.U_el1[i] = 0.5 * s.ks_ * (Xnorm-s.l0_)*(Xnorm-s.l0_)
        #Torsion
        for i, s in enumerate(self.torSprings):
            xi = self.X[s.indI_] #Pos of particle i
            xj = self.X[s.indJ_] #Pos of particle j
            xo = self.X[s.indO_] #Pos of particle k
            #Find angle between xij and xik
            a = xi-xo
            b = xo-xj
            a_u = a/np.linalg.norm(a) #unit vector
            b_u = b/np.linalg.norm(b) #unit vector
            angle = np.arccos(np.dot(a_u,b_u)) #angle between two unti vectors
            if np.isnan(angle):
                if(a_u == b_u).all():
                    angle = 0.0
                else:
                    angle = np.pi
            #calc energy
            self.U_el2[i] = 0.5 * s.ks_ * angle*angle
            
        return self.U_k.sum(), self.U_g.sum(), self.U_el2.sum(), self.U_el1.sum()
            
        
            
        
    def force(self,X,V,t):
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
        #self.F -= CONST_KD_DRAG*self.V
            
        #Add spring forces for tension springs
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
            spring_damp_force = s.kd_ * np.dot((vi-vj),deltaX/norm2)*deltaX/norm2  #Damping on spring
    
            #Add forces
            self.F[s.indI_] += spring_force +spring_damp_force
            self.F[s.indJ_] -= spring_force +spring_damp_force

        #Add spring forces for torsion springs
        for i, s in enumerate(self.torSprings): 
            xi = X[s.indI_] #Pos of particle i
            xj = X[s.indJ_] #Pos of particle j
            xo = X[s.indO_] #Pos of particle k
            #Find angle between xij and xik
            a = xi-xo
            b = xo-xj
            a_norm = np.linalg.norm(a)
            a_u = a/a_norm #unit vector
            b_norm = np.linalg.norm(b) 
            b_u = b/b_norm #unit vector
            angle = np.arccos(np.dot(a_u,b_u)) #angle between two unti vectors
            if np.isnan(angle):
                if(a_u == b_u).all():
                    angle = 0.0
                else:
                    #angle = np.pi
                    angle = 0.0
                    #print "huh?"
            if(angle > np.pi):
                print "test"
            #Check sign for angle
            cross = np.cross(a,b)
            normCro = np.linalg.norm(cross)
            if( normCro == 0.0):
                cross_hat = np.array([0.0,0.0,0.0])
            else:    
                cross_hat = cross/normCro #Torque direction vector

            #Calculate torque
            tau = s.ks_*angle*cross_hat
##            if(angle == 0):
##                print([s.indI_,s.indJ_])
##                print([self.X[s.indI_]])
##                print([self.X[s.indJ_]])
##                print([self.X[s.indO_]])

            dF1 = np.cross(tau,a)
            dF2 = np.cross(tau,b)
            #dF1 = dF1*[0.0,0.0,1.0]
            #dF2 = -dF2*[0.0,0.0,1.0]
            #Calculate force on particles
            self.F[s.indI_] += dF1
            self.F[s.indJ_] += dF2
            self.F[s.indO_] -= (dF1+dF2)

            # # # # # # # # # # # # # # # # #
            # Don't forget rotational friction!
            # tau_damp = k_d*omega
            # # # # # # # # # # # # # # # # #
            # First find the angular velocity: omega
            
            # We need to project on to plane define by cross_hat
            # Also has to be orthegonal to cross_hat and a or b
            
            # find orthegonal unit vector between cross_hat and a
            va = np.cross(a,cross_hat)
            if(np.linalg.norm(va) == 0):
                va_hat = np.array([0.0,0.0,0.0])
            else:    
                va_hat = va/np.linalg.norm(va)
            #Project Va onto v1_hat
            Va = np.dot(V[s.indI_],va_hat) * va_hat
            #Calculate angular vel
            omega_a = np.cross(a,Va)/(a_norm * a_norm)

            #Do the same for b
            vb = np.cross(b,cross_hat)
            if(np.linalg.norm(vb) == 0):
                vb_hat = np.array([0.0,0.0,0.0])
            else:
                vb_hat = vb/np.linalg.norm(vb)
            Vb = np.dot(V[s.indJ_],vb_hat) * vb_hat
            omega_b = np.cross(b,Vb)/(a_norm * a_norm)

            omega = omega_a + omega_b
            #print omega_a
            #print omega_b
            #print ""
            #print (angle - s.angle_)/t
            #omega = (angle - s.angle_)/t
            s.angle_ = angle
            tau_da = s.kd_*omega_a*cross_hat
            tau_db = s.kd_*omega_b*cross_hat

            dF_da = np.cross(tau_da,a)
            dF_db = np.cross(tau_db,b)
            self.F[s.indI_] += dF_da
            self.F[s.indJ_] += dF_db
            self.F[s.indO_] -= (dF_da+dF_db)
             
            
            

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
        self.force(self.X,self.V,stepT)
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
        self.force(Xtmp,Vtmp,stepT)
        b2 = self.F/self.mass

        #Step 3
        c1 = self.V + stepT/2*b2
        Xtmp = self.X+stepT/2*b1
        Vtmp = self.V+stepT/2*b2
        self.force(Xtmp,Vtmp,stepT)
        c2 = self.F/self.mass

        #Step 4
        d1 = self.V + stepT*c2
        Xtmp = self.X+stepT*c1
        Vtmp = self.V+stepT*c2
        self.force(Xtmp,Vtmp,stepT)
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

    def semiImplictEuler(self, stepT):
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


##c = Cloth(3,3,0.2,1.0)
##constr = np.array([0,2])
##c.constrain(constr)
##
##for i in range(0,5):
##    print""
##    #print c.Energy()
##    c.simUpdateExplicit(0.00014,explicit_method.fe)







