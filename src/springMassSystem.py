"""
Mass-spring system:

Date of creation: 10. Feb 2016
Author: Thor Staerk Stenvang
"""
import numpy as np

class Spring:
    """ Class for Springs connecting vertices in the spring system """

    def __init__(self, L0, ks, kd, I, J):
        self.l0_ = L0   # Rest length
        self.ks_ = ks    # Spring constant
        self.kd_ = kd  # Spring damping coeff
        self.indI_ = I # index to first particle connected to spring
        self.indJ_ = J # index to second ------||------
        self.sForce = np.empty(3) #Spring force
        self.dForce = np.empty(3) #Spring damping force

    def calculateSpringForce(self,p): #p: Particles
        """
        Calculates the spring force acting on particle i and j
        """
        xi = p[self.indI_].pos3D_
        xj = p[self.indJ_].pos3D_
        self.sForce = self.ks_ * (xj-xi)/np.linalg.norm(xj-xi) *(np.linalg.norm(xj-xi)-self.l0_)
        return self.sForce
    
    def calculateDampForce(self,p): #p: Particles
        """
        Calculates the spring damping force acting on particle i and j
        """
        xi = p[self.indI_].pos3D_
        xj = p[self.indJ_].pos3D_
        vi = p[self.indI_].vel3D_
        vj = p[self.indJ_].vel3D_
        self.dForce = self.kd_ * (vj-vi)*(xj-xi)/np.linalg.norm(xj-xi)
        return self.dForce
        


class Particle:
    """ Class for each particle located at each vertex """

    def __init__(self, m, pos, vel):
        self.mass_ = m   # Mass
        self.pos3D_ = pos    # Position
        self.vel3D_ = vel    # Velocity
        self.force3D_ = np.array([0.0,0.0,0.0])  # Force
        self.springIdx_ = []
        self.constrained = False

class Cloth:
    
    def __init__(self, dimX, dimY, m):
        """
            Initialize
        """
        self.sGravity = np.array([0.0,0.0,-9.8]) #Standard gravity acceleration
        #Create Particles in uniform mesh
        self.particles = []
        for i in range(0,dimY):
            for j in range(0,dimX):
                pos = np.asarray([i,j,1.0])
                p = Particle(m,pos,np.array([0.0,0.0,0.0]))
                self.particles.append(p)
        print("Particles initialized")
        #Create springs
        self.springs = []
        for i, pI in enumerate(self.particles):
            xi = pI.pos3D_
            for j, pJ in enumerate(self.particles):
                xj = pJ.pos3D_
                taxicab_dist = np.linalg.norm((xi-xj),ord=1)
                euclid_dist = np.linalg.norm((xi-xj))
                if((taxicab_dist <= 2) and i<j):
                    #Create springs!
                    if(taxicab_dist == 1 ):
                        #create stretch spring
                        spring = Spring(1.0,5.0e+4,0.8,i,j)
                        self.springs.append(spring)
                        pI.springIdx_.append(len(self.springs)-1)
                        pJ.springIdx_.append(len(self.springs)-1)
                    if(euclid_dist < 2 and taxicab_dist == 2):
                        #create shear springs
                        spring = Spring(euclid_dist,0.5e+4,0.2,i,j)
                        self.springs.append(spring)
                        pI.springIdx_.append(len(self.springs)-1)
                        pJ.springIdx_.append(len(self.springs)-1)
                    if(euclid_dist == 2.0 and taxicab_dist == 2.0):
                        #create Bend springs
                        spring = Spring(euclid_dist,0.05e+7,0.15,i,j)
                        self.springs.append(spring)
                        pI.springIdx_.append(len(self.springs)-1)
                        pJ.springIdx_.append(len(self.springs)-1)
        print("Springs initialized")
                
    def simUpdate1(self,stepT):
        """
            simulation update: Uses explicit Euler
        """
        #CALCULATE FORCES
        for s in self.springs:
            #Update spring forces
            s.calculateSpringForce(self.particles)
            s.calculateDampForce(self.particles)
            
        for i, p in enumerate(self.particles):
            #Calculate forces acting on particle
            #Spring forces
            fs = np.array([0.0,0.0,0.0])
            fd = np.array([0.0,0.0,0.0])
            for sIdx in p.springIdx_:
                #Check sign.
                if(i == self.springs[sIdx].indI_):
                    fs += self.springs[sIdx].sForce
                    fd += self.springs[sIdx].dForce
                elif(i == self.springs[sIdx].indJ_):
                    fs -= self.springs[sIdx].sForce
                    fd -= self.springs[sIdx].dForce
                else:
                    print("ERROR. This should not happen")
            #Gravity
            fg = p.mass_ * self.sGravity
            #Total force
            p.force3D_ = fs+fd+fg
            #p.force3D_ = fg
            #p.force3D_ = fs+fd
            
        #UPDATE POS AND VEL
        for p in self.particles:
            if p.constrained == False:
                p.vel3D_ += stepT*p.force3D_/p.mass_
                p.pos3D_ += stepT*p.vel3D_

    def simUpdate2(self,stepT):
        """
            simulation update: Uses second order Runga Kutta
            NOT WORKING! IS BEING DEVELOPED
        """
        #CALCULATE FORCES
        for s in self.springs:
            #Update spring forces
            s.calculateSpringForce(self.particles)
            s.calculateDampForce(self.particles)

        a1 = []
        a2 = []
        for i, p in enumerate(self.particles):
            #Calculate forces acting on particle
            #Spring forces
            fs = np.array([0.0,0.0,0.0])
            fd = np.array([0.0,0.0,0.0])
            for sIdx in p.springIdx_:
                #Check sign.
                if(i == self.springs[sIdx].indI_):
                    fs += self.springs[sIdx].sForce
                    fd += self.springs[sIdx].dForce
                elif(i == self.springs[sIdx].indJ_):
                    fs -= self.springs[sIdx].sForce
                    fd -= self.springs[sIdx].dForce
                else:
                    print("ERROR. This should not happen")
            #Gravity
            fg = p.mass_ * self.sGravity

            # update a1 and a2
            a1.append(p.vel3D_)
            a2.append((fs+fd+fg)/p.mass_)

        #UPDATE POS AND VEL
        for i, p in enumerate(self.particles):
            b1 = p.vel3D_ + stepT/2 * a1[i]
            #b2 = (p.mass_ * self.sGravity
            
            if p.constrained == False:
                p.vel3D_ += stepT*p.force3D_/p.mass_
                p.pos3D_ += stepT*p.vel3D_
            
    def constrainParticle(self,index):
        self.particles[index].constrained = True
        
"""
Test
c = Cloth(3,3,1)
print(len(c.springs))
c.simUpdate2(0.001)
c.simUpdate1(0.002)
for j in range(0,2):
##    for i in range(0,len(c.particles)):
##        print(c.particles[i].pos3D_)
##        print(c.particles[i].springIdx_)
##        print(c.springs[i].indI_)
##        print(c.springs[i].indJ_)
    c.simUpdate1(0.02)
"""

