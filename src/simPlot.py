"""
simPlot.py
Thor Stark Stenvang
Date of creation 17. feb 2016
"""
import springMassSystem as sms
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm #Color map
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec
import matplotlib.tri as tri
import time
#Data plot update function

CONST_STEP = 0.0002
CONST_FRAMES = 100

def data(i,c,surf,triang):
    global t
    global xp
    global yp1
    global yp2
    global yp3
    global yp4
    #global ypt
    
    #Obtain coordinates of points
    x = c.X[:,0]
    y = c.X[:,1]
    z = c.X[:,2]
    
    ax1.clear() # Clear plot
    ###ax_tmp = fig.gca(projection='3d')
    ###plt.hold(True)
    x_s = np.transpose(x.reshape((5,5)))
    y_s = np.transpose(y.reshape((5,5)))
    z_s = np.transpose(z.reshape((5,5)))
    ###ax_tmp.plot_surface(x_s, y_s, z_s,rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0.4)
    ###plt.show()
    #Plot tri mesh
    surf = ax1.plot_trisurf(x,y,z,triangles=triang.triangles,cmap=cm.jet,linewidth=0.4)
    #surf = ax1.plot_surface(x_s, y_s, z_s,rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0.4)
    # Limit plot for better data visualization
    ax1.set_zlim(0.0 , 1.01)
    ax1.set_xlim(-0.01 , 1.1)
    ax1.set_ylim(-0.01 , 1.16)

    #Calculate system energy
    u0, u1, u2, u3 = c.Energy()
    xp = np.append(xp,t)
    yp1 = np.append(yp1,u0)
    yp2 = np.append(yp2,u1)
    yp3 = np.append(yp3,u2)
    yp4 = np.append(yp4,u3)
    #ypt = np.append(ypt,(u0+u1+u2+u3))
    #Plot energy
    pU0.set_data(xp,yp1)
    pU1.set_data(xp,yp2)
    pU2.set_data(xp,yp3)
    pU3.set_data(xp,yp4)
    #pUT.set_data(xp,ypt)
    
    
    #Update simulation
    for i in range(0,10):
        #c.simUpdateExplicit(0.0001,sms.explicit_method.rk4)
        c.semiImplictEuler(CONST_STEP)
        t = t+CONST_STEP
    """
    for i in range(0,3):
        c.ImplictEuler(0.0003)
    """
    return surf,pU0#, pEnergy

c = sms.Cloth(5,5,0.001,0.1)
#constr = np.arange(c.dY)
#constr = np.array([0,1,2,3,4,5,6, 7,8,9,10,11,12,13])
#constr = np.array([5,6,12,13,19,20,26,27,33,34,40,41,47,48])
constr = np.array([4,24])
c.constrain(constr) #Constrain specific particles

x = []
y = []
z = []
for p in c.X:
    x.append(p[0])
    y.append(p[1])
    z.append(p[2])
#Setup which points should be connected in the trimesh
triang = tri.Triangulation(x, z)


#Create figure
fig = plt.figure(figsize=(16,10))

ax1 = fig.add_subplot(131, projection='3d')
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(235)
ax5 = fig.add_subplot(236)

surf = ax1.plot(x,y)#plot_trisurf(x, y, z,color= 'b')
ax2.set_xlim(0.0, CONST_FRAMES*CONST_STEP*10)
ax3.set_xlim(0.0, CONST_FRAMES*CONST_STEP*10)
ax4.set_xlim(0.0, CONST_FRAMES*CONST_STEP*10)
ax5.set_xlim(0.0, CONST_FRAMES*CONST_STEP*10)

yp1 = np.zeros(0)
yp2 = np.zeros(0)
yp3 = np.zeros(0)
yp4 = np.zeros(0)

xp = np.zeros(0)
t = 0
pU0, = ax2.plot(xp,yp1,'b-')
pU1, = ax3.plot(xp,yp2,'b-')
pU2, = ax4.plot(xp,yp3,'b-')
pU3, = ax5.plot(xp,yp4,'b-')

ax2.set_ylim(0.0, 0.006)
ax3.set_ylim(-2e-3, 4.9e-04)
ax4.set_ylim(0.0, 10.0)
ax5.set_ylim(0.0, 0.4e-7)

#Title
ax2.set_title("Kinetic energy")
ax3.set_title("Gravitational energy")
ax4.set_title("Torsion spring energy")
ax5.set_title("Tension spring energy")
#axt.set_title("Total energy")

#Animate!
ani = animation.FuncAnimation(fig, data, fargs=(c, surf, triang),frames=CONST_FRAMES, interval=30, blit=False,repeat=False)
#ani.save(filename='sim.mp4',fps=30,dpi=300) #Save animation

plt.show()
