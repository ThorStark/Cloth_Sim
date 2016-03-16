"""
simPlot.py
Thor Stark Stenvang
Date of creation 17. feb 2016
"""
import springMassSystem as sms
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
import time
#Data plot update function
def data(i,c,surf,triang):
    #Obtain coordinates of points
    x = c.X[:,0]
    y = c.X[:,1]
    z = c.X[:,2]
    
    ax.clear() # Clear plot
    #Plot tri mesh
    surf = ax.plot_trisurf(x,y,z,triangles=triang.triangles,color= 'gray',linewidth=0.4)
    # Limit plot for better data visualization
    ax.set_zlim(0.5 , 1.2)
    ax.set_xlim(-0.1 , 3.1)
    ax.set_ylim(-0.1 , 3.1)
    
    #Update simulation
    for i in range(0,10):
        c.simUpdateExplicit(0.00010,sms.explicit_method.fe)
    """
    for i in range(0,3):
        c.ImplictEuler(0.0003)
    """
    return surf

c = sms.Cloth(7,7,0.005,0.5)
constr = np.array([0,1,2,3,4,5,6, 7,8,9,10,11,12,13])
c.constrain(constr) #Constrain specific particles

x = []
y = []
z = []
for p in c.X:
    x.append(p[0])
    y.append(p[1])
    z.append(p[2])
#Setup which points should be connected in the trimesh
triang = tri.Triangulation(x, y)

#Create figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_trisurf(x, y, z,color= 'b')
#Animate!
ani = animation.FuncAnimation(fig, data, fargs=(c, surf, triang), interval=30, blit=False)

plt.show()
