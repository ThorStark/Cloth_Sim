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
    ax.set_zlim(-10.2 , 1.2)
    ax.set_xlim(-0.2 , 9.2)
    ax.set_ylim(-0.2 , 9.2)
    
    #Update simulation
    for i in range(0,20):
        c.simUpdateExplicit(0.00014,sms.explicit_method.runge_kutta_4)
    return surf

c = sms.Cloth(10,10,0.001)
constr = np.array([0,9,90,99])
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
