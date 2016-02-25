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
def data(i,c,surf,triang):
    x = []
    y = []
    z = []
    for p in c.X:
        x.append(p[0])
        y.append(p[1])
        z.append(p[2])
        
    ax.clear()

    surf = ax.plot_trisurf(x,y,z,triangles=triang.triangles,color= 'gray',linewidth=0.4)
    ax.set_zlim(-6.0 , 1.2)
    ax.set_xlim(-1.0 , 10)
    ax.set_ylim(-1.0 , 10)
    for i in range(0,10):
        c.simUpdateExplicit(0.0005,sms.explicit_method.forward_euler)
    return surf

c = sms.Cloth(10,10,0.29)
x = []
y = []
z = []
for p in c.X:
    x.append(p[0])
    y.append(p[1])
    z.append(p[2])

triang = tri.Triangulation(x, y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlim(0.0 , 1.2)

surf = ax.plot_trisurf(x, y, z,color= 'b')
ani = animation.FuncAnimation(fig, data, fargs=(c, surf, triang), interval=30, blit=False)

plt.show()
