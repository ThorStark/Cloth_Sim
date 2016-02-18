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
import time

def data(i,c,line):
    #time.sleep(5.0)
    x = []
    y = []
    z = []
    for p in c.particles:
        point = p.pos3D_
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])
    ax.clear()
    for s in c.springs:
        p0 = c.particles[s.indI_].pos3D_
        p1 = c.particles[s.indJ_].pos3D_
        line = ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]],color= 'g')
    line = ax.scatter(x, y, z,color= 'b')
    ax.set_zlim(-2.0 , 1.2)
    for i in range(0,20):
        c.simUpdate1(0.001)
    return line

c = sms.Cloth(10,10,0.1)
c.constrainParticle(0)
c.constrainParticle(9)
c.constrainParticle(90)
c.constrainParticle(99)
x = []
y = []
z = []
for p in c.particles:
    point = p.pos3D_
    x.append(point[0])
    y.append(point[1])
    z.append(point[2])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_zlim(0.0 , 1.2)

line = ax.scatter(x, y, z,color= 'b')
for s in c.springs:
    p0 = c.particles[s.indI_].pos3D_
    p1 = c.particles[s.indJ_].pos3D_
    line = ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]],color= 'g')

ani = animation.FuncAnimation(fig, data, fargs=(c, line), interval=30, blit=False)

plt.show()

"""
TEST
test = fabric(2,3,10.0,13.0)
test.initParticles()
a = test.getArea()
print(a)
"""
