import matplotlib.pyplot as plt
import numpy as np
import multibody.coupled_oscillators_solver as multi

mass = [10,6]
initial_position = [[0,0],[8,0]]
initial_velocity = [[0,-0],[0,1]]
dt = 0.01
total_time = 100


positions = multi.gravity(mass,initial_position,initial_velocity)
l1 = multi.lagrange1(mass,initial_position,initial_velocity)
l2 = multi.lagrange2(mass,initial_position,initial_velocity)
l3 = multi.lagrange3(mass,initial_position,initial_velocity)
l4 = multi.lagrange4(mass,initial_position,initial_velocity)
l5 = multi.lagrange5(mass,initial_position,initial_velocity)



i = 1

xbig = positions[0][0][i]
ybig = positions[0][1][i]
xsmall = positions[1][0][i]
ysmall = positions[1][1][i]
x1 = l1[i][0]
y1 = l1[i][1]
x2 = l2[i][0]
y2 = l2[i][1]
x3 = l3[i][0]
y3 = l3[i][1]
x4 = l4[i][0]
y4 = l4[i][1]
x5 = l5[i][0]
y5 = l5[i][1]


fig = plt.figure()
ax = fig.add_subplot(111)


plt.plot(xbig,ybig,'k.')
plt.plot(xsmall,ysmall,'r.')
plt.plot(x1,y1,'y.')
plt.plot(x2,y2,'y.')
plt.plot(x3,y3,'y.')
plt.plot(x4,y4,'y.')
plt.plot(x5,y5,'y.')
plt.xlim(-20,20)
plt.ylim(-20,20)
ax.set_aspect('equal', adjustable='box')
plt.show()






i = 200

xbig = positions[0][0][i]
ybig = positions[0][1][i]
xsmall = positions[1][0][i]
ysmall = positions[1][1][i]
x1 = l1[i][0]
y1 = l1[i][1]
x2 = l2[i][0]
y2 = l2[i][1]
x3 = l3[i][0]
y3 = l3[i][1]
x4 = l4[i][0]
y4 = l4[i][1]
x5 = l5[i][0]
y5 = l5[i][1]


fig = plt.figure()
ax = fig.add_subplot(111)


plt.plot(xbig,ybig,'k.')
plt.plot(xsmall,ysmall,'r.')
plt.plot(x1,y1,'y.')
plt.plot(x2,y2,'y.')
plt.plot(x3,y3,'y.')
plt.plot(x4,y4,'y.')
plt.plot(x5,y5,'y.')
plt.xlim(-20,20)
plt.ylim(-20,20)
ax.set_aspect('equal', adjustable='box')
plt.show()

