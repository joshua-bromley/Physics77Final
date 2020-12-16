import numpy as np
import matplotlib.pyplot as plt
import multibody.coupled_oscillators_solver as multi
##Implementation of 2D equation solver

def force1(position,velocty,mass):
  f = -mass[0]*mass[1]/((position[0][0] - position[1][0])*(position[0][0] - position[1][0]) + (position[0][1]-position[1][1])*(position[0][1]-position[1][1]))
  theta = np.arctan2(position[0][1]-position[1][1],position[0][0] - position[1][0])
  return f*np.cos(theta),f*np.sin(theta)

def force2(position,velocity,mass):
  f = -mass[0]*mass[1]/((position[0][0] - position[1][0])*(position[0][0] - position[1][0]) + (position[0][1]-position[1][1])*(position[0][1]-position[1][1]))
  theta = np.arctan2(position[1][1]-position[0][1],position[1][0] - position[0][0])
  return f*np.cos(theta),f*np.sin(theta)

force = [force1,force2]
mass = [4,1]
initial_position = [[0,0],[0.3,0]]
initial_velocity = [[0,-1],[0,4]]
dt = 0.005
total_time = 50


position,velocity = multi.coupled_equation_solver_2D(force,mass,initial_position,initial_velocity,dt,total_time,'euler')
time = np.arange(0,total_time,dt)

x1 = []
x2 = []
y1 = []
y2 = []
for i in range(len(position[0])):
  x1.append(position[0][i][0])
  y1.append(position[0][i][1])
  x2.append(position[1][i][0])
  y2.append(position[1][i][1])

plt.title("Euler")
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.show()

position,velocity = multi.coupled_equation_solver_2D(force,mass,initial_position,initial_velocity,dt,total_time,'taylor')
time = np.arange(0,total_time,dt)

x1 = []
x2 = []
y1 = []
y2 = []
for i in range(len(position[0])):
  x1.append(position[0][i][0])
  y1.append(position[0][i][1])
  x2.append(position[1][i][0])
  y2.append(position[1][i][1])

plt.title("Taylor Series")
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.show()

position,velocity = multi.coupled_equation_solver_2D(force,mass,initial_position,initial_velocity,dt,total_time,'stormer_verlet')
time = np.arange(0,total_time,dt)

x1 = []
x2 = []
y1 = []
y2 = []
for i in range(len(position[0])):
  x1.append(position[0][i][0])
  y1.append(position[0][i][1])
  x2.append(position[1][i][0])
  y2.append(position[1][i][1])

plt.title("Stormer Verlet")
plt.plot(x1,y1)
plt.plot(x2,y2)
plt.show()




import numpy as np
import matplotlib.pyplot as plt

##Implementation of the 1D equation solver
def force3(position,velocity,mass, time):
  return  -(position[0] - np.sin(np.pi*time)) - (position[0] - position[1])

def force4(position,velocity,mass, time):
  return -position[1] - (position[1]-position[0])

force = [force3,force4]
mass = [1,1]
initial_position = [1,0]
initial_velocity = [0,0]
dt = 0.01
total_time = 100


position,velocity,error = multi.coupled_equation_solver_1D(force,mass,initial_position,initial_velocity,dt,total_time,'euler')

time = np.arange(0,total_time,dt)

plt.title("Euler")
plt.plot(time,position[0])
plt.plot(time,error[0])
plt.show()

position,velocity,error = multi.coupled_equation_solver_1D(force,mass,initial_position,initial_velocity,dt,total_time,'taylor')

time = np.arange(0,total_time,dt)

plt.title("Taylor Series")
plt.plot(time,position[0])
plt.plot(time,error[0])
plt.show()

position,velocity,error = multi.coupled_equation_solver_1D(force,mass,initial_position,initial_velocity,dt,total_time,'stormer_verlet')

time = np.arange(0,total_time,dt)

plt.title("Stormer Verlet")
plt.plot(time,position[0])
plt.plot(time,error[0])
plt.show()

position,velocity,error = multi.coupled_equation_solver_1D(force,mass,initial_position,initial_velocity,dt,total_time,'runge_kutta')

time = np.arange(0,total_time,dt)

plt.title("Runge Kutta")
plt.plot(time,position[0])
plt.plot(time,error[0])
plt.show()

