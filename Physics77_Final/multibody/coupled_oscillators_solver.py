"""
==============================================================
Multibody/Coupled Oscillators Solver
==============================================================
This module provides analytical solutions to the multibody
or coupled oscillators problem in physics. The results by using
such methods can be graphed, as demonstrated in the test folder.
There are methods that provide solutions in 1D and 2D and uses
a runge kutta method to do so.
-------
"""
__all__ = [
    'coupled_equation_solver_1D',
    'coupled_equation_solver_2D',
    'runge_kutta',
    'gravity',
    'lagrange1',
    'lagrange2',
    'lagrange3',
    'lagrange4',
    'lagrange5'
]


def coupled_equation_solver_1D(force, mass, initial_position, initial_velocity, dt=0.1, total_time=100,
                               integration_method='euler'):
    """
    This method numerically solves coupled 2nd order ordinary differential equations in 1 dimension. This method is specifically geared toward physics, hence the mass parameter.
    The parameters are as follows
      Let the solved equations be f(x),g(x) etc
      force: A list/tuple of functions corresponding to the differential equations. The first element corresponds to the first equation ie f''(x),g''(x) etc
      The function should be take parameters position, velocity, mass, time is that order, it should return a float. The inputs will be lists/tuples except time which is a double
      Ex: def force(position,velocity,mass,time):
            return  - position[0] - (position[0] - position[1])
      mass: A list/tuple of the masses of the objects. If there is no mass or it is not relevant, enter an array of 1
      initial_position: A list/tuple of the initial values of the solved function ie f(0),g(0) etc
      initial_velocity: A list/tuple of the initial values of the derivative of the solved function ie f'(0),g'(0) etc
      dt: The time step for the numerical solution. A smaller value will yield a more accurate result but require more computing time. Default value 0.1
      total_time: The length of the interval f''(x),g''(x) etc is to be solved over. Default value 100
      integration_method: The numerical solve method used to calculate the solution. 'euler' for euler-cromer method, 'taylor' for taylor series, 'stormer_verlet' for Stormer Verlet method, 'runge_kutta" for standard Runge Kutta method (RK4)
        Note: Euler-Cromer is the simplest and therefore is the fastest. Stormer Verlet is supposed to be more accurate than Euler-Cromer and comparable in speed. Runge-Kutta is the most accurate but the slowest
        Important Note!!!!: The current implementations of Taylor Series and Runge Kutta lead to uncontrolled growth, Do Not Use!

    Return Values:
      position: a list/tuple of lists/tuples of the values for the solved equation ie f(x),g(x) position[0] gives f(x), position[1] gives g(x) etc
      velocity: a list/tuple of lists/tuples of the values for the derivative of the solved equations ie f'(x),g'(x) velocity[0] gives f'(x), velocity[1] gives g'(x)
      errorP: a list/tuple of lists/tuples of th value of the estimated error of the position at each point using the Taylor series arror bound on O(t^2)
    """
    position = []  # Initialize lists for return
    velocity = []
    errorP = []
    for i in range(len(mass)):  ##Make sure the lists have sublists so the sublists can be appended
        tempx = []
        tempx.append(initial_position[i])
        position.append(tempx)
        tempv = []
        tempv.append(initial_velocity[i])
        velocity.append(tempv)
        tempe = [0]
        errorP.append(tempe)

    for i in range(int(total_time / dt) - 1):  ##Subtracting 1 from the range keeps the array length even with np.arange(0,total_time,dt)
        current_position = []  ##Get the current values to sent to the force function
        current_velocity = []
        current_errorP = []
        for j in range(len(mass)):
            current_position.append(position[j][i])
            current_velocity.append(velocity[j][i])
            current_errorP.append(errorP[j][i])
        for j in range(len(mass)):
            a = force[j](current_position, current_velocity, mass, i * dt) / mass[j]
            velocity[j].append(current_velocity[j] + a * dt)
            errorP[j].append(current_errorP[j] + np.abs(0.5 * a * dt * dt))
            if i == 0:
                position[j].append(current_position[j] + velocity[j][
                    i + 1] * dt)  ##Must initialize second value for Stormer Verlet to work
            elif integration_method == 'euler':
                position[j].append(current_position[j] + velocity[j][i + 1] * dt)
            elif integration_method == 'taylor':
                position[j].append(current_position[j] + dt * current_velocity[j] + 0.5 * a * dt * dt)
            elif integration_method == 'stormer_verlet':
                position[j].append(2 * current_position[j] - position[j][i - 1] + a * dt * dt)
            elif integration_method == 'runge_kutta':
                k = runge_kutta(force[j], mass, current_position, current_velocity, dt, j)
                position[j].append(current_position[j] + (1 / 6) * dt * (k[0] + 2 * k[1] + 2 * k[2] + k[3]))

    return position, velocity, errorP


def runge_kutta(force, mass, position, velocity, dt, i):
    initial_velocity = velocity[i]
    initial_position = position[i]
    a = force(position, velocity, mass, i * dt) / mass[i]
    k1 = initial_velocity
    position[i] = initial_position + k1 * dt / 2
    velocity[i] = initial_velocity + a * dt / 2
    a2 = force(position, velocity, mass, i * dt) / mass[i]
    k2 = initial_velocity + a2 * dt / 2
    position[i] = initial_position + k2 * dt / 2
    velocity[i] = initial_velocity + a2 * dt / 2
    a3 = force(position, velocity, mass, i * dt) / mass[i]
    k3 = initial_velocity + a3 * dt / 2
    position[i] = initial_position + k3 * dt / 2
    velocity[i] = initial_velocity + a3 * dt / 2
    a4 = force(position, velocity, mass, i * dt) / mass[i]
    k4 = initial_velocity + a4 * dt
    return k1, k2, k3, k4

def coupled_equation_solver_2D(force,mass,initial_position,initial_velocity,dt = 0.01,total_time = 100,integration_method = 'euler'):
  position = []
  velocity = []
  for i in range(len(mass)):
    tempx = []
    tempx.append(initial_position[i])
    position.append(tempx)
    tempv = []
    tempv.append(initial_velocity[i])
    velocity.append(tempv)
  for i in range(int(total_time/dt) - 1):
    current_position = []
    current_velocity = []
    for j in range(len(mass)):
      current_position.append(position[j][i])
      current_velocity.append(velocity[j][i])
    for j in range(len(mass)):
      ax,ay = force[j](current_position,current_velocity, mass)
      ax = ax/mass[j]
      ay = ay/mass[j]
      velocity[j].append([current_velocity[j][0] + ax*dt,current_velocity[j][1] + ay*dt])
      if i == 0:
        position[j].append([current_position[j][0] + velocity[j][i+1][0]*dt,current_position[j][1] + velocity[j][i+1][1]*dt])
      elif integration_method == 'euler':
        position[j].append([current_position[j][0] + velocity[j][i+1][0]*dt,current_position[j][1] + velocity[j][i+1][1]*dt])
      elif integration_method == 'taylor':
        position[j].append([current_position[j][0] + current_velocity[j][0]*dt + 0.5*ax*dt*dt,current_position[j][1] + current_velocity[j][1]*dt + 0.5*ay*dt*dt])
      elif integration_method == 'stormer_verlet':
        position[j].append([2*current_position[j][0] - position[j][i-1][0] + ax*dt*dt,2*current_position[j][1] - position[j][i-1][1] + ay*dt*dt])

  return position,velocity


import numpy as np
import matplotlib.pyplot as plt


##Implementation of 2D equation solver


def gravity(mass, initial_position, initial_velocity, dt=0.01, total_time=100, integration_method='euler'):
    def force1(position, velocty, mass):
        f = -mass[0] * mass[1] / ((position[0][0] - position[1][0]) * (position[0][0] - position[1][0]) + (
                    position[0][1] - position[1][1]) * (position[0][1] - position[1][1]))
        theta = np.arctan2(position[0][1] - position[1][1], position[0][0] - position[1][0])
        return f * np.cos(theta), f * np.sin(theta)

    def force2(position, velocity, mass):
        f = -mass[0] * mass[1] / ((position[0][0] - position[1][0]) * (position[0][0] - position[1][0]) + (
                    position[0][1] - position[1][1]) * (position[0][1] - position[1][1]))
        theta = np.arctan2(position[1][1] - position[0][1], position[1][0] - position[0][0])
        return f * np.cos(theta), f * np.sin(theta)

    force = [force1, force2]
    position, velocity = coupled_equation_solver_2D(force, mass, initial_position, initial_velocity, dt, total_time,
                                                    integration_method)
    x1 = []
    x2 = []
    y1 = []
    y2 = []
    for i in range(len(position[0])):
        x1.append(position[0][i][0])
        x2.append(position[1][i][0])
        y1.append(position[0][i][1])
        y2.append(position[1][i][1])
    coordinates = [[x1, y1], [x2, y2]]
    return coordinates


def lagrange1(mass, initial_position, initial_velocity, dt=0.01, total_time=100, integration_method='euler'):
    position = gravity(mass, initial_position, initial_velocity, dt, total_time, integration_method)

    lagrange = []
    for i in range(len(position[0][0])):
        R = np.sqrt((position[1][0][i] - position[0][0][i]) ** 2 + (position[1][1][i] - position[0][1][i]) ** 2)
        theta = np.arctan2((position[1][1][i] - position[0][1][i]), (position[1][0][i] - position[0][0][i]))
        r = R * (mass[0] / (3 * mass[1])) ** (1 / 3)
        R1 = R - r
        lagrange.append([R1 * np.cos(theta), R1 * np.sin(theta)])
    return lagrange


def lagrange2(mass, initial_position, initial_velocity, dt=0.01, total_time=100, integration_method='euler'):
    position = gravity(mass, initial_position, initial_velocity, dt, total_time, integration_method)

    lagrange = []
    for i in range(len(position[0][0])):
        R = np.sqrt((position[1][0][i] - position[0][0][i]) ** 2 + (position[1][1][i] - position[0][1][i]) ** 2)
        theta = np.arctan2((position[1][1][i] - position[0][1][i]), (position[1][0][i] - position[0][0][i]))
        r = R * (mass[0] / (3 * mass[1])) ** (1 / 3)
        R1 = R + r
        lagrange.append([R1 * np.cos(theta), R1 * np.sin(theta)])
    return lagrange


def lagrange3(mass, initial_position, initial_velocity, dt=0.01, total_time=100, integration_method='euler'):
    position = gravity(mass, initial_position, initial_velocity, dt, total_time, integration_method)

    lagrange = []
    for i in range(len(position[0][0])):
        R = np.sqrt((position[1][0][i] - position[0][0][i]) ** 2 + (position[1][1][i] - position[0][1][i]) ** 2)
        theta = np.arctan2((position[1][1][i] - position[0][1][i]), (position[1][0][i] - position[0][0][i]))
        R1 = -R
        lagrange.append([R1 * np.cos(theta), R1 * np.sin(theta)])
    return lagrange


def lagrange4(mass, initial_position, initial_velocity, dt=0.01, total_time=100, integration_method='euler'):
    position = gravity(mass, initial_position, initial_velocity, dt, total_time, integration_method)

    lagrange = []
    for i in range(len(position[0][0])):
        R = np.sqrt((position[1][0][i] - position[0][0][i]) ** 2 + (position[1][1][i] - position[0][1][i]) ** 2)
        theta = np.arctan2((position[1][1][i] - position[0][1][i]), (position[1][0][i] - position[0][0][i]))
        theta1 = theta + np.pi / 3
        lagrange.append([R * np.cos(theta1), R * np.sin(theta1)])
    return lagrange


def lagrange5(mass, initial_position, initial_velocity, dt=0.01, total_time=100, integration_method='euler'):
    position = gravity(mass, initial_position, initial_velocity, dt, total_time, integration_method)

    lagrange = []
    for i in range(len(position[0][0])):
        R = np.sqrt((position[1][0][i] - position[0][0][i]) ** 2 + (position[1][1][i] - position[0][1][i]) ** 2)
        theta = np.arctan2((position[1][1][i] - position[0][1][i]), (position[1][0][i] - position[0][0][i]))
        theta1 = theta - np.pi / 3
        lagrange.append([R * np.cos(theta1), R * np.sin(theta1)])
    return lagrange

