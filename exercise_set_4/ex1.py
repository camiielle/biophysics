#  Problem 1: Explore the steady states of the system computationally
# Set A = 1. 
# Discretize time and space. 
# Simulate the particles position as a function of time. 
# Histogram the position and plot it together with the formula for the steady-state probability distribution derived in class.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

A=1.
# discretizing time
dt=0.001
duration=1000
iterations=int(duration/dt)
#boundaries
lower=-1
upper=1

def V(x):
    return np.cos(1.5*np.pi*x)

def dV_dx(x):
    return -np.sin(1.5*np.pi*x)*1.5*np.pi

def reflective_boundary(x):
    if x>upper:
        x=2*upper-x
    if x<lower:
        x=2*lower-x
    return x

def x_new(x):
    return x - dV_dx(x)*dt + np.sqrt(2*dt)*np.random.normal()

#simulation
x_initial=np.random.uniform(-1,1)
positions=[]
for i in range(iterations):
    x = x_new(x_initial)
    x = reflective_boundary(x)
    positions.append(x)
    x_initial=x

assert(len(positions)==iterations)

#steady state prob distribution saw in class
def steady_P(x):
    return np.exp(-V(x))

normalization, _ = quad(steady_P, -1, 1)

#plotting
plt.figure()
x_array = np.linspace(-1, 1, 1000)
plt.plot(x_array, steady_P(x_array) / normalization, color='r',label='Theoretical Distribution')
hist, _ = np.histogram(positions, bins=150)
plt.hist(positions, bins = 150, edgecolor = 'blue', alpha =0.5, density=True, label = 'Numerically Simulated Distribution')
plt.xlabel('x')
plt.ylabel('Probability Distribution')
plt.xlim(-1, 1)
plt.legend()
plt.show()