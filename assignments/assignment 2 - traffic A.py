"""assignment 3 of the Numerical MOOC.
Traffic flow
1D model solved by forward difference in time and backward difference in space 
part A
"""
import numpy as np
import matplotlib.pyplot as plt

# problem parameters
V_max = 80. # km/u
L = 11. # km
rho_max = 250. # cars/km
nx = 51
dx = L/(nx-1) # km   
dt = 0.001 # hours
nt = int(6./60./dt)
cV = 1./3.6

# initial conditions
x = np.linspace(0,L,nx)
rho0 = np.ones(nx)*10.
rho0[10:20] = 50.

def dFdrho(rho):
    return V_max*(1. - 2.*rho/rho_max)

def V(rho):
    return V_max*(1. - rho/rho_max)

# time integration
rho = rho0.copy()
for n in range(1, nt):  
    rhon = rho.copy()
    dFdrhon = dFdrho(rhon)
    rho[1:] = rhon[1:] - dFdrhon[1:]*dt/dx*(rhon[1:] - rhon[0:-1]) 
    rho[0] = 10.
    if n == 49:
        rho_3min = rho.copy()

# questions
print("V0_min", cV*np.min(V(rho0)))
#~ print("v_mean", cV*np.sum(V(rho_3min)*rho_3min)/np.sum(rho_3min))
print("v_mean alt", cV*np.mean(V(rho_3min)))
print("v_min", cV*np.min(V(rho)))

# plot of initial and final solutions
plt.plot(x, rho0, 'b:', label='rho0')
plt.plot(x, rho, 'b', label='rho')
plt.plot(x, V(rho0), 'g:', label='V0')
plt.plot(x, V(rho), 'g', label='V')
plt.legend(loc=0)
plt.ylim(0., V_max + 5.);
plt.show()


