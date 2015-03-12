"""Sod's shock tube
numerical mooc assignment of module 3
"""

import numpy as np
import matplotlib.pyplot as plt

nx = 81
dx = .25
T = 0.01
dt = .0002
gamma = 1.4
nt = int(T/dt)
x = np.linspace(-10., 10., nx)

# initial values
rho0 = np.ones(nx)
rho0[nx/2:] = 0.125
v0 = np.zeros(nx)
p0 = 100e3*np.ones(nx)
p0[nx/2:] = 10e3
#~ print(x[40])
#~ print(p0[40])
u0 = np.zeros((nx,3))
u0[:,0] = rho0
u0[:,1] = rho0*v0
u0[:,2] = p0/(gamma - 1.) + 0.5*rho0*v0**2


def computeF(u,gamma):
    """return flux"""
    F = np.zeros_like(u)
    u1 = u[:,0]
    u2 = u[:,1]
    u3 = u[:,2]
    uf = u2**2/u1
    F[:,0] = u2
    F[:,1] = uf + (gamma - 1.)*(u3 - 0.5*uf)
    F[:,2] = (u3 + (gamma - 1.)*(u3 - 0.5*uf))*u2/u1
    return F

def richtmyer(u, nt, dt, dx, gamma):
    """Compute the solution with the Richtmyer scheme"""    
    for t in range(nt):
        u_n = u.copy()
        f = computeF(u_n, gamma)
        u_star = 0.5*(u[1:,:] + u[:-1,:]) - 0.5*dt/dx*(f[1:,:] - f[:-1,:])
        f_star = computeF(u_star, gamma)
        u[1:-1,:] = u_n[1:-1,:] - dt/dx*(f_star[1:,:] - f_star[:-1,:])
    return u

u = richtmyer(u0, nt, dt, dx, gamma)

# calculate parameters from solution u 
rho = u[:,0]
v = u[:,1]/rho
p = (gamma - 1.)*(u[:,2] - 0.5*rho*v**2)

print("x", x[50])
print("v", v[50])
print("p", p[50])
print("rho", rho[50])

# show graphs
plt.figure()
plt.plot(x, rho)
plt.plot(x, rho0)
#~ plt.ylabel('rho')
plt.grid()
plt.figure()
plt.plot(x, v)
#~ plt.plot(x, v0)
plt.ylabel('v')
plt.grid()
plt.figure()
plt.plot(x, p)
#~ plt.plot(x, p0)
plt.ylabel('p')
plt.grid()
plt.show()