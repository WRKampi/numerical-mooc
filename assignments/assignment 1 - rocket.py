"""assignment 1 of the Numerical MOOC.
Vertical flight of rocket.
Forward Euler integration.
"""
import numpy as np
import matplotlib.pyplot as plt

def vmp(t):
    """return propellant burn rate as function of time"""
    if t <= 5.:
        return 20.
    else:
        return 0.

def dydt(h, v, mp, vmp, g, vc, rho, A, CD, ms):
    """return time derivative of h, v and mp"""
    dhdt = v
    dvdt = -g + (vmp*vc - 0.5*rho*v*abs(v)*A*CD)/(ms + mp)
    dmpdt = -vmp
    return dhdt, dvdt, dmpdt

# parameters
r = 0.5
params = {
    'ms' : 50.,
    'g' : 9.81,
    'rho' : 1.091,
    'A' : np.pi*r**2,
    'vc' : 325.,
    'CD' : 0.15,
    }

# time steps:
dt = 0.1
times = np.arange(0., 40., dt)

# solution arrays
hs = np.zeros_like(times)
mps = np.zeros_like(times)
vs = np.zeros_like(times)

# initial values
h0 = 0.
v0 = 0.
mp0 = 100.
hs[0] = h = h0
vs[0] = v = v0
mps[0] = mp = mp0

# time integration
for i, t in enumerate(times[1:]):
    dhdt, dvdt, dmpdt = dydt(h, v, mp, vmp(t), **params)
    h  += dt*dhdt
    v  += dt*dvdt
    mp += dt*dmpdt
    hs[i+1] = h
    vs[i+1] = v
    mps[i+1] = mp

# questions:
i = np.argwhere(times==3.2)[0][0]
print('propellant mass', mps[i])
i =  np.argmax(vs)
print('max speed', vs[i])
print('time of max speed', times[i])
print('altitude at max speed', hs[i])
i = np.argmax(hs)
print('max altitude', hs[i])
print('time of max altitude', times[i])
i = np.argmin(np.abs(hs[10:])) + 10
a = hs[i] - hs[i-1]
b = hs[i-1]
x = -b/a # linear interploation
print('impact time', times[i-1] + dt*x)
print('impact velocity', vs[i-1] + x*(vs[i] - vs[i-1]))
#~ print(vs[i])

# plot of trajectory
plt.plot(times, hs)
plt.grid()
plt.tight_layout()
plt.show()
