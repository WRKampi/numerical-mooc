"""assignment 4 of the Numerical MOOC.
Gray-Scott model of Reaction Diffusion.
2D model solved by forward time central space 
finite difference scheme.
"""
import numpy
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
import matplotlib.cm as cm
from matplotlib import animation

Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.065 ; name = 'Bacteria 1'
#~ Du, Dv, F, k = 0.00014, 0.00006, 0.035, 0.065 ; name = 'Bacteria 2'
#~ Du, Dv, F, k = 0.00016, 0.00008, 0.060, 0.062 ; name = 'Coral'
#~ Du, Dv, F, k = 0.00019, 0.00005, 0.060, 0.062 ; name = 'Fingerprint'
#~ Du, Dv, F, k = 0.00010, 0.00010, 0.018, 0.050 ; name = 'Spirals'
#~ Du, Dv, F, k = 0.00012, 0.00008, 0.020, 0.050 ; name = 'Spirals Dense'
#~ Du, Dv, F, k = 0.00010, 0.00016, 0.020, 0.050 ; name = 'Spirals Fast'
#~ Du, Dv, F, k = 0.00016, 0.00008, 0.020, 0.055 ; name = 'Unstable'
#~ Du, Dv, F, k = 0.00016, 0.00008, 0.050, 0.065 ; name = 'Worms 1'
#~ Du, Dv, F, k = 0.00016, 0.00008, 0.054, 0.063 ; name = 'Worms 2'
#~ Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.060 ; name = 'Zebrafish'

n = 192
dh = 5./(n-1)
T = 8000
dt = .9 * dh**2 / (4*max(Du,Dv))
nt = int(T/dt)

def initiate():
    """read initial values from file"""
    uvinitial = numpy.load('./uvinitial.npz')
    U = uvinitial['U']
    V = uvinitial['V']
    return U, V

def ftcs(U, V, nt, Du, Dv, F, k, dt, dx, dy):
    """2D forward time central space finite differene solver
    for Gray scott model of Reaction Diffusion with Neumann BCs"""
    for n in range(nt):
        Un = U.copy()
        Vn = V.copy()        
        U[1:-1,1:-1] = Un[1:-1,1:-1] + Du*\
            (dt/dy**2 * (Un[2:,1:-1] - 2*Un[1:-1,1:-1] + Un[:-2,1:-1]) +\
             dt/dx**2 * (Un[1:-1,2:] - 2*Un[1:-1,1:-1] + Un[1:-1,:-2])) + dt*(\
             -Un[1:-1,1:-1]*Vn[1:-1,1:-1]**2 + F*(1 - Un[1:-1,1:-1]) )
        V[1:-1,1:-1] = Vn[1:-1,1:-1] + Dv*\
            (dt/dy**2 * (Vn[2:,1:-1] - 2*Vn[1:-1,1:-1] + Vn[:-2,1:-1]) +\
             dt/dx**2 * (Vn[1:-1,2:] - 2*Vn[1:-1,1:-1] + Vn[1:-1,:-2])) + dt*(\
             Un[1:-1,1:-1]*Vn[1:-1,1:-1]**2 - (F + k)*Vn[1:-1,1:-1] )             
        # Enforce Neumann BCs
        U[-1,:] = U[-2,:] ; V[-1,:] = V[-2,:]
        U[0,:] = U[1,:]   ; V[0,:] = V[1,:]
        U[:,-1] = U[:,-2] ; V[:,-1] = V[:,-2]
        U[:,0] = U[:,1]   ; V[:,0] = V[:,1]         
    return U, V

## output for assignment questions:

U, V = initiate()
ftcs(U, V, nt, Du, Dv, F, k, dt, dh, dh)
print(U[100,::40])

## output movies:

dnt = 100       # number of time increments per frame
nF = nt//dnt    # number of frames

U, V = initiate()
fig = pyplot.figure(figsize=(8,8))
im = pyplot.imshow(U, cmap = cm.RdBu)
pyplot.xticks([]), pyplot.yticks([]);

def init():
    """return initial frame"""
    im.set_array(U)
    return im,

def animate(i):
    """return a frame each dnt time increments.
    (global variables are used)
    """
    print('frame {0} of {1}'.format(i, nF))
    ftcs(U, V, dnt, Du, Dv, F, k, dt, dh, dh)
    im.set_array(U)
    return im,

metadata = dict(title=name, artist='Matplotlib', comment='Numerical MOOC')
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=metadata)
im_ani = animation.FuncAnimation(fig, animate, init_func=init, frames=nF, blit=True)
im_ani.save(name+'.mp4', writer=writer)