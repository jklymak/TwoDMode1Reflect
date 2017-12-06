import numpy as np
import scipy
import numpy.matlib as matlib
import matplotlib.pyplot as plt
from shutil import copy, move
import shutil
from os import mkdir

import datetime
import os

def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  print (a)
  dx = dx0+np.arange(1.,n+1.,1.)*a
  return dx

fig, ax = plt.subplots()

Fr=0.13
H = 2000.
om = 2.*np.pi/12.42/3600.
f0 = 1.e-4
N0=5.2e-3
u0=0.02


outdir='../results/RunN%02dU%02dNew' % (N0*1e4, u0*1e2)
try:
  mkdir(outdir)
except:
  print (outdir+' Exists.  Backing up')
  move(outdir, outdir+'.bak%s' % datetime.datetime.utcnow().isoformat())
  mkdir(outdir)

for nm in {'input','code','build_options','analysis'}:
    to_path = outdir+'/'+nm
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree('../'+nm, outdir+'/'+nm)

mkdir(outdir + '/build')
try:
    shutil.copy('../build/mitgcmuv', outdir+'/build')
except:
    print('Warning: no mitgcmuv copied')

try:
  mkdir(outdir+'/figs')
except:
  print (outdir+'/figs Exists')
try:
  mkdir(outdir+'/indata')
except:
  print (outdir+'/indata Exists')
copy('./gendata.py',outdir+'/input/')

# These must match ../code/SIZE.h
ny = 1
nx = 16*128
print(nx)
nz = 200


# y direction:
dy = 1000
# x direction
xt = 410e3

# slope will be about 100-km wide:
nleft = 1000
dx0=100.
dx = np.zeros(nx) + dx0
for n in range(nleft, nx):
    dx[n] = dx[n - 1] * 1.003
x=np.cumsum(dx)
print('max x: %f [km]' % (x[-1] / 1000.))

with open(outdir+"/indata/delXvar.bin", "wb") as f:
	dx.tofile(f)
f.close()

dy = np.zeros(ny) + dy
y=np.cumsum(dy)

with open(outdir+"/indata/delYvar.bin", "wb") as f:
	dy.tofile(f)
f.close()

# topo
#slope of ray:

rayslope = np.sqrt( (N0**2 - om**2) / (om**2 - f0**2))
print(rayslope)
# shelf = rayslope/4
dz0 = 100. / rayslope / 10
topo = np.arange(0, 200, dz0)
dz0 = 100. / rayslope
topo = np.hstack((topo, np.arange(topo[-1], 400, dz0)))
dz0 = 100. / rayslope * 5
topo = np.hstack((topo, np.arange(topo[-1], H, dz0)))
topo = np.hstack((topo, H * (np.ones(nx-len(topo)))))
topo = np.convolve(topo, np.ones(50)/50., mode='same')
topo = -topo
topo[topo<-H]=-H

# plot
fig, ax = plt.subplots()
ax.plot(x/1e3, topo)
ax.set_ylim([-H, 0])
ax.set_xlim([0., 400.])
fig.savefig(outdir+'/figs/Topo.png')

TT0 = matlib.repmat(topo,1,ny)

with open(outdir+"/topo.bin", "wb") as f:
	TT0.tofile(f)
f.close()
# dz:
# dz is from the surface down (right?).  Its saved as positive.

dz=np.zeros(nz)+H/nz

with open(outdir+"/indata/delZvar.bin", "wb") as f:
	dz.tofile(f)
f.close()

# temperature profile...
g=9.8
alpha = 2e-4
T0 = 28+np.cumsum(N0**2/g/alpha*(-dz))

with open(outdir+"/indata/TRef.bin", "wb") as f:
	T0.tofile(f)
f.close()

# save T0 over whole domain
TT0 = matlib.repmat(T0,nx,ny).T
with open(outdir+"/indata/T0.bin", "wb") as f:
	TT0.tofile(f)


z=np.cumsum(dz)
# plot:

# plot:
fig, ax = plt.subplots()
ax.plot(T0, z)
ax.set_ylim(H, 0)
fig.savefig(outdir+'/figs/T0.png')


###############################################
# Forcing for boundaries

# get central z...
zc = np.cumsum(dz) - dz / 2.

dt=3720./2.
time = np.arange(0,24.*3720.,dt)

dzi = 100
N2i = N0**2 * np.ones(19)

import vertmodes

psi, phi, ce, zpsi =vertmodes.vertModes(N2i,dzi)

fig, ax = plt.subplots()
if ce[0] < 10:
    ind = 0
else:
    ind = 1
psi = psi[:, ind]
if psi[0] < 0:
    psi = -psi
ce = ce[ind]
dpsi = np.diff(psi)/dzi
# interp onto our grid...
psi = np.interp(zc, zpsi, psi)
dpsi = np.interp(zc, zpsi[:-1] - dzi/2, dpsi)

# horizontal wavelength...
kr = np.sqrt((om**2 - f0**2) / ce**2)
print('lambda_r=%f' % (2. * np.pi / kr))
cg = np.sqrt(om**2 - f0**2) / om * ce
print('group speed=%f' % cg)
print(ce)
ax.plot(psi, zc)
fig.savefig(outdir + '/figs/psi.png')

ntosave = 24
externalForcingPeriod=12.4*3600/ntosave
dt = externalForcingPeriod
externalForcingCycle=12.4*3600
print("externalForcingPeriod="+'%1.0f'%externalForcingPeriod)
print("externalForcingCycle="+'%1.0f'%externalForcingCycle)


kx = -kr
ky = 0.


# a bit easier to do in p:
p0 = np.abs ( (om**2 - f0**2) / (kx * om + 1j * ky * f0)) * u0
print(p0)
print(p0 * (kx * om + 1j * ky *f0) / (om**2 - f0**2))
Amp = p0 * np.exp(1j * (kx * (x-x[-1]) + ky * y))
Amp = Amp / np.max(psi)

fig, ax = plt.subplots()
ax.plot(x/1.e3, np.real(Amp))
ax.plot(x/1.e3, np.real(Amp), 'x')
fig.savefig(outdir + '/figs/px.png')

U = np.outer(psi.T,
            Amp * (kx * om + 1j * ky *f0) / (om**2 - f0**2))

fig, ax = plt.subplots()
pcm = ax.pcolormesh(x/1.e3, z, np.real(U), rasterized=True, vmin=-u0, vmax=u0)
fig.colorbar(pcm)
fig.savefig(outdir + '/figs/U0.png')

V = np.outer(psi.T,
            Amp * (ky * om - 1j * kx *f0) / (om**2 - f0**2))

fig, ax = plt.subplots()
pcm = ax.pcolormesh(x/1.e3, z, np.real(V), rasterized=True, vmin=-u0, vmax=u0)
fig.colorbar(pcm)
fig.savefig(outdir + '/figs/V0.png')
grav = 9.8
alpha = 2e-4
Tp = np.outer(dpsi,  (-Amp) / grav / alpha)

fig, ax = plt.subplots()
pcm = ax.pcolormesh(x/1.e3, z, np.real(Tp), rasterized=True)
fig.colorbar(pcm)
fig.savefig(outdir + '/figs/Tp0.png')

Tf = Tp + T0[:, np.newaxis]

print("Writing forcing files:spongeweight")

## Make the sponge indicator
spongew = len(np.where(x > (x[-1] - 300.e3))[0])
print('spongew: %d' % spongew)
aa = np.zeros((nz, nx))
aa[:,-spongew:] = np.linspace(0.0,1.,spongew)

with open(outdir+"/indata/spongeweight.bin", "wb") as f:
    aa.tofile(f)

###########################

print("Writing forcing files:")

with open(outdir+"/indata/Uforce.bin", "wb") as f:
    for tind in range(0,ntosave):
        t = tind*dt+dt/2.
        aa = np.real(U * np.exp(-1j*t*om))
        aa.tofile(f)

with open(outdir+"/indata/Vforce.bin", "wb") as f:
    for tind in range(0,ntosave):
        t = tind*dt+dt/2.
        aa = np.real(V * np.exp(-1j*t*om))
        aa.tofile(f)

with open(outdir+"/indata/Tforce.bin", "wb") as f:
    for tind in range(0,ntosave):
        t = tind*dt+dt/2.
        aa = np.real(Tf * np.exp(-1j*t*om))
        aa.tofile(f)



## Copy some other files
import shutil
if 0:
    shutil.copy('data', outdir+'input/data')
    shutil.copy('eedata', outdir+'input/')
    shutil.copy('data.kl10', outdir+'input/')
    shutil.copy('data.mnc', outdir+'input/')
    shutil.copy('data.obcs', outdir+'input/')
    shutil.copy('data.nf90io', outdir+'input/')
    shutil.copy('data.diagnostics', outdir+'input/')
    shutil.copy('data.pkg', outdir+'input/data.pkg')
# also store these.  They are small and helpful to document what we did
