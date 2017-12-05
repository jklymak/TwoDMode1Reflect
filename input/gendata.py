import numpy as np
import scipy
import numpy.matlib as matlib
from shutil import copy
from os import mkdir
import os

def lininc(n,Dx,dx0):
  a=(Dx-n*dx0)*2./n/(n+1)
  print (a)
  dx = dx0+np.arange(1.,n+1.,1.)*a
  return dx

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
  print (outdir+' Exists')
try:
  mkdir(outdir+'/figs')
except:
  print (outdir+'/figs Exists')
copy('./gendata.py',outdir)

# These must match ../code/SIZE.h
ny = 1
nx = 16*128
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

with open(outdir+"/delXvar.bin", "wb") as f:
	dx.tofile(f)
f.close()

dy = np.zeros(ny) + dy
y=np.cumsum(dy)

with open(outdir+"/delYvar.bin", "wb") as f:
	dy.tofile(f)
f.close()

# topo
#slope of ray:

rayslope = np.sqrt( (N0**2 - om**2) / (om**2 - f0**2))



boodosa

sigma = 4000. # m


topo = 1500*np.exp(-x*x/(sigma**2))-1500+h0
#topo = h0*np.exp(-x*x/(3000**2))
print(np.shape(topo))
topo[topo<0.]=0.
topo=-H+topo
topo[topo<-H]=-H

# plot
TT0 = matlib.repmat(topo,1,ny)


with open(outdir+"/topo.bin", "wb") as f:
	TT0.tofile(f)
f.close()
# dz:
# dz is from the surface down (right?).  Its saved as positive.

dz=np.zeros(nz)+H/nz

with open(outdir+"/delZvar.bin", "wb") as f:
	dz.tofile(f)
f.close()

# temperature profile...
g=9.8
alpha = 2e-4
T0 = 28+np.cumsum(N0**2/g/alpha*(-dz))

with open(outdir+"/TRef.bin", "wb") as f:
	T0.tofile(f)
f.close()

# save T0 over whole domain
TT0 = matlib.repmat(T0,nx,ny).T
with open(outdir+"/T0.bin", "wb") as f:
	TT0.tofile(f)


z=np.cumsum(dz)
# plot:

# Forcing for boundaries
dt=3720.
time = np.arange(0,12.*3720.,dt)
print (time/3600./12.4)
om = 2*np.pi/12.40/3600;
uw = u0+0.*time
ue = u0+0.*time
# plot:

# try time,nz,ny...

uen=np.zeros((np.shape(time)[0],nz,ny))
for j in range(0,ny):
  for i in range(0,nz):
    uen[:,i,j]=ue
#print(uen)

uwn=np.zeros((np.shape(time)[0],nz,ny))
for j in range(0,ny):
  for i in range(0,nz):
    uwn[:,i,j]=uw
#print(uwn)

with open(outdir+"/Ue.bin","wb") as f:
  uen.tofile(f)

with open(outdir+"/Uw.bin", "wb") as f:
  uwn.tofile(f)

t=np.zeros((np.shape(time)[0],nz,ny))
for j in range(0,ny):
	for i in range(0,nz):
		for k in range(0,np.shape(time)[0]):
			t[k,i,j]=T0[i]
with open(outdir+"/Te.bin", "wb") as f:
	t.tofile(f)
f.close()
with open(outdir+"/Tw.bin", "wb") as f:
	t.tofile(f)
f.close()

## Copy some other files
import shutil
shutil.copy('data', outdir+'/data')
shutil.copy('eedata', outdir)
shutil.copy('data.kl10', outdir)
shutil.copy('data.mnc', outdir)
shutil.copy('data.obcs', outdir)
shutil.copy('data.nf90io', outdir)
shutil.copy('data.diagnostics', outdir)
shutil.copy('data.pkg', outdir+'/data.pkg')
# also store these.  They are small and helpful to document what we did
for nm in {'input','code','build_options','analysis'}:
    to_path = outdir+'/'+nm
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree('../'+nm, outdir+'/'+nm)
shutil.copy('../build/mitgcmuv', outdir)
