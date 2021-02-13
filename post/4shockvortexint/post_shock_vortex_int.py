#!/usr/bin/env python3
import glob, os
import numpy as np
import matplotlib.pyplot as plt
import read_fields as read
import scipy.ndimage.filters as filters

x, y = read.read_grid('grid.bin')
X, Y = np.meshgrid(x, y, indexing='ij')

dx = x[1]-x[0]

list_restart = []
for file in sorted(glob.glob('restart*.bin')):
   list_restart.append(file)

rho, u, v, E = read.read_restart(list_restart[0])

gamma=1.4
prs = (gamma-1)*rho*(E - 0.5*u**2)

# Compute deltap as in Bogey et al. JCP 2009
prs_pre_shock  = prs[ 1, 1]
prs_post_shock = prs[-1,-1]

def compute_deltap(prs):
   deltap = np.asarray(prs)
   for i in range(len(x)):
      if x[i]>0:
         pnorm = prs_post_shock
      else:
         pnorm = prs_pre_shock
      deltap[i,:] = prs[i,:]/pnorm - 1.0
   return deltap

deltap   = compute_deltap(prs)
levels_p = np.linspace(-0.075, 0, 76)

isolines_p = [-0.1, -0.01, -0.005, -0.002, 0.002, 0.005, 0.01]

# Compute deltap as in Bogey et al. JCP 2009
def compute_laprho(rho):
   return filters.laplace(rho)/dx**2
laprho   = compute_laprho(prs)
levels_r = np.linspace(-np.max(laprho), np.max(laprho), 100)

# ------------------------------------------------------------------------------
# Prepare figures
plt.rc('xtick',labelsize=22)
plt.rc('ytick',labelsize=22)
plt.rc('text',usetex=True)

# list_restart = [list_restart[-1]]
for file in list_restart:
   print(file)
   rho, u, v, E = read.read_restart(file)
   prs = (gamma-1)*rho*(E - 0.5*u**2)
   # First figure - delta pressure
   deltap = compute_deltap(prs)
   plt.figure(figsize=(8,7))
   plt.xlim(-5,15)
   plt.ylim(-10,10)
   plt.axis('scaled')
   plt.contourf(X, Y, deltap, levels=levels_p, cmap='jet', extend='both')
   plt.colorbar()
   plt.contour( X, Y, deltap, levels=isolines_p, colors='k', linewidths=1.0)
   plt.tight_layout()
   plt.savefig('deltap_shock_vortex_'+file[7:16]+'.png', dpi=300)
   plt.close()
   # Second figure - rho laplacian
   laprho = compute_laprho(rho)
   plt.figure(figsize=(8,7))
   plt.xlim(-5,15)
   plt.ylim(-10,10)
   plt.axis('scaled')
   plt.contourf(X, Y, laprho, levels=levels_r, cmap='bwr', extend='both')
   plt.colorbar()
   plt.tight_layout()
   plt.savefig('laprho_shock_vortex_'+file[7:16]+'.png', dpi=300)
   plt.close()