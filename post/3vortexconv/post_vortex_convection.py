#!/usr/bin/env python3
import glob, os
import numpy as np
import matplotlib.pyplot as plt
import read_fields as read
import mesfonction as mes_fonc

x, y = read.read_grid('grid.bin')
X, Y = np.meshgrid(x, y, indexing='ij')

list_restart = []
for file in sorted(glob.glob('restart*.bin')):
   list_restart.append(file)

levels = np.linspace(0.55, 1.22, 68)
### creer dossier output 
nrk =int(input("saisir nrk:  "))
deriv_conv_order =int(input("deriv_conv_order:  "))
deriv_visc_order =int(input("deriv_visc_order:  "))
chemin = mes_fonc.cree_repertoir(nrk,deriv_conv_order,deriv_visc_order)
# ------------------------------------------------------------------------------
# Prepare figures
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rc('text',usetex=True)

# list_restart = [list_restart[-1]]
for file in list_restart:
   print(file)
   rho, u, v, E = read.read_restart(file)
   plt.figure(figsize=(8,6))
   plt.axis('scaled')
   plt.xlim(0,20)
   plt.ylim(0,20)
   plt.contourf(X, Y, rho, levels=levels, cmap='jet', extend='both')
   plt.colorbar()
   plt.tight_layout()
   plt.savefig(chemin+'vortex_conv_'+file[7:16]+'.png')
   plt.close()

rho_beg, u_beg, v_beg, E_beg = read.read_restart(list_restart[ 0])
rho_end, u_end, v_end, E_end = read.read_restart(list_restart[-1])

rho_err = np.abs(rho_end-rho_beg)/rho_beg

print("L2 Error norm:", np.linalg.norm(rho_err,2))

levels = np.linspace(-4, -1, 51)
plt.figure(figsize=(8,6))
plt.axis('scaled')
plt.xlim(0,20)
plt.ylim(0,20)
plt.contourf(X, Y, np.log10(rho_err), levels=levels, cmap='jet', extend='both')
plt.colorbar()
plt.tight_layout()
plt.savefig(chemin+'vortex_rho_error.png')
plt.close()