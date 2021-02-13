#!/usr/bin/env python3
import glob, os
import numpy as np
import matplotlib.pyplot as plt
import read_fields as read
import numpy.fft as mklfft
import tools_fft as tl
import mesfonction as mes_fonc



x, y = read.read_grid('grid.bin')
X, Y = np.meshgrid(x, y, indexing='ij')
gamma = 1.4
Pinf = 1e5
Uscale = 112.249721603218
nu = 5.612486080160912E-002
k  = 4
### creer dossier output 
nrk =int(input("saisir nrk:  "))
deriv_conv_order =int(input("deriv_conv_order:  "))
deriv_visc_order =int(input("deriv_visc_order:  "))
chemin = mes_fonc.cree_repertoir(nrk,deriv_conv_order,deriv_visc_order)

def f(t):
   return np.exp(-2*nu*k**2*t/Uscale)

def exact_sol(t):
   u_ex =  Uscale*np.sin(k*X)*np.cos(k*Y)*f(t)
   v_ex = -Uscale*np.cos(k*X)*np.sin(k*Y)*f(t)
   p_ex = Pinf + 0.25*Uscale**2*f(t)**2*(np.cos(2*k*X) + np.cos(2*k*Y))
   return u_ex, v_ex, p_ex

list_restart = []
for file in sorted(glob.glob('restart*.bin')):
   list_restart.append(file)

levels_u = np.linspace(-1, 1, 101)
levels_p = np.linspace(0.95, 1.05, 101)

# ------------------------------------------------------------------------------
# Prepare figures
plt.rc('xtick',labelsize=22)
plt.rc('ytick',labelsize=22)
plt.rc('text',usetex=True)

# list_restart= [list_restart[-1]]
for file in list_restart:
   print(file)
   rho, u, v, E = read.read_restart(file)
   # p = (gamma-1)*rho*(E - 0.5*(u**2+v**2))
   plt.figure(figsize=(8,6))
   plt.axis('equal')
   plt.contourf(X, Y, u/Uscale, levels=levels_u, cmap='jet', extend='both')
   plt.colorbar()
   plt.savefig(chemin+'tgv_u_'+file[7:16]+'.png', dpi=300)
   plt.tight_layout()
   # plt.show()
   plt.close()

   time = float(file[7:11])
   u_ex, v_ex, p_ex = exact_sol(time)

   plt.figure(figsize=(8,6))
   plt.axis('equal')
   plt.contourf(X, Y, u_ex/Uscale, levels=levels_u, cmap='jet', extend='both')
   plt.colorbar()
   plt.savefig(chemin+'tgv_uex_'+file[7:16]+'.png', dpi=300)
   plt.tight_layout()
   # plt.show()
   plt.close()


   plt.figure(figsize=(8,6))
   plt.axis('equal')
   plt.contourf(X, Y, (u-u_ex)/Uscale*100, levels=100, cmap='jet', extend='both')
   plt.colorbar()
   plt.savefig(chemin+'tgv_uerr_'+file[7:16]+'.png', dpi=300)
   plt.tight_layout()
   # plt.show()
   plt.close()

for file in list_restart:
   print(file)
   rho, u, v, E = read.read_restart(file)

   N = u.shape[0]
   k,spc = tl.power_spectra(u/Uscale,v/Uscale)
   k = k-1

   plt.figure(figsize=(8,6))
   plt.xlabel(r'$k$',fontsize=24)
   plt.ylabel(r'$E(k)$',fontsize=24)
   plt.xlim(1,120)
   plt.loglog([6,80],[1, np.exp(np.log10(1)-3.*np.log10(80)-np.log10(6))] , color='k')
   plt.ylim(5e-10,10)
   plt.text(20,0.1, r'$-3$', fontsize=22)

   plt.loglog(k[1:], spc[1:], color='r')
   plt.tight_layout()
   plt.savefig(chemin+'tgv_fft_'+file[7:16]+'.png', dpi=300)
   plt.close()