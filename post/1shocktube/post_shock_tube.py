#!/usr/bin/env python3
import glob, os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mesfonction as mes_fonc

x, ro, rou, rov, roE = np.loadtxt( fname='restart0000_2000.dat', unpack=True )
### creer dossier output 
nrk =int(input("saisir nrk:  "))
deriv_conv_order =int(input("deriv_conv_order:  "))
deriv_visc_order =int(input("deriv_visc_order:  "))
chemin = mes_fonc.cree_repertoir(nrk,deriv_conv_order,deriv_visc_order)

gamma = 1.4
rinf  = 1.0
pinf  = 1e5
uinf  = (np.sqrt(pinf/rinf))


u = rou/ro
p = (gamma-1)*(roE - 0.5*ro*u**2)

# ------------------------------------------------------------------------------
# Prepare figures
plt.rc('xtick',labelsize=22)
plt.rc('ytick',labelsize=22)
plt.rc('text',usetex=True)

x_ex, r_ex, p_ex, u_ex, c_ex = np.loadtxt(fname='REF/sod_exact.dat', skiprows=10, unpack=True)

fig1 = plt.figure(1,figsize=(8,6))
ax = fig1.add_subplot(111)
ax.plot(x_ex, r_ex, color='grey', ls='-' )
ax.plot(x_ex, p_ex, color='grey', ls='-' )
ax.plot(x_ex, u_ex, color='grey', ls='-' )

ax.plot(x, ro/rinf, color='r'     , ls='--', label=r'$\rho$')
ax.plot(x,  p/pinf, color='orange', ls='--', label=r'$p$')
ax.plot(x,  u/uinf, color='b'     , ls='--', label=r'$u$')

ax.set_xlabel(r'$x$',fontsize=24)
plt.legend(loc=0,fontsize='x-large')
plt.tight_layout()
plt.savefig(chemin+"shock_tube.png")
#plt.show()