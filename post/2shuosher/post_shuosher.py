#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import sys
import glob

gamma = 1.4
rinf  = 1.0
pinf  = 1e5
uinf  = (np.sqrt(pinf/rinf))
nrk =int(input("saisir nrk: "))
if nrk == 6 :
   chemin = 'images/rung_kutta_6/'
elif nrk == 4 :
   chemin = 'images/rung_kutta_4/'
elif nrk == 2 :
   chemin = 'images/rung_kutta_2/'
else :
   chemin = ''

# ------------------------------------------------------------------------------

x_e, r_e, p_e, u_e, e_e, c_e = np.loadtxt(fname='./REF/shu_osher.dat', skiprows=2, unpack=True)

# Prepare figures
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.rc('text',usetex=True)

list_restart = []
for file in sorted(glob.glob('restart*.dat')):
   list_restart.append(file)

# list_restart = [list_restart[-1]]  # To use only last field

for file in list_restart:
   print(file)
   x, ro, rou, rov, roE = np.loadtxt( fname=file, unpack=True )
   u = rou/ro
   p = (gamma-1)*(roE - 0.5*ro*u**2)
   u = u/uinf
   p = p/pinf

   fig1 = plt.figure(1,figsize=(22,12))
   ax1 = fig1.add_subplot(241)
   ax1.plot(x_e, r_e, color='k', ls='-' )
   ax1.plot(x  , ro, color='r', ls='-' , label=r'$\rho$')
   ax1.set_xlim(-5,5)
   ax1.set_ylim(0.5,5)
   ax1.set_title(r'$\rho$',fontsize=22)

   ax2 = fig1.add_subplot(242)
   ax2.plot(x_e, p_e, color='k', ls='-' )
   ax2.plot(x  , p, color='orange', ls='-' , label=r'$p$')
   ax2.set_xlim(-5,5)
   ax2.set_ylim(0,12)
   ax2.set_title(r'$p$',fontsize=22)

   ax3 = fig1.add_subplot(243)
   ax3.plot(x_e, u_e, color='k', ls='-' )
   ax3.plot(x  , u, color='b', ls='-' , label=r'$u$')
   ax3.set_xlim(-5,5)
   ax3.set_ylim(-0.1, 3.)
   ax3.set_title(r'$u$',fontsize=22)

   ax4 = fig1.add_subplot(244)
   ax4.plot(x_e, np.log(p_e/r_e**1.4), color='k', ls='-' )
   ax4.plot(x  , np.log(p/ro**1.4), color='g', ls='-' , label=r'$s$')
   ax4.set_xlim(-5,5)
   ax4.set_ylim(-0.4, 0.8)
   ax4.set_title(r'$s$',fontsize=22)

   ax5 = fig1.add_subplot(245)
   ax5.plot(x_e, r_e, color='k', ls='-' )
   ax5.plot(x  , ro, color='r', ls='-', marker='o', label=r'$\rho$')
   ax5.set_xlim(0.0,2.5)
   ax5.set_ylim(2.9,5.1)
   ax5.set_xlabel(r'$x$',fontsize=22)

   ax6 = fig1.add_subplot(246)
   ax6.plot(x_e, p_e, color='k', ls='-' )
   ax6.plot(x  , p, color='orange', ls='-', marker='o', label=r'$p$')
   ax6.set_xlim(-3.5,3.0)
   ax6.set_ylim(8,12)
   ax6.set_xlabel(r'$x$',fontsize=22)

   ax7 = fig1.add_subplot(247)
   ax7.plot(x_e, u_e, color='k', ls='-' )
   ax7.plot(x  , u, color='b', ls='-', marker='o', label=r'$u$')
   ax7.set_xlim(-3.5,3.0)
   ax7.set_ylim(2.45,2.85)
   ax7.set_xlabel(r'$x$',fontsize=22)

   ax8 = fig1.add_subplot(248)
   ax8.plot(x_e, np.log(p_e/r_e**1.4), color='k', ls='-' )
   ax8.plot(x  , np.log(p/ro**1.4), color='g', ls='-', marker='o', label=r'$s$')
   ax8.set_xlim( 0.5,2.5)
   ax8.set_ylim(0.20,0.75)
   ax8.set_xlabel(r'$x$',fontsize=22)

   plt.tight_layout()
   plt.savefig(chemin+'shu_osher_'+file[7:16]+'.png', dpi=300)
   #plt.show()
   plt.close()