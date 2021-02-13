#!/usr/bin/env python3
import glob, os
import numpy as np

def read_grid(filename):
   f=open(filename,'r')
   nx, ny = np.fromfile(f, dtype=np.int32, count=2)
   x = np.fromfile(f, dtype=np.float64, count=nx)
   y = np.fromfile(f, dtype=np.float64, count=ny)
   return x, y

def read_restart(filename):
   f=open(filename,'r')
   nx, ny = np.fromfile(f, dtype=np.int32, count=2)
   rho  = np.fromfile(f, dtype=np.float64, count=nx*ny).reshape((nx,ny), order='F')
   rhou = np.fromfile(f, dtype=np.float64, count=nx*ny).reshape((nx,ny), order='F')
   rhov = np.fromfile(f, dtype=np.float64, count=nx*ny).reshape((nx,ny), order='F')
   rhoE = np.fromfile(f, dtype=np.float64, count=nx*ny).reshape((nx,ny), order='F')
   return rho, rhou/rho, rhov/rho, rhoE/rho