#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import fit_lbg as fl


alpha_arr = arange(-2.4, -1.38, 0.02)
Mstar_arr = arange(-22.0, -19.0, 0.02)
logr0_arr = arange(0.1, 1.02, 0.02)
sigma_arr = arange(0.1, 1.02, 0.02)
beta_arr  = arange(0.02, 1.02, 0.02)
p2 = array([-2.078, -20.840, 0.799, 0.671, 0.278])
par_arr = [alpha_arr, Mstar_arr, logr0_arr, sigma_arr, beta_arr]

def plot_logl12():
   # plot the shapes of log(L1) and log(L2), where L1 and L2 are the likelihood function
   # of the GOODS & HUDF models, respectively, and see the differences in their shapes
   plt.figure()
   for i in range(5):
      plt.subplot(2,3,i+1)
      print i
      logl1 = zeros(len(par_arr[i]))
      logl2 = zeros(len(par_arr[i]))
      for j in range(len(logl1)):
         par = p2.copy()
         par[i] = par_arr[i][j]
         l1 = logl1[j] = fl.quickmlfunc(par, 1.0, -1.0)
         l2 = logl2[j] = fl.quickmlfunc(par, -1.0, 1.0)
         plt.plot(par_arr[i], logl1, c='blue', label='L1')
         plt.plot(par_arr[i], logl2 * (logl1[0]/logl2[0]), c='red', label='L2')
         plt.plot(par_arr[i], logl1+logl2, c='black')

   return 0


