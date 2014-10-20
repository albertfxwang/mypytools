#!/usr/bin/env python

import numpy as np
import os
from pygoods import sextractor

def read_chains(chainfile, mcmc_cat, nburn=0):
   c = sextractor(chainfile)
   nwalkers = len(np.unique(c._1))
   nstart = nburn * nwalkers  # ignore the first nburn chains
   alpha = c._2[nstart:]
   mstar = c._3[nstart:]
   logr0 = c._4[nstart:]
   sigma = c._5[nstart:]
   beta = c._6[nstart:]
   phistar = c._7[nstart:]
   f = open(mcmc_cat, 'wb')
   f.write('# 1 ALPHA\n')
   f.write('# 2 MSTAR\n')
   f.write('# 3 LOGR0\n')
   f.write('# 4 SIGMA\n')
   f.write('# 5 BETA\n')
   f.write('# 6 PHISTAR\n')
   for i in range(len(alpha)):
      f.write('%f %f ' % (alpha[i], mstar[i]))
      f.write('%f %f ' % (logr0[i], sigma[i]))
      f.write('%f %.8f ' % (beta[i], phistar[i]))
      f.write('\n')
   f.close()
