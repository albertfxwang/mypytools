#!/usr/bin/env python

from numpy import *
from plotting import plot_mcmc as pmc
import pyfits
import fitsutil
from pygoods import *

# First calculate phistar for B-dropouts

cb1 = Ftable('bdrops_mcmc_090612.fits')
cv1 = Ftable('vdrops_mcmc_090612.fits')

traces_b1 = array([cb1.alpha,cb1.mstar,cb1.logr0,cb1.sigma,cb1.beta])
traces_v1 = array([cv1.alpha,cv1.mstar,cv1.logr0,cv1.sigma,cv1.beta])
phistar_b1 = pmc.calc_phistar_parallel(traces_b1, 'b', nproc=4)
# ...wait until it finishes
fitsutil.add_columns('bdrops_mcmc_090612.fits',['phistar'],[phistar_b1],['D'])
print "Finished phistar for B-dropouts."

# Now calculate phistar for V-dropouts
phistar_v1 = pmc.calc_phistar_parallel(traces_v1, 'v', nproc=4)
fitsutil.add_columns('vdrops_mcmc_090612.fits',['phistar'],[phistar_v1],['D'])
print "Finished phistar for V-dropouts."