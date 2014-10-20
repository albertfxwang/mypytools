#!/usr/bin/env python

# run_simulations file1.sim file2.sim ... 
# This script adds fake galaxies to real images, runs sextractor
# with an association list, and outputs the SExtractor parameters
# for the matched sources.
# 
# Configuration parameters (file.sim...this list will probably grow):
#
# SEXFILE		file.sex
# PSFFILE		""
# NITER			1
# REALIMAGE		fakez_drz.fits	
# FLAGIMAGE		fakez_flag.fits	
# SAVE			no
# NGALAXIES		100
# CIRCULAR              no
# DISKFRAC              0.5
# SCALE			0.01
# MAGLOW		20 
# MAGHIGH		28 
# MAGZPT		24.961 
# RMIN			0.01	# minimum input radius arcsec
# MAX			1.5	# maximum input radius arcsec
# RUN_SEXTRACTOR        yes     
# MEASURE_PETROSIAN     yes     
# LOGNORMAL_MAG0	24.0	# Not yet implemented
# LOGNORMAL_PEAK	0.5
# LOGNORMAL_WIDTH	0.5
# 
#

# This works for GALFIT simulations. Without the GALFIT part, it can also be
# used for SExtractor simulations, but in that case one needs to specify 
# the SED of each fake galaxy.

# mkobjects creates noiseless devauc-profile galaxies that are 
# 0.13 mag too faint when dynrange=1.e5
devauc_correction = -0.00  


# from numpy import *
import numpy as np
import os, sys, glob, string, shutil, time
from pygoods import *
from pyraf import iraf
from iraf import artdata, images, stsdas
curdir = os.getcwd()
import galfit_sim

datum = time.strftime("%m%d",time.localtime())
gsigma = 10.   # sigma of Gaussian wing in pixels

def run_galfit_sim(parfile):
   gfsim = galfit_sim.galfit_sim(parfile)
   while gfsim.n < gfsim.nstop:
      print "Iteration %d:" % gfsim.n
      gfsim.insert_fake_sources()
      gfsim.run_sextractor()
      if gfsim.pixscale_galfit != gfsim.pixscale:
         print "Also generate fake galaxies for GALFIT measurement images..."
         gfsim.resample_segmap()
         gfsim.insert_fake_sources_galfit()
      gfsim.run_galfit()
      gfsim.match_galfit_catalog()
      gfsim.cleanup()
      gfsim.n = gfsim.n + 1

if __name__ == '__main__':
   # Read the configuration file and construct arguments for 
   #	simulate() and measure_petrosian()
 
   t1 = time.time()
   curdir = os.getcwd()
   
   parfile = sys.argv[1]
   run_galfit_sim(parfile)
   
   print "Finished simulation"
