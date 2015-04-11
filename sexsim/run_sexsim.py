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
from pyraf import iraf
from iraf import artdata, images, stsdas
curdir = os.getcwd()
import sextractor_sim

datum = time.strftime("%m%d",time.localtime())

def run_sexsim(parfile):
   sim = sextractor_sim.sextractor_sim(parfile)
   while sim.n < sim.nstop:
      print "Iteration %d:" % sim.n
      sim.insert_fake_sources()
      sim.run_sextractor()
      if not sim.save:
         sim.cleanup()
      sim.n = sim.n + 1

if __name__ == '__main__':
   # Read the configuration file and construct arguments for 
   #	simulate() and measure_petrosian()
 
   t1 = time.time()
   curdir = os.getcwd()
   
   parfile = sys.argv[1]
   run_sexsim(parfile)
   
   print "Finished simulation"
