#!/usr/bin/env python

import numpy as np
import os, sys, glob, string, shutil, time
from pyraf import iraf
from iraf import artdata, images, stsdas
curdir = os.getcwd()
import tphot_sim
import yaml

datum = time.strftime("%m%d",time.localtime())

def run_tphotsim(parfile):
   sim = tphot_sim.TPHOTSim(parfile)
   while sim.n <= sim.nstop:
      print "Iteration %d:" % sim.n
      tphotparfile_n = '%s_run%d.tphot.yml' % (os.path.splitext(parfile)[0], sim.n)
      sim.insert_fake_sources()
      status = sim.run_sextractor()
      if status:
         sim.write_tphot_param(tphotparfile_n)
         tpipe = sim.run_tphot(tphotparfile_n)
         sim.collect_results(tpipe)
         os.chdir(curdir)
         if not sim.save:
            tpipe.clear_all()
            sim.cleanup()
         # Also clean up symbolic links
         tpars = yaml.load(open(tphotparfile_n))
         os.chdir(tpars['hiresdir'])
         os.system('rm %s*.fits' % sim.root)
         os.system('rm *_segflg_%s_*.fits' % sim.root)
         os.system('rm -r %s_%s' % (sim.root, sim.hires_band))
         os.chdir(curdir)
         os.system('rm %s_run%d.cat' % (sim.root, sim.n))
         os.system('rm %s_run%d.tphot.yml' % (sim.root, sim.n))
      sim.n = sim.n + 1

if __name__ == "__main__":
   parfile = sys.argv[1]
   run_tphotsim(parfile)

   print "Finished TPHOT simulation!"