#!/usr/bin/env python

## Utility functions for processing GALEV model outputs
import numpy as np
import os
from pygoods import sextractor
from scipy.integrate import *
from scipy.interpolate import interp1d

pc_cm = 3.0857e18  # 1 parsec in cm
area_10pc = 4. * np.pi * (10 * pc_cm)**2

def writeTemplates(libraryFile, statFile, root=None, overwrite=True, ageUniverse=13.8, tempdir='/Users/khuang/Dropbox/codes/phot_z/eazy/GALEV'):
   """
   From a GALEV output template library *_spec.dat, write individual SED 
   template files. All files share a common root, which by default is what
   goes before the _spec.data in the library file name.
   Also read the physical property file and add stellar masses, gas masses, 
   SFR, and metallicity (Z) to each file.
   """
   libraryFile = os.path.split(libraryFile)[-1]
   statFile = os.path.split(statFile)[-1]
   if root == None:
      root = libraryFile.split('_spec.dat')[0]
   f = open(libraryFile, 'rb')
   lines = f.readlines()
   f.close()
   f2 = open(statFile, 'rb')
   lines2 = f2.readlines()
   f2.close()
   # First, read the number of age and wavelength steps
   l0 = lines[0].split()[1:]
   numAges = int(l0[0])   # number of age steps
   numWaves = int(l0[1])  # number of wavelength points
   print "numAges, numWaves:", numAges, numWaves
   # Second, read the ages
   l1 = lines[1].split()[1:]
   agetable = np.array([float(x) for x in l1])
   readWave = False
   fluxtable = np.zeros((numWaves, numAges))
   wavetable = np.zeros(numWaves)
   for i in range(numWaves):
      l = lines[i+3].split()
      # print "len(l):", len(l)
      wavetable[i] = float(l[1])
      fluxtable[i] = [float(y) for y in l[2:]]
   # Also read the physical properties
   GMass = np.zeros(numAges)  # Gas mass
   SMass = np.zeros(numAges)  # Stellar mass
   SFR = np.zeros(numAges)    # SFR
   logZ = np.zeros(numAges)   # log(Z)
   m = 0  # counter for GMass, SMass, etc.
   for l2 in lines2:
      if l2.startswith('#'):
         pass
      else:
         list2 = l2.split()  # split the columns in statFile
         GMass[m] = float(list2[1])  
         SMass[m] = float(list2[2])
         SFR[m] = float(list2[3])
         logZ[m] = float(list2[4])
         m += 1
   # Now write to individual files
   if not os.path.isdir(tempdir):
      os.mkdir(tempdir)
   previous = "%.4f" % 0.
   previousAge = 0.
   includeAge = np.zeros(numAges, 'int')
   for j in range(numAges):
      currentAge = agetable[j]
      if j % 500 == 0:
         print j
      if currentAge / 1.e9 > ageUniverse:
         print "Stop at age = %f Gyr (j = %d)" % (agetable[j-1]/1.e9, j-1)
         break
      logAge = np.log10(currentAge)
      if "%.4f" % logAge == previous:
         raise ValueError, "logAge (%.4f) is not unique!" % logAge
      previous = "%.4f" % logAge
      filename = '%s_age%.4f.dat' % (root, logAge)
      if os.path.exists(filename) and overwrite==False:
         continue
      # Only write an age step if its distance from the previous step is 
      # larger than 5 percent of the previous age
      elif (currentAge - previousAge) >= (0.05 * previousAge):
         includeAge[j] = 1
         fj = open('%s/%s_age%.4f.dat' % (tempdir, root, logAge), 'wb')
         fj.write('# Root = %s\n' % root)
         fj.write('# Age = %.5e\n' % agetable[j])
         fj.write('# GMass = %.5e\n' % GMass[j])
         fj.write('# SMass = %.5e\n' % SMass[j])
         fj.write('# SFR = %.5e\n' % SFR[j])
         fj.write('# logZ = %.5e\n' % logZ[j])
         fj.write('# angstrom   f_lam\n')
         for i in range(numWaves):
            fj.write('%6.2f %.6e\n' % (wavetable[i], fluxtable[i][j]))
         fj.close()
         previousAge = currentAge
   print np.sum(includeAge)
   # Also write EAZY template file
   ft = open('%s.spectra.param' % root, 'wb')
   k = 1
   for j in range(numAges):
      if includeAge[j] == 1:
         logAge = np.log10(agetable[j])
         filename = '%s_age%.4f.dat' % (root, logAge)
         ageGyr = agetable[j] / 1.e9
         ft.write('%5d  %s/%-30s  1.0 %.7f 1.0 \n' % (k, tempdir, filename, ageGyr))
         k = k + 1
   ft.close()
   print "Done."

def calc_phys_prop(specfile, logAge):
   # Calculate, for each model spectrum, the physical properties to be used
   # in LePhare
   s = sextractor(specfile)
   f = interp1d(s._1, s._2)  # interpolated spectrum
   age = 10. ** logAge  # in years
   # c = 3.e18  angstrom/sec
   LUV = romberg(f, 2100., 2500.) / 400. * (2300.**2 / 3.e18) * area_10pc
   LR = romberg(f, 5500., 6500.) / 1000. * (6000.**2 / 3.e18) * area_10pc
   LK = romberg(f, 21000., 23000.) / 2000. * (22000.**2 / 3.e18) * area_10pc
   LIR = -99.0
