#!/usr/bin/env python

## Utility functions for processing GALEV model outputs
import numpy as np
import os
from pygoods import sextractor
from scipy.integrate import *
from scipy.interpolate import interp1d

pc_cm = 3.0857e18  # 1 parsec in cm
area_10pc = 4. * np.pi * (10 * pc_cm)**2

def writeTemplates(specFile, statFile, root=None, overwrite=True, ageUniverse=13.8, tempdir='/Users/khuang/Dropbox/codes/phot_z/eazy/GALEV'):
   """
   From a GALEV output template library *_spec.dat, write individual SED 
   template files. All files share a common root, which by default is what
   goes before the _spec.data in the library file name.
   Also read the physical property file and add stellar masses, gas masses, 
   SFR, and metallicity (Z) to each file.
   """
   specFile = os.path.split(specFile)[-1]
   statFile = os.path.split(statFile)[-1]
   if root == None:
      root = specFile.split('_spec.dat')[0]
   f = open(specFile, 'rb')
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
   ## Since LePhare hates header lines, I need an extra file to store the 
   ## physical properties of each model SED to be used later to make the *.phys
   ## file.
   fstat = open(tempdir + '/' + root+'_stat.txt', 'wb')
   fstat.write('# 1 FILENAME\n')
   fstat.write('# 2 AGE\n')
   fstat.write('# 3 GMASS\n')
   fstat.write('# 4 SMass\n')
   fstat.write('# 5 SFR\n')
   fstat.write('# 6 logZ\n')
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
         ### Le Phare DOESN'T LIKE ANY HEADER LINES IN THE SPEC FILES... 
         # fj.write('# Root = %s\n' % root)
         # fj.write('# angstrom   f_lam\n')
         for i in range(numWaves):
            fj.write('%6.2f %.12f\n' % (wavetable[i], fluxtable[i][j]))
         fj.close()
         previousAge = currentAge
         ## write physical properties
         fstat.write('%s  ' % filename)
         fstat.write('%.5e  ' % agetable[j])
         fstat.write('%.5e  ' % GMass[j])
         fstat.write('%.5e  ' % SMass[j])
         fstat.write('%.5e  ' % SFR[j])
         fstat.write('%.5e  ' % logZ[j])
         fstat.write('\n')
   fstat.close()
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

def calc_phys_prop(specfile):
   # Calculate, for each model spectrum, the physical properties to be used
   # in LePhare
   # Manual input tau (should be improved later!!)
   s = sextractor(specfile)
   f = interp1d(s._1, s._2)  # interpolated spectrum
   # c = 3.e18  angstrom/sec
   LUV = romberg(f, 2100., 2500.) / 400. * (2300.**2 / 3.e18) * area_10pc
   LR = romberg(f, 5500., 6500.) / 1000. * (6000.**2 / 3.e18) * area_10pc
   LK = romberg(f, 21000., 23000.) / 2000. * (22000.**2 / 3.e18) * area_10pc
   # LIR = -99.0
   D4000 = romberg(f, 4050, 4250) / romberg(f, 3750, 3950)
   return [LUV, LR, LK, D4000]
   # return [age, LUV, LR, LK, LIR, SMass, SFR, Z, tau, D4000]

def write_phys_prop(specfiles, statFile, physfile, tau=1.e9):
   # Manual input tau (should be improved later!!)
   s2 = sextractor(statFile)
   ff = open(physfile, 'wb')  # the *.phys file
   i = 1  # line running count; 1st column of *.phys file
   j = 1  # model running count; 2nd column of *.phys file
   ## Le Phare seems to require a dummy line for each model with negative 
   ## values? Don't understand why...
   for f in specfiles:
      LUV, LR, LK, D4000 = calc_phys_prop(f)
      k = np.arange(len(s2))[s2.filename==f][0]  # matching index in statFile
      age = s2.age[k]
      SMass = s2.smass[k]
      SFR = s2.sfr[k]
      Z = 10.**(s2.logz[k])
      # first, write a dummy line
      # ff.write("%d  %d " % (i, j))
      # ff.write("%.6E  %.6E  %.6E  " % (-1, -99, -99))
      # ff.write("%.6E  %.6E  %.6E  " % (-99, -99, -99))
      # ff.write("%.6E  %.6E  %.6E  " % (-99, -99, tau))
      # ff.write("%.6E  " % -99.0)
      # ff.write('\n')
      # i += 1
      # Now write the real stuff
      ff.write("%d  %d  " % (i, j))
      ff.write("%.6E  %.6E  %.6E  " % (age, np.log10(LUV), np.log10(LR)))
      ff.write("%.6E  %.6E  %.6E  " % (np.log10(LK), -99.0, SMass))
      ff.write("%.6E  %.6E  %.6E  " % (SFR, Z, tau))
      ff.write("%.6E  " % D4000)
      ff.write("\n")
      i += 1
      j += 1
   ff.close()
