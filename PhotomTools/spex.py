#!/usr/bin/env python

# Read the SpeX dwarf star template spectra

import pysynphot as S
import numpy as np
import matplotlib.pyplot as plt


class spexTemplate(S.spectrum.FileSourceSpectrum):
   def __init__(self, filename):
      f = open(filename)
      lines = f.readlines()
      for l in lines:
         if 'spectral type' in l:
            try:
               l2 = l.split(':')[-1]
               self.specType = l2.split()[0]
            except:
               print "File name: ", filename
               raise IndexError
         elif 'J magnitude' in l:
            l2 = l.split('=')[-1].split()[0]
            self.jmag = float(l2)
         elif 'H magnitude' in l:
            l2 = l.split('=')[-1].split()[0]
            self.hmag = float(l2)
         elif 'Ks magnitude' in l:
            l2 = l.split('=')[-1].split()[0]
            self.ksmag = float(l2)
      f.close()
      print "Spectral type: %s" % self.specType
      S.spectrum.FileSourceSpectrum.__init__(self, filename)

def plotStar(filename, ax=None, normwave=1.6, normFnu=0.43):
   print "File name: ", filename
   if ax == None:
      fig = plt.figure()
      ax = fig.subplot(111)
   star = spexTemplate(filename)
   print "Spectral type: %s" % star.specType
   fnu = star(star.wave) * (star.wave**2 / 3e8)  # convert from flam to fnu
   fnu0 = star(normwave) * (normwave**2 / 3e8)
   fnu = (fnu / fnu0) * normFnu
   ax.plot(star.wave * 1.e4, fnu, label='class=%s' % star.specType)
   plt.draw()
   return ax
