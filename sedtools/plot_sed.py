#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pysynphot as S
from bands import filters
import os, sys
from my_mplfonts import Helvetica
from magredshift import mag_redshift

# Custom marker styles for matplotlib
downarrow = [(0,0),(0,-4),(-1,-3),(0,-4),(1,-3),(0.,-4)]
cosmopars = {'H0':70.0, 'omega_m':0.3, 'omega_l':0.73}

class plot_sed(object):
   def __init__(self, bands, abmags, magerrs, z=0., name="", waveunit='angstrom',
                cosmopars=cosmopars):
      """
      Define the photomry of this object.
      At initialization, give the AB magnitudes and filters of this object, as 
      well as its name and redshift (default to 0.)
      Input magnitudes should be 99.0 if S/N < 1. In that case, the corresponding
      magnitude errors reflect the 1-sigma magnitude limit.
      """
      fig = plt.figure()
      self.ax = fig.add_subplot(111)
      self.ax.invert_yaxis()
      self.name = name
      self.H0 = cosmopars['H0']
      self.omega_m = cosmopars['omega_m']
      self.omega_l = cosmopars['omega_l']
      abmags = np.array(abmags)
      magerrs = np.array(magerrs)
      self.abmags_plot = np.where(abmags<99.0, abmags, magerrs)
      self.magerrs_plot = np.where(abmags<99.0, magerrs, 0.)
      self.uplim = np.where(abmags<99.0, False, True)
      self.bands = bands
      self.pivot = np.zeros(len(bands))
      self.waveunit = waveunit
      self.bandindex = {}
      for i in range(len(self.bands)):
         self.bandindex[self.bands[i]] = i
      # Pivot wavelengths in angstroms
      for i in range(len(bands)):
         if bands[i] not in filters.keys():
            print "Filter %s not in database..." % bands[i]
            sys.exit()
         plam = filters[bands[i]].pivot()
         if waveunit == 'angstrom':
            pass
         elif waveunit == 'micron':
            plam = plam / 1.e4
         self.pivot[i] = plam
      self.z = z

   def set_ticklabels(self):
      self.ax.set_xticklabels(self.ax.get_xticks(), 
                              font_properties=Helvetica(14))
      self.ax.set_yticklabels(self.ax.get_yticks(), 
                              font_properties=Helvetica(14))

   def plot_phot(self, ms=10**2, mfc='red', mec='none'):
      self.ax.scatter(self.pivot[self.abmags_plot>0.], 
                      self.abmags_plot[self.abmags_plot>0.], 
                      marker='o', facecolor=mfc,
                      edgecolor=mec, s=ms)
      # Plot error bars for detections
      self.ax.errorbar(self.pivot[self.uplim==False], 
                      self.abmags_plot[self.uplim==False],
                      yerr=self.magerrs_plot[self.uplim==False], fmt=None, 
                      ecolor='black', elinewidth=1.5, capsize=np.sqrt(ms)-2)
      # Plot upper limit symbols
      self.ax.scatter(self.pivot[self.uplim==True],
                   self.abmags_plot[self.uplim==True],
                   marker=downarrow, s=(np.sqrt(ms)+2)**2, color='black')
      self.ax.set_ylabel('Magnitude (AB)', font_properties=Helvetica(18))
      self.ax.set_xlabel('Wavelength (%s)' % self.waveunit, 
                         font_properties=Helvetica(18))
      self.ax.set_title(self.name.upper(), font_properties=Helvetica(24))
      plt.draw()

   def plot_modelSED(self, sp, normfilter="", label="", lw=2.0):
      """
      Plot a model SED along with the photometric points.
      Convert the model SED into units of Fnu if it is not already it.
      normfilter specifies the filter at which the model SED normalizes to 
      match the photometry of the object. If empty, it defaults to the 
      reddest filter.
      """
      # wave = sp.wave.copy()
      if self.z > 0:
         sp = mag_redshift(sp, self.z, filters[self.bands[-1]], H0=self.H0, 
                           omega_m=self.omega_m, omega_l=self.omega_l)[1]
      if sp.fluxunits.name == 'fnu':
         flux = sp.flux.copy()
      elif sp.fluxunits.name == 'flam':
         flux = sp.flux * sp.wave**2  # miss a factor of c here, but it's OK
      if normfilter in self.bands:
         normlam = self.pivot[self.bandindex[normfilter]]
      else:
         normfilter = self.bands[-1]
         normlam = self.pivot[-1]
         # assume that self.bands is in increasing order in wavelength
      wave = sp.wave.copy()
      spmag = -2.5 * np.log10(flux)
      normflux = sp.sample(normlam) * normlam**2
      normmag = -2.5 * np.log10(normflux)
      normfactor = self.abmags_plot[self.bandindex[normfilter]] - normmag
      # print len(wave), len(spmag)
      self.ax.plot(wave, spmag + normfactor, linestyle='-', color='0.5', 
                   lw=lw, label=label)
      xmin = self.pivot[0] - filters[self.bands[0]].rectwidth() / 2.
      xmax = self.pivot[-1] + filters[self.bands[-1]].rectwidth() / 2.
      self.ax.set_xlim(xmin, xmax)
      ymin = self.abmags_plot[self.abmags_plot>0.].max() + 2.
      ymax = self.abmags_plot[self.abmags_plot>0.].min() - 2.
      self.ax.set_ylim(ymin, ymax)
      self.set_ticklabels()
      plt.draw()
      # return wave, spmag

