#!/usr/bin/env python

import numpy as N
from pygoods import sextractor, Ftable
import cPickle
from bivariate import zdist
import cosmoclass
import matplotlib.pyplot as plt
from dropout_selection import lbg_colorcrit as lcc

filter_dic = {'acs_f435w':'b','acs_f606w':'v','acs_f775w':'i',
   'acs_f850lp':'z','wfc3_f125w':'j','wfc3_f160w':'h'}
zeropoints = {'wfc3_f160w':25.960,'acs_f435w':25.673,'acs_f606w':26.505,
   'acs_f775w':25.678,'acs_f850lp':24.867,
   'wfc3_f125w':26.250}
def_bands = ['acs_f435w','acs_f606w','acs_f775w','acs_f850lp', 
             'wfc3_f125w', 'wfc3_f160w']

class zdist_idrops_mag(zdist.zdgrid):
   """
   A kludge class for calculating 1-D P(z) in apparent magnitude bins.
   """
   def __init__(self, M0, M1, dM, z0, z1, dz, area):
      self.M0 = M0
      self.M1 = M1
      self.dM = dM
      self.M1500arr = N.arange(M0, M1, dM)
      # turn magnitude into strings
      #self.magarr_str = map(lambda x:'%.1f'%x, self.magarr)
      self.z0 = z0
      self.z1 = z1
      self.dz = dz
      self.zarr = N.arange(z0, z1, dz)
      # This is for the bin edges of N.histogram
      self.zarr_edges = N.concatenate([self.zarr,[z1+dz]])
      self.area = area
      self.calc_dVdz()
      self.zdist = {}

   def __call__(self, simcatalog, bands=def_bands, key='Hua13',
                sn_lolim={'wfc3_f160w':5.0,'wfc3_f125w':3.5,'acs_f850lp':2.0},
                sn_hilim={'acs_f435w':2.0}, expand=[0.,0.]):
      ## Most stuff is the same as the 2D case:
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc('f775w_drop', key)
      c = Ftable(simcatalog)
      # Now initialize the attributes
      for b in bands:
         bs = filter_dic[b]  # short name of filter b
         flux_iso = getattr(c,'%s_flux_iso'%bs)
         fluxerr_iso = getattr(c,'%s_fluxerr_iso'%bs)
         mag_iso = getattr(c,'%s_mag_iso'%bs)
         # calculate 1-sigma upper limits on sources with S/N < 1
         mag_iso = N.where(flux_iso/fluxerr_iso>1.0,mag_iso,
                           zeropoints[b]-2.5*N.log10(fluxerr_iso))
         setattr(self,b+'_mag',mag_iso)  
         # MAG_ISO --> for color calculation
         setattr(self,b+'_flux',getattr(c,'%s_flux_auto'%bs))  
         # FLUX_AUTO
         setattr(self,b+'_flux_aper',getattr(c,'%s_flux_aper_2'%bs))
         # Aperture flux within 0.18" aperture (2.94 pixels)
         setattr(self,b+'_fluxerr',getattr(c,'%s_fluxerr_auto'%bs))  
         # FLUXERR_AUTO
         setattr(self,b+'_fluxerr_aper',getattr(c,'%s_fluxerr_aper_2'%bs))
         # Aperture flux errors within 0.18" aperture
         setattr(self,b+'_sn',
            getattr(self,'%s_flux'%b)/getattr(self,'%s_fluxerr'%b))  
         # S/N calculated using AUTO APERTURE 
         setattr(self,b+'_sn_aper',
            getattr(self,'%s_flux_aper'%b)/getattr(self,'%s_fluxerr_aper'%b))
         # S/N calculated using 0.18" APERTURE
      # Now construct S/N criteria
      #self.sncrit = N.ones(len(c.d),'bool') 
      self.sn_lolim_crit = N.ones(len(c.d), 'bool')
      self.sn_hilim_crit = N.ones(len(c.d), 'bool') 
      for b in sn_lolim.keys():  # enforce S/N lower limits
         self.sn_lolim_crit = self.sn_lolim_crit & \
               (getattr(self,b+'_sn')>=sn_lolim[b])
      for b in sn_hilim.keys():  
         # enforce S/N upper limits (veto bands)
         # but use fixed aperture S/N
         self.sn_hilim_crit = self.sn_hilim_crit & \
               (getattr(self,b+'_sn_aper')<sn_hilim[b])
      print "Total number of objects satisfying the S/N criteria:", \
         sum((self.sn_hilim_crit)&(self.sn_lolim_crit))
      self.sncrit = (self.sn_hilim_crit==True) & (self.sn_lolim_crit==True)
      print "Do selections..."
      self.color1 = self.acs_f775w_mag-self.acs_f850lp_mag
      self.color2 = self.acs_f850lp_mag-self.wfc3_f125w_mag
      self.lcc.select(self.color1,self.color2)  # do color selection!
      self.colorcrit = self.lcc.crit.copy()
      self.lcc.crit = self.lcc.crit & self.sncrit  # enforce S/N criteria
      self.dropcrit = self.lcc.crit & (c.detect==True)  # just in case
      self.detect = c.detect.copy()
      print "Selection done."
      print "Total number of objects in the catalog: %d" % len(c.d)
      print "Total number selected as i-dropouts: %d" % (sum(self.dropcrit))
      for M in self.M1500arr:
         # For each magnitude bin, calculate completeness P(z) 
         bincrit = (c.m1500_in >= M) & (c.m1500_in < (M+self.dM))
         z_in_bin = c.z_in[bincrit==True]
         z_in_bin_detect = c.z_in[(bincrit==True)&(c.detect==True)]
         zhist_bin_input = N.histogram(z_in_bin, self.zarr_edges)[0]
         zhist_bin_detect = N.histogram(z_in_bin_detect, self.zarr_edges)[0]
         self.zdist['%.1f'%M] = zhist_bin_detect.astype('float') / \
               N.maximum(zhist_bin_input.astype('float'), 1.0)
         self.zdist['%.1f'%M] = self.zdist['%.1f'%M]

   def calc_Pk(self, modellimits, pixdx, mc):
      #print modellimits
      npix = N.around((modellimits[1]-modellimits[0])/pixdx).astype('int')
      self.npix_model = npix
      mag_arr = N.arange(modellimits[0],modellimits[1],pixdx)
      self.Pk = N.zeros((len(self.zarr), npix))
      for i in range(len(self.zarr)):
         z = self.zarr[i]
         dk = mc(z) 
         for j in range(len(self.M1500arr)):
            M0 = self.M1500arr[j]
            M1 = M0 + self.dM
            if (M0 + dk) >= modellimits[1]:
               pass
            elif (M1 + dk) < modellimits[0]:
               pass
            else:
               k0 = N.searchsorted(mag_arr, M0+dk)
               k1 = N.searchsorted(mag_arr, M1+dk)
               self.Pk[i][k0:k1] = self.zdist['%.1f'%M0][i]*self.dVdz[i]

   def write(self, outname):
      zdist.write_zdgrid(self, outname)


class zdgrid_idrops(zdist.zdgrid,object):
   """
   A special class for i-dropouts.
   Selected with i-z v.s. z-J colors.
   """
   def __init__(self, M0, M1, dM, logR0, logR1, dlogR, z0, z1, dz, area):
   	zdist.zdgrid.__init__(self, M0, M1, dM, logR0, logR1, dlogR, z0, z1, dz, 
   	                      'i', area)

      
   def __call__(self,simcatalog,field,plot_diag=False,ktest=None,
                bands=def_bands, key='Hua13',
                sn_lolim={'wfc3_f160w':5.0,'wfc3_f125w':3.5,'acs_f850lp':2.0},
                sn_hilim={'acs_f435w':2.0},
                interpolate=False,dznew=0.1,expand=[0,0]):
      """
      simcatalog --- the simulation catalog as a FITS table.
      """
      # Initialize i-dropout color criteria
      print "dznew=", dznew
      print "interpolate?", interpolate
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc('f775w_drop', key)
      c = Ftable(simcatalog)
      # Now initialize the attributes
      for b in bands:
         bs = filter_dic[b]  # short name of filter b
         flux_iso = getattr(c,'%s_flux_iso'%bs)
         fluxerr_iso = getattr(c,'%s_fluxerr_iso'%bs)
         mag_iso = getattr(c,'%s_mag_iso'%bs)
         # calculate 1-sigma upper limits on sources with S/N < 1
         mag_iso = N.where(flux_iso/fluxerr_iso>1.0,mag_iso,
                           zeropoints[b]-2.5*N.log10(fluxerr_iso))
         setattr(self,b+'_mag',mag_iso)  
         # MAG_ISO --> for color calculation
         setattr(self,b+'_flux',getattr(c,'%s_flux_auto'%bs))  
         # FLUX_AUTO
         setattr(self,b+'_flux_aper',getattr(c,'%s_flux_aper_2'%bs))
         # Aperture flux within 0.18" aperture (2.94 pixels)
         setattr(self,b+'_fluxerr',getattr(c,'%s_fluxerr_auto'%bs))  
         # FLUXERR_AUTO
         setattr(self,b+'_fluxerr_aper',getattr(c,'%s_fluxerr_aper_2'%bs))
         # Aperture flux errors within 0.18" aperture
         setattr(self,b+'_sn',
            getattr(self,'%s_flux'%b)/getattr(self,'%s_fluxerr'%b))  
         # S/N calculated using AUTO APERTURE 
         setattr(self,b+'_sn_aper',
            getattr(self,'%s_flux_aper'%b)/getattr(self,'%s_fluxerr_aper'%b))
         # S/N calculated using 0.18" APERTURE
      # Now construct S/N criteria
      #self.sncrit = N.ones(len(c.d),'bool') 
      self.sn_lolim_crit = N.ones(len(c.d), 'bool')
      self.sn_hilim_crit = N.ones(len(c.d), 'bool') 
      for b in sn_lolim.keys():  # enforce S/N lower limits
         self.sn_lolim_crit = self.sn_lolim_crit & \
               (getattr(self,b+'_sn')>=sn_lolim[b])
      for b in sn_hilim.keys():  
         # enforce S/N upper limits (veto bands)
         # but use fixed aperture S/N
         self.sn_hilim_crit = self.sn_hilim_crit & \
               (getattr(self,b+'_sn_aper')<sn_hilim[b])
      print "Total number of objects satisfying the S/N criteria:", \
         sum((self.sn_hilim_crit)&(self.sn_lolim_crit))
      self.sncrit = (self.sn_hilim_crit==True) & (self.sn_lolim_crit==True)
      print "Do selections..."
      self.color1 = self.acs_f775w_mag-self.acs_f850lp_mag
      self.color2 = self.acs_f850lp_mag-self.wfc3_f125w_mag
      self.lcc.select(self.color1,self.color2)  # do color selection!
      self.colorcrit = self.lcc.crit.copy()
      self.lcc.crit = self.lcc.crit & self.sncrit  # enforce S/N criteria
      self.dropcrit = self.lcc.crit & (c.detect==True)  # just in case
      self.detect = c.detect.copy()
      print "Selection done."
      print "Total number of objects in the catalog: %d" % len(c.d)
      print "Total number selected as i-dropouts: %d" % (sum(self.dropcrit))

      # Now invoke zdist.zdgrid.__call__ to calculate P(z)
      zdist.zdgrid.__call__(self,c,interpolate=interpolate,dznew=dznew,
                            plot_diag=plot_diag,ktest=ktest,expand=expand)

   def write(self,outname):
   	zdist.write_zdgrid(self,outname)



