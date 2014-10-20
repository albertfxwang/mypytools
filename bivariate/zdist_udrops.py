#!/usr/bin/env python

import numpy as N
from pygoods import sextractor, Ftable
import cPickle
from bivariate import zdist
from bivariate import udropsim_uflux as uu
import matplotlib.pyplot as plt
from dropout_selection import lbg_colorcrit as lcc

filter_dic = {
   'vimos_u':'u',
   'acs_f435w':'b',
   'acs_f606w':'v',
   'acs_f775w':'i',
   'acs_f850lp':'z',
   'wfc3_f125w':'j',
   'wfc3_f160w':'h',
   'wfc3_f105w':'y',
   'wfc3_f098m':'y'}
zeropoints = {
   'wfc3_f160w':25.960,
   'acs_f435w':25.673,
   'acs_f606w':26.505,
   'acs_f775w':25.678,
   'acs_f850lp':24.867,
   'wfc3_f125w':26.250,
   'vimos_u':26.158,
   'wfc3_f105w':26.270,
   'wfc3_f098m':26.270}

class zdgrid_udrops(zdist.zdgrid, object):
   """
   A special class for U-dropouts, because it uses U-band photometry.
   Selected with U-B v.s. B-Y colors.
   """
   def __init__(self, M0, M1, dM, logR0, logR1, dlogR, z0, z1, dz, area, 
                umag_kpdf, umag_1sig_kpdf):
      zdist.zdgrid.__init__(self, M0, M1, dM, logR0, logR1, dlogR, z0, z1, 
                            dz, 'u', area)
      self.umag_kpdf = umag_kpdf
      self.umag_1sig_kpdf = umag_1sig_kpdf  
      # the file containing the PDF of the 1-sigma magnitude limits

   def __call__(self, simcatalog, field, plot_diag=False,ktest=None,n_repeat=1,
      sn_lolim={'wfc3_f160w':5.0,'wfc3_f105w':5.0,'acs_f435w':3.0},
      mode='ubvy', interpolate=False, dznew=0.1, redraw=False,
      drawnflux=None, testumag=False, testdraw=False, mag0=22.0,
      expand=[0.,0.], mag_in_col='m1500_in', re_in_col='re_in'):
      """
      simcatalog --- the simulation catalog as a FITS table.
      Because the U-band magnitudes are randomly drawn from each object, we 
      would like to repeat the magnitude drawing process for a few times to 
      average out any statistical anomalies. 
      mode --- the bands to use for selection (uby v.s. ubvy)
      """
      # Initialize U-dropout color criteria
      print "n_repeat", n_repeat
      print "dznew", dznew
      print "interpolate?", interpolate
      print "Redraw?", redraw
      self.lcc = lcc.colorcrit()
      if mode=='uby':
         if field=='ers':
            self.lcc = self.lcc('uvimos_drop','UBY098')
         else:
            self.lcc = self.lcc('uvimos_drop','UBY105')
      elif mode=='ubvy':
         if field=='ers':
            self.lcc = self.lcc('uvimos_drop','UBVY098')
         else:
            self.lcc = self.lcc('uvimos_drop','UBVY105')
      c = Ftable(simcatalog)
      #self.c = c
      ## First: draw U-band magnitudes, and initialize correct attributes of c
      ## Tile the arrays by `n_repeat` times
      if field=='ers':
         bands=['acs_f435w','acs_f606w','wfc3_f098m','wfc3_f160w']  
         # HST bands only
      else:
         bands=['acs_f435w','acs_f606w','wfc3_f105w','wfc3_f160w']
      for b in bands:
         bs = filter_dic[b]  # short name of filter b
         flux_iso = getattr(c,'%s_flux_iso'%bs)
         fluxerr_iso = getattr(c,'%s_fluxerr_iso'%bs)
         mag_iso = getattr(c,'%s_mag_iso'%bs)
         # calculate 1-sigma upper limits on sources with S/N < 1
         mag_iso = N.where((flux_iso / fluxerr_iso) > 1.0, mag_iso, 
                           zeropoints[b] - 2.5 * N.log10(fluxerr_iso))
         setattr(self,b + '_mag', N.tile(mag_iso, n_repeat))  
         # MAG_ISO --> for color calculation
         # setattr(self,b + '_flux', 
   #               N.tile(getattr(c, '%s_flux_auto'%bs),n_repeat))  
         setattr(self,b + '_flux', 
                 N.tile(getattr(c, '%s_flux_iso'%bs),n_repeat))
         # FLUX_AUTO
         # setattr(self, b + '_fluxerr', 
   #               N.tile(getattr(c,'%s_fluxerr_auto'%bs),n_repeat))  
         setattr(self, b + '_fluxerr', 
                 N.tile(getattr(c,'%s_fluxerr_iso'%bs),n_repeat)) 
         # FLUXERR_AUTO
         setattr(self, b + '_sn',
                 getattr(self, '%s_flux'%b) / getattr(self, '%s_fluxerr'%b))  
         # S/N calculated using FLUX_ISO

      print "Step 1: draw U-band magnitudes and 1-sigma limits from files \
            %s and %s" % (self.umag_kpdf,self.umag_1sig_kpdf)
      # Draw U-band magnitude
      # Calculate the fraction that have U-band S/N >= 1
      if field == 'udf':
         cu = Ftable('/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/udrops_fitting/simcatalogs/run1_udf_130625.fits')
         umag_bins = N.arange(25., 38.5, 0.5)
      else:
         cu = Ftable('/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/udrops_fitting/simcatalogs/run2_tfitsim_130322.fits')
         umag_bins = N.arange(25., 38.5, 0.5)
      umag_in = cu.uvimos_mag_in[cu.uvimos_fitquerr>0.]
      u_ston = (cu.uvimos_fitqty/cu.uvimos_fitquerr)[cu.uvimos_fitquerr>0.]
      n1 = N.histogram(umag_in[u_ston>=1.], umag_bins)[0]
      n2 = N.histogram(umag_in, umag_bins)[0]
      u_detect_frac = n1.astype('float') / N.maximum(1., n2.astype('float'))

      if testumag==False:
         # If not testing anything...
         if redraw == True:
            print "Repeat the drawing process %d times" % n_repeat
            self.vimos_u_mag = uu.draw_umag(N.tile(c.u_mag_in, n_repeat),
                                      self.umag_kpdf, mag0=mag0)  
            self.vimos_u_mag_1sig = uu.draw_ulimmag(n_repeat*len(c.u_mag_in),
                                          self.umag_1sig_kpdf)
            # Draw 1-sigma magnitude upper limit
            # Now replace U-band magnitude with 1-sigma upper limit by using 
            # the S/N > 1 fraction
            # index_detfrac = (c.u_mag_in - c.u_mag_in.min()) - (c.u_mag_in % 0.5)
            # index_detfrac = N.minimum(len(u_detect_frac)-1, index_detfrac)
            # detfrac = u_detect_frac.take(index_detfrac.astype('int'))
            # # Now randomly draws a reference level to compare with the detection
            # # fraction
            # rdn_level = N.random.uniform(0., 1., size=len(c.d))
            # self.vimos_u_mag = N.where(detfrac>=rdn_level, self.vimos_u_mag,
            #                            self.vimos_u_mag_1sig)
            self.vimos_u_mag = N.minimum(self.vimos_u_mag,
                                        self.vimos_u_mag_1sig)
         else:
            print "Read U-band magnitudes from %s." % drawnflux
            # Read the drawn U-band magnitudes from file drawnflux
            df = cPickle.load(open(drawnflux))
            self.vimos_u_mag_1sig = df.vimos_u_mag_1sig
            self.vimos_u_mag = df.vimos_u_mag
      else:
         self.vimos_u_mag = c.u_mag_in
         self.vimos_u_mag_1sig = c.u_mag_in
      # Now construct S/N criteria
      #self.sncrit = N.ones(len(c.d)*n_repeat,'bool')
      self.sn_lolim_crit = N.ones(len(c.d)*n_repeat, 'bool')
      self.detect = c.detect.copy()
      #print len(self.sncrit), len(self.acs_f435w_sn)
      if len(sn_lolim.keys()) > 0:
         print "sn_lolim", sn_lolim
         for b in sn_lolim.keys():
            self.sn_lolim_crit = self.sn_lolim_crit & \
               (getattr(self,b+'_sn')>=sn_lolim[b])
            #self.sncrit = self.sncrit & (getattr(self,'%s_sn'%b)>=sn_lolim[b])
      self.sncrit = (self.sn_lolim_crit==True)
      ## Second: calculate dropcrit (do color selection)
      #import select_udrops_uby as suu
      print "Do selections..."
      if field=='ers':
         self.color1 = self.vimos_u_mag-self.acs_f435w_mag
         if mode=='uby':
            self.color2 = self.acs_f435w_mag-self.wfc3_f098m_mag
         elif mode=='ubvy':
            self.color2 = self.acs_f606w_mag-self.wfc3_f098m_mag
         self.lcc.select(self.color1,self.color2)

      else:
         self.color1 = self.vimos_u_mag-self.acs_f435w_mag
         if mode=='uby':
            self.color2 = self.acs_f435w_mag-self.wfc3_f105w_mag
         elif mode=='ubvy':
            self.color2 = self.acs_f606w_mag-self.wfc3_f105w_mag
         self.lcc.select(self.color1,self.color2)
      self.colorcrit = self.lcc.crit.copy()
      self.lcc.crit = self.lcc.crit & self.sncrit  # Fold in S/N criteria
      #if field=='ers':
      #   self.udrops = suu.udrops_ubvy098(self.c,
      #      bands=['acs_f435w','acs_f606w','wfc3_f098m','wfc3_f160w'],
      #      sn_lolim={'wfc3_f160w':5.0,'acs_f435w':3.0,'wfc3_f098m':5.0})
      #else:
      #   self.udrops = suu.udrops_ubvy105(self.c,
      #      bands=['acs_f435w','acs_f606w','wfc3_f105w','wfc3_f160w'],
      #      sn_lolim={'wfc3_f160w':5.0,'acs_f435w':3.0,'wfc3_f105w':5.0})
      #dropcrit = self.udrops.crit   # the dropout-selection array
      dropcrit = self.lcc.crit & (N.tile(c.detect,n_repeat)==True)
      #self.dropcrit = dropcrit[:len(c.d)]
      self.dropcrit = dropcrit
      self.detect = c.detect.copy()
      print "Selections done."
      print "Total number of objects in the catalog: %d" % len(c.d)
      print "Total number selected as U-dropouts: %d" % (sum(dropcrit)/float(n_repeat))
      ## Third: call make_Pz to calculate P(z)
      
      if not testdraw:
         # Now invoke zdist.zdgrid.__call__ to calculate P(z)
         zdist.zdgrid.__call__(self,c,interpolate=interpolate,dznew=dznew,
                            plot_diag=plot_diag,ktest=ktest,n_repeat=n_repeat,
                            expand=expand, mag_in_col=mag_in_col,
                            re_in_col=re_in_col)
      if self.n_repeat > 30:
         delattr(self, 'color1')
         delattr(self, 'color2')
         delattr(self, 'sncrit')


   def write(self,outname):
      self.filename = outname
      zdist.write_zdgrid(self,outname)

