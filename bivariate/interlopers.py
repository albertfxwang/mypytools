#!/usr/bin/env python

from numpy import *
from pygoods import *
import KPDFadaptnumpy as KPDF
import cPickle
from transfunc.bivmodel import bivmodel

# calculate the contribution of interlopers to the total number counts
# add the contribution to the "clean" bivariate model to get the expected number counts
# this is an effort to take into account the interlopers statistically


interloper_num_goods = {'v':[0.,0.,0.,0.,0.,0.0263,0.0614,0.0394,0.0468,0.061,0.079],
                         'b':[0.,0.,0.,0.,0.0076,0.0032,0.0207,0.0216,0.0276,0.0498,0.0728]}
interloper_mag_goods = arange(21., 26.5, 0.5)
interloper_frac_udf = {'v':[0.,0.,0.,0.,0.0988,0.0291,0.0546,0.0179,0.0206,0.0394,0.0542],
                       'b':[0.,0.,0.0087,0.,0.0564,0.0242,0.0509,0.0503,0.0537,0.0822,0.1164]}
interloper_mag_udf = arange(23.0, 28.5, 0.5)
interloper_frac_goods_all = {'v':0.0707, 'b':0.0550}
interloper_frac_udf_all = {'v':0.0377, 'b':0.076}


#def create_interloper_frac(drop, field):
#   # return the fraction of interlopers from simulations
#   #simdir = '/Users/khuang/dropbox/research/bivariate/phot_scatter_sim'
#   if field == 'goods':
#      if drop == 'b':
#         c = sextractor('bdrops_goods_interloper.txt')
#      elif drop == 'v':
#         c = sextractor('vdrops_goods_interloper.txt')
#      mbins = c._1
#      n_drops = c._6
#      frac_int = c._4 / maximum(n_drops, 1)
#      place(frac_int, n_drops<1, 0.)
#      place(frac_int, mbins<24.5, 0.)  
#      # only consider interloper fractions fainter than 24.5
#      return mbins, frac_int
#   elif field == 'udf':
#      if drop == 'b':
#         cu = sextractor('bdrops_udf_interloper.txt')
#      elif drop == 'v':
#         cu = sextractor('vdrops_udf_interloper.txt')
#      mbins = cu._1
#      n_drops = cu._6
#      frac_int = cu._4 / maximum(n_drops, 1)
#      place(frac_int, n_drops<1, 0.)
#      place(frac_int, mbins<25.5, 0.)
#      return mbins, frac_int
#   else:
#      print "No interloper data exists."
#      return 1

class sizedist(object): 
   def __init__(self, mbins, logrbins=arange(-2.,3.,0.2)):
      self.mbins = mbins
      self.logrbins = logrbins
      #self.field = field
      # First, filter out unreasonable sizes
      
      #reout_err = cs[0].reout_err
      #if len(cs) > 1:
      #   for i in range(1, len(cs)):
      #      magauto = concatenate((magauto, cs[i].mag_auto))
      #      reout = concatenate((reout, cs[i].reout))
      #      reout_err = concatenate((reout_err, cs[i].reout_err))
      self.dm = mbins[1] - mbins[0]
      self.dlogr = logrbins[1] - logrbins[0]
      self.sd = zeros([len(mbins),len(logrbins)])

   def __call__(self, mag_array, size_array, h=None):
      # mag_array: array of all magnitudes in the catalog cs
      # size_array: array of all sizes in the catalog cs (in pixels)
      size_array = size_array[size_array>0]
      mag_array = mag_array[size_array>0]
      logr_array = log10(size_array)
      logr_bins = concatenate([self.logrbins, [self.logrbins[-1]+self.dlogr]])
      #print "logr_bins", logr_bins
      #qflag = (reout_err / reout <= 0.6)
      for i in range(len(self.mbins)):
         #re = compress(qflag & (magauto>=mbins[i]) & (magauto<mbins[i]+dm), reout)
         #logre = log10(re)
         bincrit = (mag_array>=self.mbins[i])&\
                   (mag_array<(self.mbins[i]+self.dm))
         logr_thisbin = logr_array[bincrit==True]
         #print "sum(bincrit)", sum(bincrit)
         if h==None:
            self.sd[i] = histogram(logr_thisbin, logr_bins)[0]
            print "len(logr_thisbin)", len(logr_thisbin)
            self.sd[i] = self.sd[i] / sum(self.sd[i])
         else: 
            if h==-1:
               h = KPDF.UPDFOptimumBandwidth(logr_thisbin)
            p = KPDF.UPDFEpanechnikov(logr_thisbin, self.logrbins, h)
            self.sd[i] = p / sum(p)  
         # normalize to 1.0
         #mname = "%.1f" % mbins[i]
         #sddict[mname] = sd
      #sddict['lograrr'] = lograrr
      

#def interloper_model(sdfile, limits, drop, field, model0, pixdx=array([0.02,0.02])):
class interloper_model(bivmodel):
   # construct a component of interloper number density on the desired model grid
   # model0: the bivariate model distribution of dropouts only (calculated from the parameters)
   # Strategy: we have the interloper fractions in magnitude bins, so we normalize the
   # interloper contribution to the total number counts according to the fractions and the
   # number counts given by model0. So model0 provides a standard of normalization for the
   # interloper model. The input model0, on the other hand, has the unit of number **density**
   # so it already includes the dropout selection effective volume.
   # N: the total number of dropouts in the given field
   def __init__(self, limits, model0, intfrac, sd, 
                pixdx=array([0.02,0.02]), mag_lolim=23.0):
      """
      model0: a bivmodel instance of the uncorrected RL distribution.
      sdfile: a pickled file containing the sizedist class instance.
      mag_lolim: the bright limit over which do not include interlopers,
                 because it's easier to identify interlopers at the bright end.
      sd: a sizedist class instance for the size distribution
      intfrac: a FITS table with interloper fractions
      """
      npix = (limits[:,1] - limits[:,0]) / pixdx
      npix = around(npix).astype('int')
      intmodel = zeros(npix, 'float')
      bivmodel.__init__(self, limits, pixdx, intmodel)
      # Read interloper fractions (from a FITS table)
      #intfrac = Ftable(intfrac_file)
      # Read size distribution in the field in this filter
      #f = open(sdfile)
      #sd = cPickle.load(f)
      #f.close()
      #lograrr = sd['lograrr']
      #dlogr = sd.lograrr[1] - sd.lograrr[0]
      #mbin, frac_int = create_interloper_frac(drop, field)
      #print mbin, frac_int
   
      #dm = mbin[1] - mbin[0]
      # figure out the size distribution in each magnitude bin
      # the magnitude bin limits of the interloper numbers MUST match the bin limits
      # in the size distribution... can I write a more generalized way to make it not the case?
      #for i in range(len(sd.mbin)):
      dm = intfrac.mag0[1] - intfrac.mag0[0]
      dlogr = sd.dlogr
      for i in range(len(intfrac.mag0)):
         if intfrac.mag0[i] >= mag_lolim:
            mag_str = "%.1f" % intfrac.mag0[i]
            try:
               ii = map(lambda x:"%.1f"%x, sd.mbins).index(mag_str)
            except:
               print "%.1f not in sd.mbins" % intfrac.mag0[i]
               continue
            if (intfrac.mag0[i]+dm < self.limits[0,0]): continue
            if (intfrac.mag0[i] > self.limits[0,1]): continue
            #if frac_int[i] <= 1.e-5: continue
            #mkey = "%.1f" % mbin[i]
            #logrdist = sd[mkey] / (sum(sd[mkey])*pixdx[1])
            x0 = (intfrac.mag0[i]-self.limits[0,0])/self.pixdx[0] 
            x0 = int(round(x0))
            x1 = (intfrac.mag0[i]-self.limits[0,0]+sd.dm)/self.pixdx[0]; 
            x1 = int(round(x1))
            for j in range(len(sd.logrbins)):
               if (sd.logrbins[j]+dlogr < self.limits[1,0]): continue
               if (sd.logrbins[j] > self.limits[1,1]): continue
               # calculate the appropriate indexes
               y0 = (sd.logrbins[j]-self.limits[1,0])/self.pixdx[1]
               y0 = int(round(y0))
               y1 = (sd.logrbins[j]-self.limits[1,0]+dlogr)/self.pixdx[1]
               y1 = int(round(y1))
               # assign size distribution value to this mag-size bin
               self.model[x0:x1,y0:y1] = sd.sd[ii,j]
            # normalize the size distribution
            #intmodel[x0:x1] = intmodel[x0:x1] / (sum(intmodel[x0:x1].ravel()) * pixdx[1])
            # use interloper fraction and N_0 to calculate the expected numbers
            # where N_0 is the number density of REAL dropouts
            #print frac_int[i]
            N_0 = sum(model0.model[x0:x1].ravel())*self.pixdx[0]*self.pixdx[1]
            N_int = N_0 * intfrac.interloper_frac_drops[i] / \
                          (1. - intfrac.interloper_frac_drops[i])
            # normalize the interloper model in this magnitude bin to N_int
            self.model[x0:x1] = self.model[x0:x1] * N_int / \
                           (sum(self.model[x0:x1])*self.pixdx[0]*self.pixdx[1])
            #print N_0, N_int, sum(intmodel[x0:x1].ravel())*pixdx[0]*pixdx[1]
            # N_int = f_int * N_tot 
            # N_tot = N_0 + N_int = N_0 + f_int * N_tot 
            # => N_0 = N_tot * (1 - f_int) => N_tot = N_0 / (1 - f_int)
            # N_int = f_int * N_tot = N_0 * f_int / (1 - f_int)
            #intmodel[x0:x1,y0:y1] = factor * (sd[mkey][j]/sum(sd[mkey])) * \
            #   model0[x0:x1,y0:y1]

   #n_int_tot = sum(intmodel.ravel()) * pixdx[0] * pixdx[1]
   #n0 = sum(model0.ravel()) * pixdx[0] * pixdx[1]
   #print "total num of interlopers:", n_int_tot
   #print "total num of real dropouts:", n0
   #print n_int_tot / n0 
   #return intmodel
