#!/usr/bin/env python

from numpy import *
from pygoods import *
import os
import scipy

def draw_indices(N,indexarray):
   """draw N galaxies randomly from indexarray
      len(indexarray) must be greater than N"""
   output_index = zeros(N,'int')
   for i in range(N):
      ntot = len(indexarray)

      # draw a random number from index
      x = random.uniform(-0.5,ntot-0.5)   
      x = int(round(x))
      output_index[i] = indexarray[x]
      delete(indexarray,x)  # remove x-th element from indexarray

   return output_index

def read_ndist(ndistfile):
   c = sextractor(ndistfile)
   sersicn = c._1
   number = c._2
   return sersicn, number

def draw_ndist(diskfrac,ntot,gtype_array,verbose=0,disk=0):
   """gtype == disk (default is 1) if disk, otherwise devauc
      draw a total of ntot galaxies (ntot*diskfrac of them
      are disks) and return their indices from gtype_array
   """

   indexarray = arange(len(gtype_array))
   if verbose:
      print "In input, number of disks = ",sum(gtype_array==disk)
      print "In input, number of devauc = ",sum(gtype_array!=disk)
   disk_indexarray = compress(gtype_array==disk,indexarray)
   devauc_indexarray = compress(gtype_array!=disk,indexarray)

   # calculate number of each type to be drawn
   ndisk = int(round(ntot * diskfrac))  # number of disks to be drawn
   ndevauc = ntot - ndisk  # number of devauc objects
   if (ndisk>len(disk_indexarray)) | (ndevauc>len(devauc_indexarray)):
      raise ValueError, "Too many disk/devauc objects to be drawn"
   if verbose:
      print "ndisk, ndevauc = ",ndisk,ndevauc

   # draw indices from respective index array
   disk_outindex = draw_indices(ndisk,disk_indexarray)
   devauc_outindex = draw_indices(ndevauc,devauc_indexarray)

   # joins two index arrays
   outindex = concatenate([disk_outindex,devauc_outindex])
   return outindex

def sim_ndist(diskfrac,ntot,gtype,nout,magerr,reout,reerr,recovered,
   chisqnu,dmag=1.0,reerrrmax=1.0,chinu_max=0.45,nbins=arange(0.,9.),verbose=0):
   """gtype == 1 if devauc or 0 if disk
      draw a total of ntot galaxies (ntot*diskfrac of them
      are disks) and return their n-distribution"""
   # draw points 
   outindex = draw_ndist(diskfrac,ntot,gtype,nout,nbins=nbins,verbose=verbose)

   # Retain only points recovered by GALFIT & qflags==1
   nout = nout.take(outindex)
   gtype = gtype.take(outindex)
   recovered = recovered.take(outindex)
   magerr = magerr.take(outindex)
   reout = reout.take(outindex)
   reerr = reerr.take(outindex)
   chisqnu = chisqnu.take(outindex)
   qflags = (magerr<=dmag) * ((reerr/reout)<=reerrrmax) * (chisqnu<=chinu_max)
   qflags = qflags * (reout>0) * (reerr!=inf)
   crit = recovered * qflags
    
   # read output Sersic n
   nout_drawn = compress(crit,nout)
   ndist_out = histogram(nout_drawn,bins=nbins)
   gtype_drawn = compress(crit,gtype)
   diskfrac_drawn = float(sum(gtype_drawn==0))/float(len(gtype_drawn))
   outindex = compress(crit,outindex)
   print "Final number of drawn objects: %d, final diskfrac = %.1f%%"%(len(outindex),
     diskfrac_drawn*100.)

   # return both the drawn n array and histogram
   return nout_drawn,ndist_out,outindex
   
def diffsq_ndist(diskfrac,gtype,nout,magerr,reout,reerr,recovered,ndistfile,
   chisqnu,dmag=1.0,reerrrmax=1.0,chinu_max=0.45,ntot=50000):
   """Calculates the squared differences between drawn n-dist and 
   real LBG n-dist"""
   nbins,ndist_real = read_ndist(ndistfile)
   totn_real = sum(ndist_real)
   nbins = append(nbins,nbins[-1]+(nbins[1]-nbins[0]))
   
   # draw ndist with diskfrac
   nout_drawn,ndist_out,outindex = sim_ndist(diskfrac,ntot,gtype,nout,
      magerr,reout,reerr,recovered,chisqnu,dmag=dmag,reerrrmax=reerrrmax,
      chinu_max=chinu_max,nbins=nbins)

   # normalize ndist_out to total number of ndist_real
   ndist_out_norm = ndist_out[0]*float(totn_real)/float(sum(ndist_out[0]))

   # calculate squared difference between 2 distributions
   diffsq = sum((ndist_out_norm-ndist_real)**2)
   return diffsq

def fit_ndist(diskfrac_guess,ar_gtype,ar_nout,ndistfile,ntot):
   args = [ar_gtype,ar_nout,ndistfile,ntot]
   diskfrac_best,others = scipy.optimize.fmin(diffsq_ndist,diskfrac_guess,args=args)
   print "best-match diskfrac to %s:"%ndistfile,diskfrac_best
   
