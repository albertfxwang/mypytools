#!/usr/bin/env python

from numpy import *
from pygoods import *
import bivRL as bl
import mlutil
import fit_lbg as fl
import bivariate_fit as bf
import time


def reduced_chi2(par, drop, mbins1, rbins1, mbins2, rbins2, chisqnulim=(0.4,5.0),
   zlo=-1):
   """
   Calculate the reduced chi2 of the best-fit model
   mbins, rbins should include the upper limits.
   """
   if drop == 'b':
      cat1 = 'bdrops_gf_v2.cat'
      cat2 = 'bdrops_udf_gf_v2.cat'
      kgrid1 = mlutil.readkgrid('kernel_I.p')
      kgrid2 = mlutil.readkgrid('kernel_I_udf.p')
      zdf1 = 'zdgrid_bdrops.p'
      zdf2 = 'zdgrid_bdrops_udf.p'
      mcfile = 'M1500_to_i.txt'
      mc = bl.mconvert(mcfile)
      mk = mc(4.0)
   elif drop == 'v':
      cat1 = 'vdrops_gf_v2.cat'
      cat2 = 'vdrops_udf_gf_v2.cat'
      kgrid1 = mlutil.readkgrid('kernel_Z.p')
      kgrid2 = mlutil.readkgrid('kernel_Z_udf.p')
      zdf1 = 'zdgrid_vdrops.p'
      zdf2 = 'zdgrid_vdrops_udf.p'
      mcfile = 'M1500_to_z.txt'
      mc = bl.mconvert(mcfile)
      mk = mc(5.0)
   cullflags = [0,1,2,3,4,12,13,14,18,19]
   limits1 = array([[21.0,26.5],[-2.0,3.0]])
   limits2 = array([[23.0,28.5],[-2.0,3.0]])
   pixdx = array([0.02, 0.02])
   modshape1 = (limits1[:,1]-limits1[:,0])/pixdx; modshape1=around(modshape1).astype('int')
   modshape2 = (limits2[:,1]-limits2[:,0])/pixdx; modshape2=around(modshape2).astype('int')
   # bin the points & tally the counts
   mag1, re1, crit1 = cleandata(cat1, chisqnulim=chisqnulim[0], magautolim=26.5,
      cullflags=cullflags, limits=limits1, zlo=zlo)
   mag2, re2, crit2 = cleandata(cat2, chisqnulim=chisqnulim[1], magautolim=28.5,
      cullflags=cullflags, limits=limits2, zlo=zlo)
   bincounts1 = histogram2d(mag1, log10(re1), bins=[mbins1, rbins1])[0].astype('float')
   bincounts2 = histogram2d(mag2, log10(re2), bins=[mbins2, rbins2])[0].astype('float')
   #print bincounts1, bincounts2

   # calculate the best-fit models
   model1 = bl.bivariate_lf(par, limits1, pixdx, kgrid=kgrid1, zdgridfile=zdf1, \
      mcfile=mcfile, drop=drop, field='goods', meankcorr=mk, add_interloper=True)
   model2 = bl.bivariate_lf(par, limits2, pixdx, kgrid=kgrid2, zdgridfile=zdf2, \
      mcfile=mcfile, drop=drop, field='udf', meankcorr=mk, add_interloper=True)
   phistar_mod = phistar(par, drop, zlo=zlo)
   model1.model = phistar_mod * model1.model
   model2.model = phistar_mod * model2.model
   #model1.model = ones(modshape1)/(modshape1[0]*modshape1[1]*pixdx[0]*pixdx[1])*len(mag1)
   #model2.model = ones(modshape2)/(modshape2[0]*modshape2[1]*pixdx[0]*pixdx[1])*len(mag2)
   print sum(model1.model.ravel())*pixdx[0]*pixdx[1]
   
   chi2tot = 0.
   nbins = 0

   mindex1 = (mbins1 - 21.0) / 0.02; mindex1 = around(mindex1).astype('int')
   rindex1 = (rbins1 - (-2.0)) / 0.02; rindex1 = around(rindex1).astype('int')
   mindex2 = (mbins2 - 23.0) / 0.02; mindex2 = around(mindex2).astype('int')
   rindex2 = (rbins2 - (-2.0)) / 0.02; rindex2 = around(rindex2).astype('int')

   num_exp1 = []   # number of expected 
   num_exp2 = []
   num_obs1 = bincounts1.ravel()[bincounts1.ravel()>=5]
   num_obs2 = bincounts2.ravel()[bincounts2.ravel()>=5]
   # iterate through all bins and calculate the chi2
   for i in range(len(mbins1)-1):
      for j in range(len(rbins1)-1):
         if bincounts1[i,j] >= 5:
            num_mod = sum(model1.model[mindex1[i]:mindex1[i+1],rindex1[j]:rindex1[j+1]].ravel())
            num_mod = num_mod * pixdx[0] * pixdx[1]
            num_exp1 += [num_mod]
            #chi2 = (bincounts1[i,j] - num_mod)**2 / num_mod
            #print bincounts1[i,j], num_mod
            #chi2tot += chi2
            #nbins += 1
            #if bincounts1[i,j] < nmin: nmin = bincounts1[i,j]
   for i in range(len(mbins2)-1):
      for j in range(len(rbins2)-1):
         if bincounts2[i,j] >= 5:
            num_mod = sum(model2.model[mindex2[i]:mindex2[i+1],rindex2[j]:rindex2[j+1]].ravel())
            num_mod = num_mod * pixdx[0] * pixdx[1]
            num_exp2 += [num_mod]
            #chi2 = (bincounts2[i,j] - num_mod)**2 / num_mod
            #chi2tot += chi2
            #nbins += 1
   print "nbins", nbins
   
   
   # Run chi-square test
   num_exp = concatenate([num_exp1, num_exp2])
   num_obs = concatenate([num_obs1, num_obs2])
   ndeg = len(num_exp) - 1   # degree of freedom = num. of contributing bins - 1
   print "ndeg", ndeg
   chi2, pval = stats.mstats.chisquare(num_obs, f_exp=num_exp)
   print chi2, pval
   #print num_exp, num_obs
   #chi2nu = chi2tot / float(ndeg)
   return chi2, pval, num_exp, num_obs


def mc_goodness_fit(par, drop, niter, zlo=-1.):
   # Determine goodness of fit given the parameters.
   # Draw the same number of points from the model as observed, then calculate the 
   # loglikelihood of the drawn points. Calculate the probability that a simulated 
   # observation has lower loglikelihood than the observed points.
   if drop == 'b':
      cat1 = 'bdrops_gf_v2.cat'
      cat2 = 'bdrops_udf_gf_v2.cat'
      kgrid1 = mlutil.readkgrid('kernel_I.p')
      kgrid2 = mlutil.readkgrid('kernel_I_udf.p')
      zdf1 = 'zdgrid_bdrops.p'
      zdf2 = 'zdgrid_bdrops_udf.p'
      mcfile = 'M1500_to_i.txt'
      mc = bl.mconvert(mcfile)
      mk = mc(4.0)
      chisqnulim = [0.4, 5.0]
   elif drop == 'v':
      cat1 = 'vdrops_gf_v2.cat'
      cat2 = 'vdrops_udf_gf_v2.cat'
      kgrid1 = mlutil.readkgrid('kernel_Z.p')
      kgrid2 = mlutil.readkgrid('kernel_Z_udf.p')
      zdf1 = 'zdgrid_vdrops.p'
      zdf2 = 'zdgrid_vdrops_udf.p'
      mcfile = 'M1500_to_z.txt'
      mc = bl.mconvert(mcfile)
      mk = mc(5.0)
      chisqnulim = [0.5, 5.0]
   cullflags = [0,1,2,3,4,12,13,14,18]
   limits1 = array([[21.0,26.5],[-2.0,3.0]])
   limits2 = array([[23.0,28.5],[-2.0,3.0]])
   #limits1 = bl.limits1
   #limits2 = bl.limits2
   pixdx = array([0.02, 0.02])
   mag1, re1, crit1 = fl.cleandata(cat1, chisqnulim=chisqnulim[0], magautolim=26.5,
      cullflags=cullflags, limits=limits1, zlo=zlo, drop=drop)
   mag2, re2, crit2 = fl.cleandata(cat2, chisqnulim=chisqnulim[1], magautolim=28.5,
      cullflags=cullflags, limits=limits2, zlo=zlo, drop=drop)
   data1 = array([mag1, log10(re1)])
   data2 = array([mag2, log10(re2)])
   N1 = len(mag1)
   N2 = len(mag2)
   print N1, N2, N1+N2
   model1 = bl.bivariate_lf(par, limits1, pixdx, drop, 'goods', kgrid=kgrid1, zdgridfile=zdf1,\
      mcfile=mcfile, meankcorr=mk, add_interloper=True, norm=-1.)
   model2 = bl.bivariate_lf(par, limits2, pixdx, drop, 'udf', kgrid=kgrid2, zdgridfile=zdf2,\
      mcfile=mcfile, meankcorr=mk, add_interloper=True, norm=-1.)
   sum1 = sum(model1.model.ravel()) * pixdx[0] * pixdx[1]
   sum2 = sum(model2.model.ravel()) * pixdx[0] * pixdx[1]
   phistar_mod = float(N1+N2)/(sum1 + sum2)
   print phistar_mod
   model1.model = phistar_mod * model1.model
   model2.model = phistar_mod * model2.model
   logl_ref = bf.loglikelihood(data1,model1,floor=0.)+bf.loglikelihood(data2,model2,floor=0.)
   #logl_ref = bf.mlfunc(par, data1, data2, limits1, limits2, pixdx, kgrid1, kgrid2,
   #   1.0, 1.0, -21.0, zdf1, zdf2, mcfile, 1, mk, 0, -1., 'phistar', drop)
   print "logl_ref", logl_ref
   simlogl_arr = zeros(niter)  # actually -1*logL...
   print "Start drawing simulated observations..."
   t1 = time.time()
   for i in range(niter):
      if i % 1000 == 0: print i
      simdata1 = mlutil.draw_from_pdf(N1, model1, model1.limits)
      simdata2 = mlutil.draw_from_pdf(N2, model2, model2.limits)
      simlogl1 = bf.loglikelihood(simdata1, model1)
      simlogl2 = bf.loglikelihood(simdata2, model2)
      simlogl = simlogl1 + simlogl2
      simlogl_arr[i] = simlogl
   t2 = time.time()
   dt = t2 - t1
   dtmin = int(floor(dt)) / 60
   dtsec = dt % 60
   n_worse = sum(simlogl_arr > logl_ref)
   print "%d iterations took %d min %.1f sec" % (niter, dtmin, dtsec)
   print "Percentage of simulated observations with lower log-likelihood: %.2f %%" % (
      100.*float(n_worse)/float(niter))
   return logl_ref, simlogl_arr

