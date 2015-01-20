#!/usr/bin/env python

import sextractor_sim, fake_galaxies
import time, os, sys, glob
import pyfits
import numpy as np
from PhotomTools import SkyEstimate
from pygoods import sextractor, angsep
import matplotlib.pyplot as plt
from stats import robust, gauss
import copy

class SkyRMSSim(sextractor_sim.SExtractorSim):
   """
   A class that reads a text file containing input galaxy positions instead of 
   randomly generate input positions, in order to calculate a median sky RMS.
   """
   def __init__(self, parfile):
      super(SkyRMSSim, self).__init__(parfile)
      self.posfile = self.c['POSFILE']
      self.skyaper = self.c['SKYAPER']
      self.fit_sky = self.c['FIT_SKY']

   def calc_input_xy(self, N, flagimage, output='inputxy.txt'):
      # calculate the input positions for sky RMS simulations. Avoid real 
      # sources.
      inputx = np.zeros(N)
      inputy = np.zeros(N)
      flagimg = pyfits.getdata(flagimage)
      i = 0
      while i < N:
         xn = np.random.uniform(0, flagimg.shape[1])
         yn = np.random.uniform(0, flagimg.shape[0])
         if not flagimg[yn,xn] > 0:
            inputx[i] = xn
            inputy[i] = yn
            i += 1
      f = open(output, 'wb')
      f.write('# 1 X\n')
      f.write('# 2 Y\n')
      for i in range(N):
         f.write('%.2f  %.2f  ' % (inputx[i], inputy[i]))
         f.write('\n')
      f.close()


   def insert_fake_sources(self):
      """
      Generate attributes for artificial galaxies, and then insert artificial 
      galaxies into real images.
      ** Only insert fake galaxies into the detection image. **
      """
      ### Need to be updated for multi-band, input-SED case
      # if no specified list of fake galaxies
      # make artdata file of input fake galaxies
      if glob.glob('*.list'):
         os.system('rm *.list')
      args = [self.realimages, self.flagimages, self.bands, self.posfile]
      kwargs = dict(ngal=self.ngal, diskfrac=self.diskfrac, magfile=self.magfile,
                    rdist=self.rdistfunc, mag0=self.maglow, mag1=self.maghigh,
                    logr0=self.logrmin, logr1=self.logrmax, lognormal_beta=self.lognormal_beta,
                    lognormal_mag0=self.lognormal_mag0, lognormal_peak=self.lognormal_peak,
                    lognormal_sigma=self.lognormal_sigma)
      self.fg = fake_galaxies.FakeGalaxiesCustomXY(*args, **kwargs)
      # Will update later to use fake_galaxies_sexsim!
      b = self.detect_band
      self.fg.makegals_multiband(self.flagimages[b], bands=[b])  
      # write input file for mkobjects
      self.makenoiselessimage(b, self.fg.artfiles[b], self.zeropoints[b],
                              self.psffiles[b])
      self.addsimulated(b)
      self.fake_detect_image = copy.copy(self.fakeimages[b])
      print "self.fake_detect_image", self.fakeimages[b]
      # Now zero-out the pixels outside of footprint...
      h = pyfits.open(self.fake_detect_image, mode='update')
      flg_detect = pyfits.getdata(self.flagimages[b])
      temp = np.where(flg_detect > 0, 0., h[0].data)
      h[0].data = temp
      h.flush()
      h.close()
      # now use images *without* fake sources as fake images
      for b in self.bands:
         self.fakeimages[b] = self.realimages[b]
         self.fg.artfiles[b] = self.fg.artfiles[self.detect_band]

   def fit_sky(self, band, niter, skyaper=0.4, width=64, plot=False):
      # Fit the local sky fluctuation...
      c = sextractor('%s_%s/run%d.cat' % (self.root, band, niter))
      print "number of regions to fit: ", len(c)
      if len(c) > 0:
         sky_sigma = np.zeros(len(c))
         sky_sigma_aper = np.zeros(len(c))
         # update segmentation map for sky fitting
         h = pyfits.open(self.segimages[band], mode='update')
         seg = h[0].data
         flg = pyfits.getdata(self.flagimages[band])
         h[0].data = np.where((seg==0)&(flg>0), seg.max()+1, seg)
         h.flush()
         h.close()
         for i in range(len(c)):
            skyest = SkyEstimate.SkyEstimate(self.realimages[band], 
                                             self.segimages[band], band)
            bkgd, sigma = skyest.skyest_xy(c.x_image[i], c.y_image[i], 
                                           width=width, plot=plot)
            sky_sigma[i] = sigma
            sky_sigma_aper[i] = skyest.calc_aperture_noise(skyaper)
         c.addcolumn('sky_sigma', sky_sigma, '%f')
         c.addcolumn('sky_sigma_aper', sky_sigma_aper, '%f')
         c.writeto(c._fname, clobber=True)

   def merge_runs(self, band, refcat='run0.cat'):
      curdir = os.getcwd()
      os.chdir('%s_%s' % (self.root, band))
      cr = sextractor(refcat)
      columns = cr._colnames
      allcats = glob.glob('run*.cat')
      for ac in allcats:
         if ac == refcat:
            continue
         else:
            ca = sextractor(ac)
            for col in columns:
               assert hasattr(ca, col), "Column %s does not exist in %s..." % (col, ac)
               newcol = np.concatenate([getattr(cr,col), getattr(ca,col)])
               setattr(cr, col, newcol) 
      cr.writeto('%s_%s_output.cat' % (self.root, band), clobber=True)
      os.chdir(curdir)

datum = time.strftime("%m%d",time.localtime())

def merge_sim_catalogs(catalogs, output):
   """
   Merge simulated catalogs from each run. Assume that they all have the same
   columns!!
   """
   with open(output, 'wb') as f:
      head = 0
      for cat in catalogs:
         with open(cat, 'rb') as c:
            lines = c.readlines()
         if head == 0:
            for l in lines:
               f.write(l)
            head = 1
         else:
            for l in lines:
               if l[0] != '#':
                  f.write(l)
   print "Done."


def run_skyrms_sim(parfile, plot=False):
   sim = SkyRMSSim(parfile)
   while sim.n < sim.nstop:
      print "Iteration %d:" % sim.n
      sim.insert_fake_sources()
      sim.run_sextractor()
      if sim.fit_sky:
         # Now estimate local sky at the location of each source
         for b in sim.bands:
            sim.fit_sky(b, sim.n, skyaper=sim.skyaper, plot=plot)
      sim.cleanup()
      sim.n = sim.n + 1

def collect_fluxerr(band, root='skyrms', apercol='aper_1', fit_sky=False, ra_center=-99, dec_center=-99, radius=20):
   """
   Collect the SExtractor flux error and Gaussian-fitted flux error, within 
   the same aperture
   Optionally, the function can return just the fluxes around a given
   coordinate within some given radius (in arcsec). By default, return fluxes 
   across the entire image.
   """
   allruns = glob.glob('%s_%s/run*.cat' % (root, band))
   SEx_flux_aper = np.zeros(0)
   SEx_fluxerr_aper = np.zeros(0)
   sky_fluxerr_aper = np.zeros(0)
   if os.path.exists('%s_%s/%s_%s_output.cat' % (root,band,root,band)):
      catalog = '%s_%s/%s_%s_output.cat' % (root,band,root,band)
      print "Reading %s..." % catalog
      c = sextractor(catalog)
      within = np.ones(len(c), 'bool')
      if (ra_center>-99) and (dec_center>-99):
         angdist = angsep.angsep(c.alpha_j2000, c.delta_j2000, 
                                 ra_center, dec_center)
         within = ((angdist * 3600.) <= radius)
      SEx_flux_aper = getattr(c, 'flux_%s' % apercol)
      SEx_fluxerr_aper = getattr(c, 'fluxerr_%s' % apercol)
      SEx_flux_aper = SEx_flux_aper[(c.imaflags_iso==0)&within]
      SEx_fluxerr_aper = SEx_fluxerr_aper[(c.imaflags_iso==0)&within]
   else:
      for r in allruns:
         c = sextractor(r)
         if (ra_center>-99) and (dec_center>-99):
            angdist = angsep.angsep(c.alpha_j2000, c.delta_j2000, 
                                    ra_center, dec_center)
            within = ((angdist * 3600.) <= radius)
         else:
            within = np.ones(len(c), 'bool')
         flux_aper_col = getattr(c, 'flux_%s' % apercol)[(c.imaflags_iso==0)&within]
         fluxerr_aper_col = getattr(c, 'fluxerr_%s' % apercol)[(c.imaflags_iso==0)&within]
         within = np.ones(len(c), 'bool')
         SEx_flux_aper = np.concatenate([SEx_flux_aper, flux_aper_col])
         SEx_fluxerr_aper = np.concatenate([SEx_fluxerr_aper, 
                                           fluxerr_aper_col])
         if fit_sky:
            sky_fluxerr_aper = np.concatenate([sky_fluxerr_aper, 
                               c.sky_sigma_aper[(c.imaflags_iso==0)&within]])
   if fit_sky:
      return SEx_flux_aper, SEx_fluxerr_aper, sky_fluxerr_aper
   else:
      return SEx_flux_aper, SEx_fluxerr_aper

def calc_local_fluxerr(band, ra, dec, magzero, radius=20., root='skyrms', apercol='aper_1'):
   """
   Calculates the 1-sigma flux measured at empty regions within the aperture 
   provided by the user. The 1-sigma flux is derived from simulated sources
   centered at (ra, dec) within radius (in arcsec).
   """
   flux, fluxerr = collect_fluxerr(band, root=root, apercol=apercol, 
                                   ra_center=ra, dec_center=dec, radius=radius)
   flux_stats = gauss.fitgauss(flux[fluxerr>0], clip=True, disp=0)
   print "A total of %d sources used in %s for statistics." % (np.sum(fluxerr>0), band.upper())
   print "Mean & sigma of the flux around (RA, DEC) = (%.7f, %.7f) are %.3e, %.3e" % (ra, dec, flux_stats[0], flux_stats[1])
   print "This corresponds to a 1-sigma magnitude limit of %.3f." % (magzero - 2.5 * np.log10(flux_stats[1]))

def calc_local_scaling_factor(band, ra, dec, magzero, radius=20., root='skyrms', apercol='iso', fit_sky=False):
   sex_flux, sex_fluxerr = collect_fluxerr(band, root=root, apercol=apercol,
                                           ra_center=ra, dec_center=dec, 
                                           radius=radius, 
                                           fit_sky=fit_sky)
   assert len(sex_flux) > 0, "No simulated sources found..."
   print "Number of simulated sources: ", len(sex_flux)
   flux_stats = gauss.fitgauss(sex_flux[(sex_fluxerr>0)], clip=True, disp=0)
   fluxerr_stats = gauss.fitgauss(sex_fluxerr[(sex_fluxerr>0)], clip=True, disp=0)
   factor = flux_stats[1] / fluxerr_stats[0]
   print "Scaling factor in %s from median(fluxerr) to sigma(flux_empty) is %.3f" % (band.upper(), factor)
   return factor

def calc_scaling_factor(band, root='skyrms', apersize=0.4, apercol='iso', fit_sky=False):
   # Use ISO aperture to calculate the scaling factor, because we know the 
   # number of pixels (isoarea_image) SExtractor uses to calculate flux errors
   sex_flux, sex_fluxerr = collect_fluxerr(band, root=root, apercol=apercol,
                                           fit_sky=fit_sky)
   # fig = plt.figure()
   # ax = fig.add_subplot(111)

   ### This is how I calculated the scaling factor, which is comparing SExtractor
   ### flux error directly with sqrt(A * sky_sigma**2), where A is the total
   ### area of the aperture. However this is problematic b/c SExtractor also
   ### includes shot noise (from object flux) into their flux error estimate,
   ### so this comparison misses one component --- the object flux.
   # factor = sky_fluxerr[sex_fluxerr>0] / sex_fluxerr[sex_fluxerr>0]
   # print factor.min(), factor.max()
   # madn = robust.MADN(factor)
   # med_factor = np.median(factor)
   # xmin = med_factor - 10. * madn
   # xmax = med_factor + 10. * madn
   # h = ax.hist(factor, bins=np.linspace(xmin, xmax, 40), histtype='step', 
   #             lw=2.0)
   # print "Total number of sources: %d" % len(sex_fluxerr)
   # print "Median scaling factor is %.3f" % med_factor
   # ax.plot([med_factor] * 2, [0, h[0].max()], ls='--', lw=3.0, 
   #          c='black')
   # ax.set_title('Scaling factor in %s within %.1f" diameter aperture' % (band, apersize), size=22)
   # ax.set_xlabel('sky_fluxerr / SEx_fluxerr', size=16)
   # ax.text(0.95, 0.95, 'N=%d\nmedian=%.3f' % (len(sex_fluxerr), med_factor), 
   #         ha='right', va='top', transform=ax.transAxes, size=20,
   #         bbox=dict(boxstyle='round',facecolor='wheat',alpha=0.5))
   # ax.set_xlim(h[1].min(), h[1].max())
   # fig.savefig('%s_sim_%s.png' % (root, band))
   ### Now I'm just gonna follow the recipe in Trenti et al. 2011, where I 
   ### compare median(SEx_fluxerr) and sigma(SEx_flux), where SEx_flux are the
   ### fluxes at empty regions. Sigma(SEx_flux) is what flux errors should match.
   flux_stats = gauss.fitgauss(sex_flux[(sex_fluxerr>0)], clip=True, disp=0)
   fluxerr_stats = gauss.fitgauss(sex_fluxerr[(sex_fluxerr>0)], clip=True, disp=0)
   factor = flux_stats[1] / fluxerr_stats[0]
   print "Scaling factor in %s from median(fluxerr) to sigma(flux_empty) is %.3f" % (band.upper(), factor)


if __name__ == '__main__':
   # Read the configuration file and construct arguments for 
   #  simulate() and measure_petrosian()
 
   t1 = time.time()
   curdir = os.getcwd()
   
   parfile = sys.argv[1]
   run_sexsim(parfile)
   
   print "Finished simulation"