#!/usr/bin/env python

import numpy as np
from sexsim import sextractor_sim
import os, glob
import pyfits, pywcs
from galfit_measure import merge_gf_results as mgr
from galsimfit import galsimfit
from pygoods import sextractor
import segmap_30mas
from pyraf import iraf
iraf.artdata()
from hconvolve import imhconvolve

"""
Run GALFIT measurements after running SExtractor on images with fake sources.
"""

gsigma = 10.   # sigma of Gaussian wing in pixels

class galfit_sim(sextractor_sim.sextractor_sim):
   def __init__(self, parfile):
      sextractor_sim.sextractor_sim.__init__(self, parfile)
      self.magfile = None
      if self.c.has_key('GALFIT_BANDS'):
         self.galfit_bands = self.c['GALFIT_BANDS']
      else:
         self.galfit_bands = self.c['BANDS']
      print "GALFIT Bands are: ", self.galfit_bands
      ## Add capability to handle different pixel sizes
      self.mask = self.c['MASK']
      self.pixscale_galfit = self.c['SCALE_GALFIT']
      self._multires_galfit = False
      if self.pixscale_galfit != self.pixscale:
         self._multires_galfit = True
         self.realimages_galfit = dict(zip(self.bands, self.c['REALIMAGE_GALFIT']))
         self.rmsimages_galfit = dict(zip(self.bands, self.c['RMSIMAGE_GALFIT']))
         self.flagimages_galfit = dict(zip(self.bands, self.c['FLAGIMAGE_GALFIT']))
         self.psffiles_galfit = dict(zip(self.bands, self.c['PSFFILE_GALFIT']))
         self.fakeimages_galfit = {}
         hdr = pyfits.getheader(self.realimages_galfit[self.galfit_bands[0]])
         self.xmax_galfit = hdr['naxis1']
         self.ymax_galfit = hdr['naxis2']
         assert self.xmax_galfit != self.xmax, "Image dimension of new GALFIT image should be different!"
         self.zeropoints_galfit = dict(zip(self.bands, self.c['MAGZPT_GALFIT']))

   def makenoiselessimage_galfit(self, band, galfile, magz, psffile, save=0, 
                                 gain=1.0):
      """
      Creates a noiseless convolved image
      """
      assert os.path.exists(psffile), "PSF image %s does not exist." % psffile
      broot = self.root + '_%s' % band
      outfile = broot + '_sim.fits'
      outfile_nonoise = broot + '_sim_noiseless.fits'
      if os.path.exists(outfile):
         os.remove(outfile)
      # print outfile,xmax,ymax,galfile,magz,gain
      iraf.unlearn('artdata')
      iraf.unlearn('mkobjects')
      iraf.artdata.dynrange=1.e5
      print "Running iraf.mkobject in %s for GALFIT..." % band
      iraf.mkobjects(outfile_nonoise, output="", title="", ncols=self.xmax_galfit, 
         nlines=self.ymax_galfit, header="", background=0.0, objects=galfile,
         xoffset=0., yoffset=0., star="gaussian", radius=0.1,
         beta=2.5, ar=1., pa=0., distance=1., exptime=1., 
         magzero=magz, gain=gain, rdnoise=0., poisson=0,
         seed=2, comments=1)
      imhconvolve(outfile_nonoise, psffile, outfile, overwrite=True)
      self.noiselessimages[band] = outfile
      os.remove(outfile_nonoise)

   def addsimulated_galfit(self, band, save=0):
      """
      Add the noiseless images of artificial galaxies to the real images. 
      Specifically for GALFIT measurement images with different pixel scales 
      as the detection image.
      """
      # simulation = root+'_sim.fits'
      assert os.path.exists(self.noiselessimages[band]), \
         "Noiseless image with artificial galaxies not calculated."
      broot = self.root + '_%s' % band
      outimage = broot + '_galfit.fits'
      if os.path.exists(outimage):
         os.remove(outimage)
      noiseless_img = pyfits.getdata(self.noiselessimages[band])
      realimage_img = pyfits.getdata(self.realimages_galfit[band])
      hdr = pyfits.getheader(self.realimages_galfit[band])
      simulated_img = realimage_img + noiseless_img
      pyfits.append(outimage, simulated_img, hdr)
      self.fakeimages_galfit[band] = outimage
      if not save:
         os.remove(self.noiselessimages[band])

   def insert_fake_sources_galfit(self):
      """
      Re-insert fake sources into GALFIT measurement images, because GALFIT
      measurement images have different pixel scales than the detection images.
      """
      assert hasattr(self, 'fg'), "Generate attributes for fake sources first."
      assert self.c.has_key('REALIMAGE_GALFIT'), "REALIMAGE_GALFIT not provided."
      assert self.c.has_key('RMSIMAGE_GALFIT'), "RMSIMAGE_GALFIT not provided."
      assert self.c.has_key('FLAGIMAGE_GALFIT'), "FLAGIMAGE_GALFIT not provided."
      assert self.c.has_key('PSFFILE_GALFIT'), "PSFFILE_GALFIT not provided."
      # If pixel scales are different, one needs to recalculate the (x, y) 
      # positions of fake sources (assuming that WCS information is consistent),
      # and one needs to re-calculate the input galaxy sizes.
      hdr_detect = pyfits.getheader(self.realimages[self.detect_band])
      hdr_measure = pyfits.getheader(self.realimages_galfit[self.galfit_bands[0]])
      wcs_detect = pywcs.WCS(hdr_detect)
      wcs_measure = pywcs.WCS(hdr_measure)
      radec = wcs_detect.wcs_pix2sky(np.array([self.fg.x, self.fg.y]).T, 1)
      xy_measure = wcs_measure.wcs_sky2pix(radec, 1)
      self.fg.x_detect = self.fg.x.copy()
      self.fg.y_detect = self.fg.y.copy()
      self.fg.x = xy_measure[:,0]
      self.fg.y = xy_measure[:,1]
      self.fg.re_detect = self.fg.re.copy()
      # scale the input Re to match the new pixel scale
      self.fg.re = self.fg.re_detect * (self.pixscale / self.pixscale_galfit)
      # Now insert galaxies
      # makegals_multiband is inherited from class fake_galaxies
      self.fg.makegals_multiband(self.flagimages_galfit[self.galfit_bands[0]], 
                                 bands=self.galfit_bands)
      for b in self.galfit_bands:
         self.makenoiselessimage_galfit(b, self.fg.artfiles[b], 
                                 self.zeropoints_galfit[b],
                                 self.psffiles_galfit[b])
         self.addsimulated_galfit(b)

   def resample_segmap(self):
      """
      Resample the seg-map from SExtractor run, to make a segmap in the same
      pixel resolution as the GALFIT measurement image.
      """
      assert self._finished_sextractor, "Has not run SExtractor yet!"
      print "Resampling segmentation map..."
      newseg = segmap_30mas.segmap_to_30mas(self.segimages[self.detect_band])
      for b in self.galfit_bands:
         self.segimages[b] = newseg

   def make_maskimg(self, band):
      broot = self.root + '_%s' % band
      self.maskimg = "%s_masknew.fits" % broot
      if os.path.exists(self.maskimg):
         os.remove(self.maskimg)
      segnew = pyfits.getdata(self.segimages[band])
      print segnew.shape
      maxseg = segnew.ravel().max()
      os.system('cp %s %s' % (self.segimages[band], self.maskimg))
      if self._multires_galfit:
         flagmap = pyfits.getdata(self.flagimages_galfit[band])
      else:
         flagmap = pyfits.getdata(self.flagimages[band])
      h = pyfits.open(self.maskimg, mode='update')
      masknew = np.where(flagmap.ravel() > 0, 1000000, segnew.ravel())
      masknew = masknew.reshape(segnew.shape)
      h[0].data = masknew
      h.flush()
      h.close()

   def run_galfit(self):
      assert self._finished_sextractor, "Has not run SExtractor yet!"
      if len(glob.glob('mask*.fits')):
         os.system('rm mask*.fits')
      for b in self.bands:
         if b in self.galfit_bands:
            # will delegate the work of Gaussian smoothing 
            # to galfit_measure.py
            # Now incorporate the flag map into the mask image
            # Set flagged pixel to 1 + maximum segmentation number
            # so the flagged pixels won't be zero-ed out in galfitpl.py 
            broot = self.root + '_%s' % b
            if self.mask:
               self.make_maskimg(b)
            gfile = "glart_%s.list" % (b)  # artdata output
            try:
               c = sextractor(broot+'/run%d.cat' % self.n)
            except:
               print "No galaxies detected in this iteration."
               break
            # os.system("mv %s %s/run%d.cat" % (self.catalogs[b], broot, self.n))
            # os.system("mv %s %s/run%d.newcat" % (self.newcatalogs[b], broot, self.n)) 
            # os.system("mv %s %s/glart%d.list" % (gfile, broot, self.n))                
            # os.system("cp %s %s" % (self.psffiles[b], broot))
            #---------------------------------------=
            # Run galfit using single sersic function
            #----------------------------------------
            os.chdir(broot)
            print "%s" % broot
            gfitfile = "gfit%d.cat" % (self.n)
            # nsexfile = "run%d.joincat" % (self.n)
            # gnsortfile = "run%d.nosortcat" % (self.n)
            c1 = sextractor('run%d.cat' % (self.n))
            cnew = sextractor('run%d.newcat' % (self.n))
            f1 = open(gfitfile, "w")
            f1.write('# 1 GALFIT_NUMBER\n')
            f1.write('# 2 FLUX_AUTO\n')
            f1.write('# 3 FLUXERR_AUTO\n')
            f1.write('# 4 MAG_AUTO\n')
            f1.write('# 5 MAGERR_AUTO\n')
            f1.write('# 6 FLUX_RADIUS_1\n')
            f1.write('# 7 BACKGROUND\n')
            f1.write('# 8 ISOAREA_IMAGE\n')
            f1.write('# 9 X_IMAGE\n')
            f1.write('# 10 Y_IMAGE\n')
            f1.write('# 11 ALPHA_J2000\n')
            f1.write('# 12 DELTA_J2000\n')
            f1.write('# 13 THETA_IMAGE\n')
            f1.write('# 14 ELLIPTICITY\n')
            f1.write('# 15 FWHM\n')
            f1.write('# 16 IMAFLAGS_ISO\n')
            f1.write('# 17 CLASS_STAR\n')
            f1.write('# 18 NUMBER\n')
            fluxauto = c1.flux_auto
            fluxautoerr = c1.fluxerr_auto
            magauto = c1.mag_auto
            magautoerr = c1.magerr_auto
            halflight_rad = c1.flux_radius_1
            bkg = c1.background
            isoarea = c1.isoarea_image
            ximage = c1.x_image
            yimage = c1.y_image
            alpha  = c1.alpha_j2000
            delta  = c1.delta_j2000
            theta  = c1.theta_image
            ellip = c1.ellipticity
            fwhm = c1.fwhm_image
            flag = c1.imaflags_iso
            class_star = c1.class_star
            # Find the SEx ID for each fake galaxy in *.newcat
            objid = np.zeros(len(c1), 'int')
            for i in range(len(objid)):
               dist = np.sqrt((cnew.x_image-c1.x_image[i])**2+(cnew.y_image-c1.y_image[i])**2)
               j = np.argsort(dist)[0]
               objid[i] = cnew.number[j]
            # for i in range(len(objid)):
               line =  "%4d %6.4f %6.4f %8.3f %8.3f %5.3f %8.5f %5.2f %8.3f %8.3f %8.5f %8.5f %5.2f %5.3f %5.3f %2d %5.3f %d\n" % (i,
                     fluxauto[i],fluxautoerr[i],magauto[i],
                     magautoerr[i],halflight_rad[i],bkg[i],isoarea[i],
                     ximage[i],yimage[i],alpha[i],delta[i],theta[i],
                     ellip[i],fwhm[i],flag[i],class_star[i],objid[i])
               f1.write(line)
            f1.close()
            os.chdir('../')
            catname=broot+"/"+gfitfile
            newcatname = broot + '/run%d.newcat' % self.n
            consname3='constraints3'
            outname1=broot+"/gfit%d_s.cat" % (self.n)
            sc = sextractor(broot+'/run%d.cat' % self.n)
            nsex = len(sc)
            os.system('rm galfit_output/obj*.fits')
            os.system('rm galfit_input/obj*')
            if self._multires_galfit:
               galsimfit(image=self.fakeimages_galfit[b], 
                         sigimg=self.rmsimages_galfit[b],
                         sexcat_target=catname, sexcat_all=newcatname,
                         band=b, startnum=0, stopnum=len(sc),
                         psffile=self.psffiles_galfit[b], constraint=consname3,
                         magzpt=self.zeropoints_galfit[b], 
                         plate=self.pixscale_galfit, nofat=True, 
                         mask=self.mask, maskimg=self.maskimg, 
                         gsigma=gsigma)
            else:
               galsimfit(image=self.fakeimages[b], sigimg=self.rmsimages[b], 
                      sexcat_target=catname, sexcat_all=newcatname, 
                      band=b, startnum=0, stopnum=len(sc), 
                      psffile=self.psffiles[b], constraint=consname3, 
                      magzpt=self.zeropoints[b], plate=self.pixscale, 
                      nofat=True, mask=self.mask, maskimg=self.maskimg, 
                      gsigma=gsigma)
            # Merge results from GALFIT output image headers
            # readgfheader.write_gf_output('galfit_output', 'obj*_out.fits', outname1)
            mgr.merge_galfit(outname1, 'galfit_output')
            mgr.add_crashes(outname1, 'galfit_input')
            mgr.add_objid(outname1)

   def match_galfit_catalog(self):
      for b in self.galfit_bands:
         broot = self.root + '_%s' % b
         os.chdir(broot)
         c1 = sextractor('run%d.cat' % self.n)  # SExtractor catalog
         c2 = sextractor('gfit%d_s.cat' % self.n)  # GALFIT catalog
         ncolumns1 = len(c1._colnames)
         ncolumns2 = len(c2._colnames)
         ncolumns_tot = ncolumns1 + ncolumns2
         header = ""
         for i in range(ncolumns1):
            header += "#  %d %s\n" % ((i+1), c1._colnames[i].upper())
         for i in range(ncolumns2):
            header += "#  %d %s\n" % ((i+ncolumns1+1), c2._colnames[i].upper())
         f = open('gfit%d_s_matched.cat' % self.n, 'wb')
         f.write(header)
         for i in range(len(c1)):
            j = np.arange(len(c2))[c2.objid==c1.number[i]][0]
            f.write(' '.join(c1._colentries[i]))
            f.write(' ')
            f.write(' '.join(c2._colentries[j]))
            f.write(' ')
            f.write('\n')
         print "Made matched catalog in %s for iteration %d." % (b, self.n)
         f.close()
         os.chdir('../')

   def cleanup(self):
      if not self.save:
         os.system('rm galfit_output/obj*.fits')
         print "removing %s*.fits" % self.root
         os.system('rm %s*.fits' % self.root)
         os.system('rm subpanels/*.fits')

