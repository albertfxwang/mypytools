#!/usr/bin/env python
# This script came from newgalsimfit.py, but with several important changes:
# - does not split the mosaic into subpanels (for better tracking of object
#   locations)
# - gets rid of any perl remnants (no part of the simulation uses perl anymore)
# 
# K-H Huang, 2013/11/22

from pyraf import iraf
from iraf import images, stsdas
import os, glob, subprocess, sys
import numpy as np
import pyfits
from galfitsim_util import galfit_measure

def galsimfit(image, sigimg, sexcat_target, sexcat_all, startnum, stopnum, 
              psffile, constraint,
              band, magzpt=22., plate=0.03, nofat=True, mask=True, 
              maskimg='mask.fits', gsigma=10.):
    # image: Input image to be fitted
    # sigimg: Sigma image to use
    # sexcat: Sextractor catalog name
    # startnum: Start from which object number in the catalog?
    # stopnum: Stop at which number in the catalog?
    # psfloc: PSF location (directory)
    # constraint: Parameter constraint file
    # fitfunc: Single sersic or B/D decomposition?
    # region: (A) Whole image catalog or, (B) region of image?
    # boundary: Boundary of region if (B) (xmin,xmax,ymin,ymax)
    # skyval: Sky background to add back to Sextractor value
    # mrange: Acceptable magnitude range (min, max)
    # stargal: Range of good galaxies: 0.0 (gal) to 1.0 (star)
    # fwhmpsf: FWHM of the PSF for cosmic ray rejection
    # magzpt: Photometric zeropoint
    # plate: Image plate scale
    # exist: Rid all existing temp-*.fits images before run?
    # nofat: Rid all temp images, galfit input, after run?
    # nsplit: Split image into how many parts along one axis?
    # mask: Use mask image in galfitting?
    # maskimg: Mask image name
    print "Mask image is %s" % maskimg
    # BUFFER = 200

    # Write input parameter file for galfit_measure.py
    f = open('galfitsim.pars.yml', 'wb')
    f.write('startnum: %d\n' % startnum)
    f.write('stopnum: %d\n' % stopnum)
    f.write('drzimage: %s\n' % image)
    f.write('rmsimage: %s\n' % sigimg)
    f.write('segimage: %s\n' % maskimg)
    f.write('psfimage: %s\n' % psffile)
    f.write('magzpt: %f\n' % magzpt)
    f.write('consfile: %s\n' % constraint)
    f.write('pixscale: %.2f\n' % plate)
    f.write('galfitpath: /Users/khuang/bin \n')  # hard-coded for now; should change later
    f.write('target_catalog: %s\n' % sexcat_target)
    f.write('all_catalog: %s\n' % sexcat_all)
    f.write('band: %s\n' % band)
    f.write('id_column: number\n')
    f.write('ra_column: alpha_j2000\n')
    f.write('dec_column: delta_j2000\n')
    f.write('re_column: flux_radius_1\n')
    f.write('mask: %s\n' % mask)
    f.write('mag_column: mag_auto\n')
    f.write('x_column: x_image\n')
    f.write('y_column: y_image\n')
    f.write('ell_column: ellipticity\n')
    f.write('theta_column: theta_image\n')
    f.write('isoa_column: isoarea_image\n')
    f.write('overwrite: True\n')
    f.write('nblocks_x: 1\n')
    f.write('nblocks_y: 1\n')
    f.write('padding: 100\n')
    f.write('split_image: no\n')
    f.close()

    # Run GALFIT
    galfit_measure.run_galfit('galfitsim.pars.yml')

    #  Delete all temporary images if the user wants this

    if nofat:
        os.system('rm obj[0-9]*')
        os.system('rm galfit.[0-9]*')
    #try:
    #    os.rmdir("galfit-tempfits")
    #except OSError:  
    #    pass
    
