#!/usr/bin/env python
# This is a wrapper script for galfit.pl.  ACS images are huge, and it takes
# more time for galfit to read in an image than it takes to fit a galaxy.
# Thus, this script takes an ACS image and splits it up into little chunks.
# Then, it outputs a temporary "perl.in" file which is read in by the
# "galfit.pl" Perl script.  The Perl script reads in a Sextractor file
# for galaxy information, formats it into galfit input file, figures out
# which galaxy falls into which panel, then finally runs galfit.
#
# Once all the galaxies have been fit, we return to this ***PYTHON*** script
# This script then reformats the headers of all the postage stamp images to
# reflect the correct information relative to the big ACS image.
#
# Parameter explanations:
#   skyval - The value of the sky background in the image being fitted that
#            SExtractor knows nothing about.  This value gets added to 
#            background determined by SExtractor when creating the galfit 
#            input file, and doesn't actually change the image being fitted
#            or the catalog itself.  This is needed because sometimes people
#            run Sextractor on images that have the background removed, only
#            to realize they want GALFIT to create a sigma map internally,
#            which needs the sky background intact.  If this is the case,
#            the user has to add the sky background back into the image
#            manually.

#############################################################################
#   Re-written from newgalsimfit.cl by K-H Huang (12/17/08)                 #
#############################################################################

from pyraf import iraf
from iraf import images, stsdas
import os, glob, subprocess, sys
from numpy import *
import pyfits

def newgalsimfit(image, sigimg, sexcat, startnum, stopnum, psfloc, constraint,
    fitfunc='sersic',
    region='A', boundary="(0,0,0,0)", skyval=0., mrange="(0.,99.)", stargal="(0.,1.)",
    fwhmpsf=0.06, magzpt=22., plate=0.03, exist=False, nofat=True, nsplit=10,
    mask=False, maskimg='mask.fits',gsigma=10.):
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
    BUFFER = 200

    img = image
    sig = sigimg
    sex = sexcat
    snum = startnum
    stnum = stopnum

    if not os.path.exists("galfit-tempfits"):
        os.mkdir("galfit-tempfits")

    print ""
    print "Putting all sub-panel images into the directory: galfit-tempfits."
    print ""

################  Split image into nsplit panels along a side  ################
    imhead = pyfits.getheader(img)
    # Read the header of the image

    ncol = imhead['naxis1']
    nrow = imhead['naxis2']

    nsubx = int(ncol / nsplit)
    nsuby = int(nrow / nsplit)

    # Get rid of the last "/" in the directory name:

    dum = psfloc[-1]
    if dum == "/":
        psfloc = psfloc[:-1]

    # Check to see if we should first delete existing sub-panels

    if exist:
        os.system("rm galfit-tempfits/temp-*.fits")
        os.system("rm galfit-tempfits/sig-*.fits")
        os.system("rm mask-*.fits")

    for j in range(1, nsplit+1):
        for i in range(1, nsplit+1):
            if (nsubx * (i-1) - BUFFER <= 1.):
                xlo = 1
            else:
                xlo = nsubx * (i-1) - BUFFER + 1
            if (nsuby * (j-1) - BUFFER <= 1.):
                ylo = 1
            else:
                ylo = nsuby * (j-1) - BUFFER + 1
            if (nsubx * i + BUFFER >= ncol):
                xhi = ncol
            else:
                xhi = nsubx * i + BUFFER
            if (nsuby * j + BUFFER >= nrow):
                yhi = nrow
            else:
                yhi = nsuby * j + BUFFER

            tempout = "galfit-tempfits/temp-%d-%d.fits" % (i,j)
            sigout = "galfit-tempfits/sig-%d-%d.fits" % (i,j)
            maskout = "mask-%d-%d.fits" % (i,j)
            if (os.path.exists(img)==True) & (os.path.exists(tempout)==False):
                iraf.imcopy(input=img+"[%d:%d,%d:%d]" % (xlo,xhi,ylo,yhi),
                            output=tempout)
            if (os.path.exists(sig)==True) & (os.path.exists(sigout)==False):
                iraf.imcopy(input=sig+"[%d:%d,%d:%d]" % (xlo,xhi,ylo,yhi),
                            output=sigout)
            if (os.path.exists(maskimg)==True) & (os.path.exists(maskout)==False):
                iraf.imcopy(input=maskimg+"[%d:%d,%d:%d]" % (xlo,xhi,ylo,yhi),
                            output=maskout)

#######################  Set up and run Perl script  ##########################

    f = open("perl.in","w")

    print >>f, img, "   # The ACS image name to fit"
    print >>f, sex, "   # The sextractor catalog"
    print >>f, psfloc, "   # PSF directory"
    print >>f, snum, "   # Starting from object number"
    print >>f, stnum, "  # Stop at object number"
    print >>f, constraint, "   # Parameter constraint file"
    print >>f, "galfit-tempfits   # Temporary working directory"
    print >>f, fitfunc, "   # Fit a single sersic or do B/D decomposition?"
    print >>f, region, "    # (A) Whole image catalog or, (B) region of the image"
    print >>f, boundary, "   # Boundary of the object region if (B)"
    print >>f, skyval, "    # Sky value to add to Sextractor determined value"
    print >>f, mrange, "    # Acceptable magnitude range"
    print >>f, stargal, "    # Range of good galaxies: 0.0 (galaxy) to 1.0 (star)"
    print >>f, fwhmpsf, "    # FWHM of the PSF for cosmic ray rejection"
    print >>f, magzpt, "     # Photometric zeropoint"
    print >>f, plate, "     # Plate scale"
    print >>f, ncol, nrow, "    # ncol, nrow of original ACS image"
    print >>f, "temp    # Temporary sub-panel file names"
    print >>f, nsplit, "    # Number of subpanels along a side to split ACS image"
    print >>f, BUFFER, "    # Annulus width to pad around the split image"
    print >>f, mask, "   # To use mask image or not"
    print >>f, gsigma, "    # sigma of Gaussian wing in masks"

    f.close()

    pf = open("psf.temp","w")
    psflist = glob.glob(psfloc+"/psf*.fits")
    for psfimg in psflist:
        print >>pf, psfimg
    pf.close()

######### Enter into Perl(or python??) script, which runs GALFIT

    galfit = subprocess.Popen(['galfitpl.py','perl.in'],stdout=sys.stdout) 
    retcode2 = galfit.wait()
    print retcode2

######### After GALFIT has been run on ALL the galaxies in the database 
######### modify the header info of the image blocks.

    fileptr = open("galfit-tempfits/offsets","r")
    lines = fileptr.readlines()
    for l in lines:
        if l[0] == '#': continue
        objn,stampname,offx,offy,section,xcboxcent,ycboxcent,ncomp = l.split()
        if os.path.exists(stampname+".fits"):
            h = pyfits.open(stampname+".fits","update")
            if len(h) != 4: pass
            else:
               h[2].header.update("objnum",objn)
               h[1].header.update("object",(img+section))
               h[2].header.update("datain",img)
               h[2].header.update("noise",sig)
               h[2].header.update("fitsect",section)
               h[2].header.update("cboxcent","%s, %s" % (xcboxcent,ycboxcent))
               h[2].header.update("ncomp",ncomp)
               xmin, xmax = section[1:-1].split(',')[0].split(':')
               xmin = int(xmin); xmax = int(xmax)
               h[2].header.update("xmin", xmin)
               h[2].header.update("xmax", xmax)
               ymin, ymax = section[1:-1].split(',')[1].split(':')
               ymin = int(ymin); ymax = int(ymax)
               h[2].header.update("ymin", ymin)
               h[2].header.update("ymax", ymax)
            try:
               h.flush(output_verify='fix')
            except:
               h.flush(output_verify='ignore')

#  Delete all temporary images if the user wants this

    if nofat:
        #os.system("rm galfit-tempfits/temp-*.fits")
        #os.system("rm galfit-tempfits/sig-*.fits")
        #os.system("rm galfit-tempfits/offsets")
        #os.system('rm obj[0-9]*')
        os.system('rm galfit.[0-9]*')
        #os.system("rm perl.in")
        os.system("rm psf.temp")
    #try:
    #    os.rmdir("galfit-tempfits")
    #except OSError:  
    #    pass
    
