#!/usr/bin/env python

import os, glob, sys, subprocess, time, signal
from pyraf import iraf
from iraf import stsdas, images
from numpy import *
import numpy as np
import pyfits
from scipy import ndimage  # to use its Gaussian filter capability

PI = 3.1415926535
MAGDIFF = 1

#######################################################
#    Re-written from galfit.pl to incorporate         #
#    selective masking by K-H Huang 12/17/08          #
#######################################################

#################### Example of galfit.pl control parameters  #################

# $ACSIMG = "s1z09.fits"  # The ACS image name to fit.
# $SEXCAT = "s1z09.cat"   # The Sextractor image catalogue name.
# $PSFLOC = "psfs"        # The directory storing all the PSFs.
# $STARTNUM = 1       # Object starting number to fit.
# $STOPNUM = 999       # Object starting number to fit.
# $CONSFILE = "constraints"  # The constraint file name.
# $TEMPLOC = "galfit-tempfits"  # The directory storing image sub-tiles.
# $FITFUNC = "sersic"     # Fit a single sersic or do B/D decomposition?
# $REGION = "A"           # (A) Whole image catalog or, (B) region of the image;
# $BOUNDARY = "  "        # Boundary of the object region if (B);
# $SKYVAL = 0.           # Sky value added back to an image that SExtractor knows nothing about.
# $MRANGE = "0. 99."      # Acceptable magnitude range;
# $STARGAL = "0.0 1.0"    # Range of good galaxies: 0.0 (galaxy) to 1.0 (star);
# $FWHMPSF = 0.05         # FWHM of the PSF for cosmic ray rejection.
# $MAGZPT = 23.64         # Photometric zeropoint.
# $PLATE = 0.03           # Plate scale
# $NCOL, $NROW = "512 512"  # The image size of the sub-tiles.
# $TNAME = "temp"      # The prefix name of the image sub-tiles.
# $BUFFER = "200"      # The padding size around the image edges

##############################################################################
##############################################################################

def components(objn,centobj,curobj,xoff,yoff):
    # Calculate the (x,y) position of the current object relative to
    # the tile in which it lives.
    global x,y,FITFUNC,mag,Re,axrat,Out1

    xpos = x[curobj] - xoff
    ypos = y[curobj] - yoff

    if (FITFUNC == 'sersic') or (centobj != curobj):
        n = 1.5
        msersic = mag[curobj]
        Re_sersic = Re[curobj]
        ar_sersic = axrat[curobj]
    elif centobj == curobj:
        n = 2.5
        msersic = mag[curobj] + 0.5
        ar_sersic = 0.8
        Re_sersic = Re[curobj] * 0.7
        mexpdisk = mag[curobj] + 0.5
        Rs = Re[curobj]

    print >>OUT1, "# Object number: %d " % objn
    print >>OUT1, " 0)     sersic         #    Object type "
    print >>OUT1, " 1) %f  %f  1 1    #   position x, y " % (xpos,ypos)
    print >>OUT1, " 3) %f   1       #  total magnitude  " % msersic
    print >>OUT1, " 4) %f 1       #      R_e " % Re_sersic
    if (FITFUNC == 'sersic') and (centobj == curobj):
        print >>OUT1, " 5) 1.5   1  #  exponent (de Vaucouleurs = 4) "
    else:
        print >>OUT1, " 5) %f    1  #  exponent (de Vaucouleurs = 4) " % n
    print >>OUT1, " 9) %f 1       #  axis ratio (b/a) " % ar_sersic
    print >>OUT1, "10) %f 1       #  position angle (PA) " % theta[curobj]
    print >>OUT1, " Z) 0                  #  Output option (0 = residual, 1 = Don't subtract) "
    print >>OUT1, ""

    if curobj == centobj and FITFUNC == "BD":
        objn += 1
        print >>OUT1, "# Object number: %d " % objn 
        print >>OUT1, " 0)     expdisk         #    Object type "
        print >>OUT1, " 1) %f  %f  1 1    #   position x, y " % (xpos,ypos)
        print >>OUT1, " 3) %f  1       #  total magnitude  " % mexpdisk
        print >>OUT1, " 4) %f        1       #      R_e " % Rs
        print >>OUT1, " 9) %f 1       #  axis ratio (b/a) " % axrat[curobj]
        print >>OUT1, "10) %f 1       #  position angle (PA) " % theta[curobj]
        print >>OUT1, " Z) 0                  #  Output option (0 = residual, 1 = Don't subtract) "
    # zero-out the objects included as components in the mask image
    return objn,xpos,ypos

def sky(objn,bkgd):

    print >>OUT1, "# Object number: %d " % objn
    print >>OUT1, " 0)        sky         #    Object type "
    print >>OUT1, " 1) %f     1       #  sky background " % bkgd
    print >>OUT1, " 2) 0.     0       #  sky gradient in x"
    print >>OUT1, " 3) 0.     0       #  sky gradient in y"
    print >>OUT1, " Z) 0                #  Output option (0 = residual, 1 = Don't subtract) "
    print >>OUT1, ""
    print >>OUT1, "================================================================================"

##############################################################################
##############################################################################


#subroutine deg2cel.precise  (version 06/09/02) -- re-written to python 12/17/08
#----------
def deg2celestial(ra,dec):
    #Call w/ ($ra_celestial,$dec_celestial) = deg2celestial($ra,$dec);
    #returns celestial coordinates from input of ra and dec in decimal degree
    #format.  This subroutine returns all significant digits in seconds of
    #RA and Dec; thus simpler version of deg2celestial.sub w/out roundoff.
    # NOTE: Run this subroutine in a loop, only does one transf. at a time.
    # 6/9/02: fixed bug for -1<Dec<0.

    ## RA portion
    hr = ra / 15.
    rah,junk = str(hr).split('.')
    min = 60. * float("0.%s" % junk)
    ram, junk = str(min).split('.')
    ras = 60. * float("0.%s" % junk)

    # always want double digits
    if int(rah) < 10:
        rah = "0" + rah
    if int(ram) < 10:
        ram = "0" + ram
    if int(ras) < 10:
        ras = "0" + str(ras)
    ra_celestial = rah + ":" + ram + ":" + str(ras)   # in the format hh:mm:ss

    ## DEC portion
    degree,junk = str(dec).split('.')
    if float(degree) >= 0.0:
        sign = "p"
        try:
            junk1,deg = degree.split('+')
        except:
            deg = degree
    if float(degree) <= 0.0 and dec < 0.0:
        sign = "m"
        junk1,deg = degree.split('-')

    dmin = 60. * float("0.%s" % junk)
    amin,junk = str(dmin).split('.')
    asec = 60. * float("0.%s" % junk)

    # always want double digits
    if int(deg) < 10:
        deg = "0" + deg
    if int(amin) < 10:
        amin = "0" + amin
    if int(asec) < 10:
        asec = "0" + str(asec)

    # put sign back on
    if sign == "m":
        dec_celestial = "-" + deg + ":" + amin + ":" + str(asec)
    if sign == "p":
        dec_celestial = "+" + deg + ":" + amin + ":" + str(asec)

    return (ra_celestial,dec_celestial)    # strings returned

#############################################################################

def uniqname(ra_celest,dec_celest):
    ra1,ra2,ra3 = ra_celest.split(":")
    ra3 = "%.2f" % float(ra3)     # rounds sec to 0.01
    if float(ra3) < 10.:
        ra3 = "0" + ra3
    ra_name = ra1 + ra2 + ra3
    d1,d2,d3 = dec_celest.split(":")
    d3 = "%.1f" % float(d3)       # rounds asec to 0.1
    if float(d3) < 10.:
        d3 = "0" + d3
    if float(d1) < 0.:
        sign = "m"
    if float(d1) > 0.:
        sign = "p"
    d1 = d1[1:]    # deletes sign
    dec_name = sign + d1 + d2 + d3
    out_name = "GEMS" + ra_name + dec_name
    return out_name

##############################################################################

# Choose the appropriate convolution PSF.  Given a list of PSF names with
# names PSF-x-y.fits, located at position (x, y), figure out for a galaxy
# at position xgal, ygal, which PSF is the closest to this object.
#

def psfdist(xg,yg):
    psfs = []
    PSFFILE = open("psf.temp")
    lines = PSFFILE.readlines()
    for l in lines:
        l = l[:-1]    # chomp
        psfs += [l]

    neardist = 1.e6
    for i in range(len(psfs)):
        dum,xp,yp = psfs[i].split('-')
        yp = yp.split('.fits')[0]
        dist = sqrt((xp-xg)**2 + (yp-yg)**2)
        if dist < neardist:
            neardist = dist
            usepsf = psfs[i]
    
    PSFFILE.close()
    return usepsf


############  New subroutine to create GALFIT input files  ###################

def run_galfit():
    global OUT1,OUT2,OUT3,OUT4
    global totobj
    global objid,flux,fluxerr,mag,magerr,Re,bkg,isoa,x,y,alpha,delta
    global fwhm,axrat,flag,index,theta,sexid

    NRe = 2.

    for i in range(totobj):
        # Iterate through ALL objects
        sys.stdout.flush()
        if (objid[i] >= STARTNUM and objid[i] <= STOPNUM and (mag[i] <= ftcut and mag[i] >= brtcut)):
            # if the object number is in the to-be-fitted range

            sobj = objid[i]  # the ID of object to be fitted
            print "Fitting object ID %d" % sobj
            sys.stdout.flush()
            ncomp = 1
            parfile = "obj%d" % sobj
            OUT1 = open(parfile,"w")

            # Experimental image size
            majx = abs(sqrt(isoa[i]) * 4. * sin(theta[i]/180.*PI))   
            # convert degree into radian
            majy = abs(sqrt(isoa[i]) * 4. * cos(theta[i]/180.*PI))
            minx = axrat[i] * sqrt(isoa[i]) * 4. * abs(cos((theta[i]+
              PI/2.)/180.*PI))
            miny = axrat[i] * sqrt(isoa[i]) * 4. * abs(sin((theta[i]+
              PI/2.)/180.*PI))
            Re_max = 20
            if Re[i] > Re_max:
                Re_max = Re[i]

            xsize = majx
            if minx >= majx:
                xsize = minx
            if xsize <= 30.:
                xsize = 30.
            ysize = majy
            if miny >= majy:
                ysize = miny
            if ysize <= 30.:
                ysize = 30.
            print "Experimental image size: (%d, %d)" % (int(round(xsize)), int(round(ysize)))

            # Figure out which subpanel this current object falls into
            ix = int(x[i] / xpansz) + 1
            iy = int(y[i] / ypansz) + 1
            fitname = TNAME + "-%d-%d.fits" % (ix,iy)
            # rmsname = "rms-%d-%d.fits" % (ix,iy)
            rmsname = "sig-%d-%d.fits" % (ix,iy)
            maskname = "mask-%d-%d.fits" % (ix,iy)

            # Now make sure the object coordinate in the subpanel is correct
            if ix == 1:
                xbuffer = 0
            else:
                xbuffer = BUFFER
            if iy == 1:
                ybuffer = 0
            else:
                ybuffer = BUFFER

            # Figure out which PSF is closest to object center
            # psffile = psfdist(x[i],y[i])    
            psffile = glob.glob(PSFLOC+'/psf*.fits')[0]    
            # settle for this just for now

            # Construct postage stamp name
            #ra_cel,dec_cel = deg2celestial(alpha[i],delta[i])   
            #outname = uniqname(ra_cel,dec_cel)   
            outname = "obj%d" % sobj

            ####################################
            #  Create GALFIT input file header #
            ####################################

            print >>OUT1, "================================================================================"
            print >>OUT1, "# IMAGE PARAMETERS "
            print >>OUT1, "A) %s/%s   # Input Data image (FITS file) " % (TEMPLOC,fitname)
            if MASK == "True":
                print >>OUT1, "B) %s-wmask-out.fits   # Output data image block " % outname
            else:
                print >>OUT1, "B) %s-out.fits   # Output data image block " % outname
            print >>OUT1, 'C) %s/%s   # Noise image name (made from data if blank or "none") ' % (TEMPLOC,rmsname)
            print >>OUT1, "D) %s            # Input PSF image for convolution (FITS file) " % psffile
            print >>OUT1, "E) 1                   # PSF oversampling factor relative to data "
            if MASK == "True":
                print >>OUT1, "F) mask.fits           # Bad pixel mask (FITS image or ASCII coord list) " 
            else:
                print >>OUT1, "F) none                # Bad pixel mask (FITS image or ASCII coord list) "
            print >>OUT1, "G) %s           # Parameter constraint file " % CONSFILE
            print >>OUT1, "J) %f             # Magnitude photometric zeropoint " % MAGZPT
            print >>OUT1, "K) %f %f       # Plate scale (dx dy). " % (PLATE,PLATE)
            print >>OUT1, "O) regular             # Display type (regular, curses, both) "
            print >>OUT1, "P) 0                   # Create ouput only? (1=yes; 0=optimize) "
            print >>OUT1, "S) 0                   # Modify/create objects interactively? "
            print >>OUT1, ""
            print >>OUT1, "# INITIAL FITTING PARAMETERS "
            print >>OUT1, ""
            print >>OUT1, "#   For object type, allowed functions are: sersic, nuker, "
            print >>OUT1, "#                       expdisk, devauc, moffat, gaussian. "
            print >>OUT1, ""
            print >>OUT1, "# Objtype:      Fit?         Parameters "
            print >>OUT1, ""

            # Calculate the (x,y) position of the current object relative to
            # the tile in which it lives.

            offx = (ix-1) * xpansz - xbuffer
            offy = (iy-1) * ypansz - ybuffer
            xfit = x[i] - offx
            yfit = y[i] - offy
            #bkgnd = bkg[i] + SKYVAL
            bkgnd = SKYVAL

            # Calculate fitting box needed to plug into galfit header:
            # Should I impose an upper limit in xsize and ysize?

            xmin = int(xfit - xsize)
            if xmin <= 0:
                xmin = 1
            xmax = int(xfit + xsize)
            ymin = int(yfit - ysize)
            if ymin <= 0:
                ymin = 1
            ymax = int(yfit + ysize)


            # Collect positions of all components
            xposarray = []
            yposarray = []

            # Figure out the nearest neighbors that needed to be 
            # deblended together.

            # component of interest
            sidarray = []
            print offx, offy
            ncomp,xpos,ypos = components(ncomp,i,i,offx,offy) 
            xposarray += [xpos]
            yposarray += [ypos]
            sidarray += [sexid[i]]     # ID's that will be zero-ed out

            # other components
            for j in range(totobj):
                if i != j:
                    dx = x[i] - x[j]
                    dy = y[i] - y[j]
                    dist = sqrt(dx**2 + dy**2)

                    #  Fit all objects that fall within a certain 
                    #  radius of the central obj.:

                    radck = sqrt(isoa[j])
                    angle = arctan2(dy,dx)
                    phi = PI/2. - (theta[j]/180.*PI + PI/2. - angle)
                    thetap = arctan2(axrat[j] * tan(phi), 1.)
                    isorad = radck * sqrt((axrat[j]*cos(thetap))**2 + sin(thetap)**2)

                    # If object is within the fitted image region...
                    if (((x[j]-offx >= xmin) and (x[j]-offx <= xmax) and 
                      (y[j]-offy >= ymin) and (y[j]-offy <= ymax) and
                    # but is greater than 2 Re away, and is not more than 3 mag fainter, or
                      #((dist > 2 * Re[i] and mag[j]-3. <= mag[i]) or 
                      ((dist > 4 * Re[i] and mag[j]-3. <= mag[i]) or
                    # is less than 2 Re away and not more than 5 mag fainter, or
                      #(dist < 2 * Re[i] and mag[j]-5. <= mag[i]))) or
                      (dist < 4 * Re[i]))) or 
                    # is outside image, but is about 3 Re away and 1 mag brighter.
                      ((x[j]-offx >= xmin-isorad) and (x[j]-offx <= xmax+isorad) and 
                      (y[j]-offy >= ymin-isorad) and (y[j]-offy <= ymax+isorad) and
                      mag[j] + MAGDIFF <= mag[i])):
                        print "x[j], offx, xmin, xmax", x[j], offx, xmin, xmax
                        print "y[j], offy, ymin, ymax", y[j], offy, ymin, ymax
                        print "dist, Re[i], mag[j]", dist, Re[i], mag[j]
                        ncomp += 1
                        ncomp,xpos,ypos = components(ncomp,i,j,offx,offy)
                        xposarray += [xpos]
                        yposarray += [ypos]
                        sidarray += [sexid[j]]
                        if Re[j] > Re_max:
                            Re_max = Re[j]
                    elif (dist < (4 * Re[i])):
                        print "Object within 4 Re but not included:", objid[j]
    
            # Zero-out the included components in the mask image

            if MASK == "True":
                ## The condition to make mask
                ## maskexpr = ""
                ## for j in range(len(sidarray)):
                ##     inid = sidarray[j]
                ##     maskexpr = maskexpr + "(a==%d)" % inid
                ##     if j != range(len(sidarray))[-1]:
                ##         maskexpr = maskexpr + "||"
                ## if len(glob.glob("mask?.fits")):
                ##     os.system("rm mask?.fits")
                ## if len(glob.glob("mask.fits")):
                ##     os.system("rm mask.fits")
                ## #if iraf.access("mask.fits"):
                ## #    iraf.imdelete("mask.fits")

                ## # zero-out the objects in sidarray
                ## iraf.imexpr(maskexpr+" ? 0 : a","mask2.fits",a=maskname)

                ## # Now flatten the mask image and convolve with Gaussian
                ## # profile
                ## iraf.imexpr("a>0 ? 100. : 0.","mask3.fits",
                ##     a="mask2.fits")
                ## iraf.gauss(input="mask3.fits",output="mask3.fits",
                ##     sigma=GSIGMA)   
                ## iraf.imexpr("a<=10. ? 0. : a","mask.fits",
                ##     a="mask3.fits")
                ## iraf.imdelete("mask2.fits")
                ## iraf.imdelete("mask3.fits")
                # Use scipy.ndimage to do Gaussian filtering instead of IRAF
                if os.path.exists('mask.fits'):
                    os.remove('mask.fits')
                mask0 = pyfits.getdata(maskname)
                mask0_flat = mask0.flatten()
                mask2_flat = np.where(np.in1d(mask0_flat, sidarray), 0, 
                                      mask0_flat)
                mask2 = mask2_flat.reshape(mask0.shape)
                mask3 = np.where(mask2 > 0, 100., 0.)
                # normalize mask value to 100
                mask4 = ndimage.filters.gaussian_filter(mask3, GSIGMA, order=0)
                mask4 = np.where(mask4 <= 10., 0, mask4)
                hdu = pyfits.PrimaryHDU(mask4)
                hdu.writeto('mask.fits')

        
            # Determine image fitting size and convolution box info 
            #into the galfit template file
            Refrac = 0.7
            if xmin >= min(xposarray) - Refrac*Re_max:
                xmin = min(xposarray) - Refrac*Re_max
            if xmin < 1:
                xmin = 1
            if ymin >= min(yposarray) - Refrac*Re_max:
                ymin = min(yposarray) - Refrac*Re_max
            if ymin < 1:
                ymin = 1
            if xmax <= max(xposarray) + Refrac*Re_max:
                xmax = max(xposarray) + Refrac*Re_max
            if ymax <= max(yposarray) + Refrac*Re_max:
                ymax = max(yposarray) + Refrac*Re_max
            cbox = 60

            print >>OUT1, "H) %d %d %d %d # Image region to fit (xmin xmax ymin ymax) " % (xmin,xmax,ymin,ymax)
            print >>OUT1, "I) %d   %d     # Size of the convolution box (x y) " % (cbox,cbox)

            # Make sky component

            sky(ncomp,bkgnd)
            

            # Calculate information needed to plug back into the image 
            # header at the end of the fit.

            xmin = xmin + offx   # The [xmin:xmax,ymin:ymax] of the box
            xmax = xmax + offx   # relative to the big ACS image from which
            ymin = ymin + offy   # the current tile was extracted
            ymax = ymax + offy

            if MASK == "True":
                print >>OUT2, "%d %s-wmask-out   %d   %d  [%d:%d,%d:%d] %d %d %d " % (sobj,outname,offx,offy,xmin,xmax,ymin,ymax,x[i],y[i],ncomp)
            else:
                print >>OUT2, "%d %s-out   %d   %d  [%d:%d,%d:%d] %d %d %d " % (sobj,outname,offx,offy,xmin,xmax,ymin,ymax,x[i],y[i],ncomp)

            OUT1.close()
            t1 = time.time()
            try:
                galfitrun = subprocess.Popen(["galfit",
                   parfile],stdout=sys.stdout)    # Runs GALFIT!!
                # ======================================================= #
                # Try to time GALFIT and kill it after spending more than 
                # 10 min on one object

                timelimit = 10.*60.   # time limit on one object: 10 minutes
                timeint = 10.    # check time interval
                maxcheck = int(timelimit/timeint)
                kill = 1

                for j in range(maxcheck):
                    time.sleep(timeint)   # sleep for check time interval
                    if galfitrun.poll() == None:   # if it's still running:
                        pass   # haven't reached time limit; let it run
                    else:
                        kill = 0  # finished before time limit; don't need to kill
                        break   # jump out of process checking

                if kill:   # if reached time limit and still running:
                    os.kill(galfitrun.pid, signal.SIGKILL)   # kill GALFIT

                retcode = galfitrun.returncode
                print "retcode", retcode
                t2 = time.time()
                print >>OUT3, "returned code is %d" % retcode
                dt = t2 - t1
                if retcode != 7:    # not sure about the behavior when GALFIT crashes...!!!
                    print >>OUT3, "%d  %s" % (sobj,parfile)
                else:
                    if MASK == "True":
                        print >>OUT4, "%s-wmask-out.fits  %s  %.1f seconds" % (outname,PSFLOC[-1],dt)
                    else:
                        print >>OUT4, "%s-out.fits  %s   %.1f seconds" % (outname,PSFLOC[-1],dt)
            except:
                "Error in running GALFIT..."

        
        #i += 1
        #print "i", i




############# Read in the control parameters as defined above  ###############

if __name__ == "__main__":
    INFILE = sys.argv[1]
    try:
        IN1 = open(INFILE,'r')
    except:
        raise IOError, "Can't open %s !" % INFILE
    l = IN1.readline()
    ACSIMG = l.split()[0]        # The ACS image name to fit.
    l = IN1.readline()
    SEXCAT = l.split()[0]        # The sextractor catalog.
    l = IN1.readline()
    PSFLOC = l.split()[0]        # The directory storing all the PSFs.
    l = IN1.readline()
    STARTNUM = int(l.split()[0])    # Object starting number to fit.
    l = IN1.readline()
    STOPNUM = int(l.split()[0])        # Object stopping number to fit.
    l = IN1.readline()
    CONSFILE = l.split()[0]        # The constraint file name.
    l = IN1.readline()
    TEMPLOC = l.split()[0]        # The directory storing image sub-tiles.
    l = IN1.readline()
    FITFUNC = l.split()[0]        # Fit a single sersic or do B/D decomposition?
    l = IN1.readline()
    REGION = l.split()[0]        # (A) Whole image catalog or, (B) region of an image.
    l = IN1.readline()
    BOUNDARY = l.split()[0]        # Boundary of the object region if (B).
    l = IN1.readline()
    SKYVAL = float(l.split()[0])    # Sky value added back to an image that SExtractor knows nothing about.
    l = IN1.readline()
    MRANGE = l.split()[0]        # Acceptable magnitude range.
    l = IN1.readline()
    STARGAL = l.split()[0]        # Range of good galaxies: 0.0 (galaxy) to 1.0 (star).
    l = IN1.readline()
    FWHMPSF = float(l.split()[0])   # FWHM of the PSF for cosmic ray rejection.
    l = IN1.readline()
    MAGZPT = float(l.split()[0])    # Photometric zeropoint.
    l = IN1.readline()
    PLATE = float(l.split()[0])        # The image plate scale.
    l = IN1.readline()
    NCOL,NROW = int(l.split()[0]),int(l.split()[1])   # The image size of the sub-tiles.
    l = IN1.readline()
    TNAME = l.split()[0]        # The prefix name of the image sub-tiles.
    l = IN1.readline()
    NSPLIT = int(l.split()[0])        # The number of subtiles along one axis.
    l = IN1.readline()
    BUFFER = int(l.split()[0])        # The padding size around the image edges.
    l = IN1.readline()
    MASK = l.split()[0]            # To use mask or not, added by K-H Huang
    l = IN1.readline()
    GSIGMA = float(l.split()[0])    # Sigma of Gaussian wing of masks, added by K-H Huang

    IN1.close()

    if REGION == "B":
        x1,y1,x2,y2 = BOUNDARY[1:-1].split(',')
        #bound = BOUNDARY[1:-1].split(',')
        #for j in range(len(bound)):
        #    bound[j]=int(bound[j])
        #    x1,y1,x2,y2 = bound
    mr = MRANGE[1:-1].split(',')
    brtcut = float(mr[0])
    ftcut = float(mr[1])
    stg = STARGAL[1:-1].split(',')
    cut1 = float(stg[0])
    cut2 = float(stg[1])
    minfwhm = (0.95*FWHMPSF)/0.05   # rids Crays
    xpansz = int(NCOL / NSPLIT)
    ypansz = int(NROW / NSPLIT)

##############################################################################

    print "============================================================="
    print " This script takes the SExtractor output, generated by Dan   "
    print " McIntosh and formats it for GALFIT.  The sextractor output  "
    print " has 18 columns.                            "
    print "============================================================="

    try:
        IN2 = open(SEXCAT)
    except:
        raise IOError, "Can't open %s !" % SEXCAT

###############  Read in sextractor data, weed out bad objs  #################

    i = 1    # good array element
    ii = 0    # overall counting

    OUT2 = open(TEMPLOC+'/offsets','wb')
    OUT2.write('#  ID  output_image  offx  offy  xyrange  x  y  ncomp\n')
    lines = IN2.readlines()
    objid=[]
    flux=[]
    fluxerr=[]
    mag=[]
    magerr=[]
    Re=[]
    bkg=[]
    isoa=[]
    x=[]
    y=[]
    alpha=[]
    delta=[]
    fwhm=[]
    axrat=[]
    flag=[]
    index=[]
    theta=[]
    sexid=[]
    for l in lines:
        l = l.split('\n')[0]    # as chomp in perl
        (num,f,fe,mg,mgerr,rk,bkgd,ia,xx,yy,ra,dec,pa,ell,fw,flg,
            idx,sid)=l.split()
        if num != "#":
            num = int(num)
            f,fe = float(f),float(fe)
            mg,mgerr = float(mg),float(mgerr)
            rk,bkgd = float(rk),float(bkgd)
            ia = float(ia)
            xx,yy = float(xx),float(yy)
            ra,dec = float(ra),float(dec)
            pa,ell = float(pa),float(ell)
            fw = float(fw)
            flg = int(flg)
            idx = float(idx)  # CLASS_STAR
            sid = int(sid)
        
            ii += 1
            # std. flag and cosmic ray culling
            # objid, etc. are arrays for ALL objects in SEXCAT, not
            # just the artificial ones
            #if flg < 4 and fw > minfwhm:
            if idx <= cut2 and idx >= cut1:
                if REGION == "A" or (REGION == "B" and xx >= x1 and xx <= x2 and yy >= y1 and yy <= y2):
                    objid += [num]
                    flux += [f]
                    fluxerr += [fe]
                    mag += [mg]
                    magerr += [mgerr]
                    Re += [0.3*sqrt(ia)]
                    bkg += [bkgd]
                    isoa += [ia]
                    x += [xx]
                    y += [yy]
                    alpha += [ra]
                    delta += [dec]
                    fwhm += [fw]
                    axrat += [1.-ell]
                    flag += [flg]
                    index += [idx]
                    sexid += [sid]
                    # Sextractor to galfit orientation
                    theta += [pa - 90.]
                    i += 1
    objid = array(objid)
    flux = array(flux)
    fluxerr = array(fluxerr)
    mag = array(mag)
    magerr = array(magerr)
    Re = array(Re)
    bkg = array(bkg)
    isoa = array(isoa)
    x = array(x)
    y = array(y)
    alpha = array(alpha)
    delta = array(delta)
    fwhm = array(fwhm)
    axrat = array(axrat)
    flag = array(flag)
    index = array(index)
    sexid = array(sexid)
    theta = array(theta)
    totobj = len(objid)

    OUT3 = open("galfit.crashes","w")
    OUT4 = open("galfit.fitted","a")

    run_galfit()

    OUT4.close()
    OUT3.close()
    OUT2.close()
    try:
        IN2.close()
    except:
        raise IOError, "Can't close %s !" % SEXCAT



