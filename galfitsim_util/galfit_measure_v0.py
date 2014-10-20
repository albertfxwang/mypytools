#!/usr/bin/env python

import numpy as N
from pygoods import sextractor, Ftable
from pyraf import iraf
from iraf import images, stsdas
import os, glob, subprocess, sys, time, signal
import pyfits
import pywcs

class galfit_measure(object):
    """
    Runs GALFIT measurements for the specified objects in a given catalog.
    Makes image cutouts (including RMS maps and seg maps) automatically.
    """
    def __init__(self,id2fit,sexcat,sci_image,rms_image,seg_image,
        sexcatfmt='fits',idcol='id_1',racol='ra_1',deccol='dec_1',
        magcol='mag_auto',recol='flux_radius_2',axratiocol='axratio',thetacol='theta_image',
        isoareacol='isoarea_image',root='images',pixscale=0.06,psffile="",
        constraint_file="constraints3",magzpt=25.96,cbox=160):
        """
        id2fit: array of ID numbers to run GALFIT on.
        sexcat: the SExtractor catalog where the object IDs are specified.
        NOTE: sexcat needs to contain ALL objects in the field, not just the objects to 
        be fit by GALFIT!!
        sci_image: the science (drizzled) mosaic
        rms_image: the RMS image
        seg_image: the segmentation map
        racol: column name that contains RA
        deccol: column name that contains DEC
        magcol: column name that contains magnitude (default to MAG_AUTO)
        recol: column name for Re (half-light radius)
        axratiocol: column name for axis ratio
        thetacol: column name for position angle (theta)
        root: the directory name for all the GALFIT cutouts and output images
        pixscale: the pixel scale in arcsec
        """
        self.id2fit = id2fit
        if sexcatfmt=='fits':
            self.sexcat = Ftable(sexcat)
            self.nobj_all = len(self.sexcat.d)
        else:
            self.sexcat = sextractor(sexcat)
            self.nobj_all = len(self.sexcat)
        if not os.path.exists(root):
            os.mkdir(root)
        elif not os.path.isdir(root):
            raise OSError, "%s already exists and is not a directory." % root
        self.ra_array = getattr(self.sexcat,racol)
        self.dec_array = getattr(self.sexcat,deccol)
        self.id_array = getattr(self.sexcat,idcol)
        self.mag_array = getattr(self.sexcat,magcol)
        self.Re_array = getattr(self.sexcat,recol)
        self.axratio_array = getattr(self.sexcat,axratiocol)
        self.theta_array = getattr(self.sexcat,thetacol)
        self.isoarea_array = getattr(self.sexcat,isoareacol)

        self.sci_image = sci_image
        self.rms_image = rms_image
        self.seg_image = seg_image
        # figure out the WCS
        hdr = pyfits.getheader(self.sci_image)
        self.wcs = pywcs.WCS(hdr)
        # figure out the X, Y positions of each source... will need later
        skycrd = N.array([self.ra_array,self.dec_array])
        skycrd = skycrd.swapaxes(0,1)
        pixcrd = self.wcs.wcs_sky2pix(skycrd,1)  # pixel coordinate starts at 1
        pixcrd = pixcrd.swapaxes(0,1)
        self.x_array = pixcrd[0]
        self.y_array = pixcrd[1]
        self.crashes = "galfit.crashes"
        self.fitted = "galfit.fitted"
        self.root = root  # a directory containing all the cutouts
        self.xscale = pixscale
        self.yscale = pixscale
        self.psffile = psffile
        self.constraint_file = constraint_file
        self.magzpt = magzpt
        self.cbox = cbox  # convolution box size --- should be comparable to the PSF image size?

    def fitall(self):
        """
        Fit all objects in self.id2dit
        Only works if all objects have the same PSF file and magnitude zero-point...
        """
        for objectid in self.id2fit:
            self.fitone(objectid,self.psffile,self.magzpt)

    def fitmany(self, id_array):
        """
        Fit many objects.
        """
        for objectid in id_array:
            self.fitone(objectid,self.psffile,self.magzpt)
            
    def fitone(self, objectid, psffile=None, magzpt=None, exptime=1.0):
        """
        A driver function that fits one object, given its ID.
        Most of the work is done here.
        Provide a PSF file and magnitude zero point for each object!
        """
        if psffile==None: psffile=self.psffile
        if magzpt==None: magzpt=self.magzpt
        idcrit = (self.id_array==objectid)
        ra = self.ra_array[idcrit][0]
        dec = self.dec_array[idcrit][0]
        mag = self.mag_array[idcrit][0]
        Re = self.Re_array[idcrit][0]
        axratio = self.axratio_array[idcrit][0]
        theta = self.theta_array[idcrit][0]
        isoarea = self.isoarea_array[idcrit][0]
        # Now figure out image coordinates
        x = self.x_array[idcrit][0]
        y = self.y_array[idcrit][0]
        print "x, y =", x, y
        OUT1 = open(self.crashes,"ab")
        OUT2 = open(self.fitted,"ab")

        # now create arrays that will be included in this GALFIT run
        id_fit = [objectid]
        mag_fit = [mag]
        Re_fit = [Re]
        x_fit = [x]
        y_fit = [y]
        theta_fit = [theta]
        axratio_fit = [axratio]
        
        # step 1: determine temporary image size
        # Experimental image size
        majx = abs(N.sqrt(isoarea) * 4. * N.sin(theta/180.*N.pi))  # convert degree into radian
        majy = abs(N.sqrt(isoarea) * 4. * N.cos(theta/180.*N.pi))
        minx = axratio * N.sqrt(isoarea) * 4. * abs(N.sin((theta+N.pi/2.)/180.*N.pi))
        miny = axratio * N.sqrt(isoarea) * 4. * abs(N.sin((theta+N.pi/2.)/180.*N.pi))
        #Re_max = 20.
        Re_max = 60
        if Re > Re_max:
            Re_max = Re
        xsize = majx
        # determine the offset of this cutout relative to the entire mosaic
        if minx >= majx:
            xsize = minx
        if xsize <= 30.:
            xsize = 30.
        ysize = majy
        if miny >= majy:
            ysize = miny
        if ysize <= 30.:
            ysize = 30.
        xmin = int(x - xsize)
        if xmin <= 0:
            xmin = 1
        xmax = int(x + xsize)
        ymin = int(y - ysize)
        if ymin <= 0:
            ymin = 1
        ymax = int(y + ysize)
        # construct postage stamp name
        outname = self.root+'/obj%d_out.fits' % objectid
        
        # step 2: determine which objects to be included
        for jj in range(self.nobj_all):  # iterate through other objects
            if self.id_array[jj]!=objectid:
                dx = x - self.x_array[jj]
                dy = y - self.y_array[jj]
                dist = N.sqrt(dx**2 + dy**2)
                # Fit all objects that fall within a certain 
                #  radius of the central object:
                radck = N.sqrt(self.isoarea_array[jj])
                angle = N.arctan2(dy,dx)
                phi = N.pi/2. - (self.theta_array[jj]/180.*N.pi + N.pi/2. - angle)
                thetap = N.arctan2(self.axratio_array[jj] * N.tan(phi), 1.)
                isorad = radck * N.sqrt((self.axratio_array[jj]*N.cos(thetap))**2 + N.sin(thetap)**2)
                # If object is within the fitted image region...
                if (((self.x_array[jj] >= xmin) and (self.x_array[jj] <= xmax) and
                  (self.y_array[jj] >= ymin) and (self.y_array[jj] <= ymax) and
                  # but is greater than 2 Re away, and is not more than 3 mag fainter, or
                  ((dist > 2 * Re and self.mag_array[jj]-3. <= mag) or
                  # is less than 2 Re away and not more than 5 mag fainter, or
                  (dist < 2. * Re and self.mag_array[jj]-5. <= mag))) or
                  # is outside image, but is about 3 Re away and 1 mag brighter
                  ((self.x_array[jj] >= (xmin-isorad) and (self.x_array[jj] <= (xmax+isorad) and
                  (self.y_array[jj] >= (ymin - isorad)) and (self.y_array[jj] <= (ymax+isorad)) and
                  (self.mag_array[jj] + 1.) <= mag)))):
                    id_fit += [self.id_array[jj]]
                    x_fit += [self.x_array[jj]]
                    y_fit += [self.y_array[jj]]
                    mag_fit += [self.mag_array[jj]]
                    Re_fit += [self.Re_array[jj]]
                    axratio_fit += [self.axratio_array[jj]]
                    theta_fit += [self.theta_array[jj]]
                    if self.Re_array[jj] > Re_max:
                        Re_max = self.Re_array[jj]
        
        # step 3: determine fitting size & copy images
        # Determine image fitting size and convolution box info
        # into the galfit template file
        Refrac = 0.7
        #Refrac = 2.0
        if xmin >= min(x_fit) - Refrac*Re_max:
            xmin = min(x_fit) - Refrac*Re_max
        if xmin < 1:
            xmin = 1
        if ymin >= min(y_fit) - Refrac*Re_max:
            ymin = min(y_fit) - Refrac*Re_max
        if ymin < 1:
            ymin = 1
        if xmax <= max(x_fit) + Refrac*Re_max:
            xmax = max(x_fit) + Refrac*Re_max
        if ymax <= max(y_fit) + Refrac*Re_max:
            ymax = max(y_fit) + Refrac*Re_max
        #cbox = 60
        imgxsize = xmax-xmin+1
        imgysize = ymax-ymin+1
        if os.path.exists(outname):
            os.remove(outname)
        # copy drizzled images
        if os.path.exists(self.root+'/obj%d_drz.fits' % objectid):
            os.system('rm %s/obj%d_drz.fits' % (self.root,objectid))
        drzimg = self.sci_image+'[%d:%d,%d:%d]'%(xmin,xmax,ymin,ymax)
        iraf.imcopy(drzimg,self.root+'/obj%d_drz.fits' % objectid)
        h = pyfits.open(self.root+'/obj%d_drz.fits' % objectid,'update')
        h[0].header.update('xmin',xmin)
        h[0].header.update('xmax',xmax)
        h[0].header.update('ymin',ymin)
        h[0].header.update('ymax',ymax) 
        h[0].header.update('exptime',exptime)
        h.flush()
        h.close()
        # copy sigma image
        if os.path.exists(self.root+'/obj%d_rms.fits' % objectid):
            os.system('rm %s/obj%d_rms.fits' % (self.root,objectid))
        iraf.imcopy(self.rms_image+'[%d:%d,%d:%d]'%(xmin,xmax,ymin,ymax),self.root+'/obj%d_rms.fits' % objectid)
        # copy segmentation images
        if os.path.exists(self.root+'/obj%d_seg.fits' % objectid):
            os.system('rm %s/obj%d_seg.fits' % (self.root,objectid))
        iraf.imcopy(self.seg_image+'[%d:%d,%d:%d]'%(xmin,xmax,ymin,ymax),self.root+'/obj%d_seg.fits' % objectid)
        # make mask image -- zero-out the included components in the mask image
        if os.path.exists("mask2.fits"):
            iraf.imdelete("mask2.fits")
        if os.path.exists("mask3.fits"):
            iraf.imdelete("mask3.fits")
        maskexpr = ""
        maskname = self.root+"/obj%d_seg.fits" % objectid
        print "id_fit", id_fit
        for j in range(len(id_fit)):
            inid = id_fit[j]
            maskexpr = maskexpr + "(a==%d)" % inid
            if j != len(id_fit)-1:
                maskexpr = maskexpr + "||"
        print maskexpr   
        # Now make the masks 
        # In GALFIT, good pixels have mask value 0, while bad pixels have mask values > 0
        iraf.imexpr(maskexpr+" ? 0 : a","mask2.fits",a=maskname)
        # After this step, the only non-zero pixels will be those belonging to the objects that are NOT
        # included in the list id_fit.
        # Now flatten the mask image and convolve with a Gaussian, in order to grow the masked region a little bit
        # to account for wings of objects.
        iraf.imexpr("a>0 ? 100. : 0.","mask3.fits",a="mask2.fits")
        iraf.gauss(input="mask3.fits",output="mask3.fits",sigma=5)
        if os.path.exists(self.root+'/mask_obj%d.fits'%objectid):
            os.remove(self.root+'/mask_obj%d.fits'%objectid)
        iraf.imexpr("a<=10. ? 0. : a",self.root+"/mask_obj%d.fits" % objectid,
          a="mask3.fits")
        os.remove("mask2.fits")
        os.remove("mask3.fits")
        
        # step 4: write GALFIT input file
        parfile = "obj%d" % objectid
        f = open(parfile,'w')
        # Write image parameters
        print >>f, "==============================================================================="
        print >>f, "# IMAGE PARAMETERS"
        print >>f, "A) %s/obj%d_drz.fits   # Input Data image (FITS file) " % (self.root,objectid)
        print >>f, "B) %s/obj%d_out.fits   # Output data image block " % (self.root,objectid)
        print >>f, 'C) %s/obj%d_rms.fits   # Noise image name (made from data if blank or "none") ' % (self.root,objectid)
        print >>f, "D) %s         # Input PSF image for convolution (FITS file) " % psffile
        print >>f, "E) 1                   # PSF oversampling factor relative to data "
        print >>f, "F) %s/mask_obj%d.fits     # Bad pixel mask (FITS image or ASCII coord list)" % (self.root,objectid) 
        print >>f, "G) %s           # Parameter constraint file " % self.constraint_file
        print >>f, "H) %d %d %d %d  # Image region to fit " % (1,imgxsize,1,imgysize)
        print >>f, "I) %d  %d       # Size of the convolution box (x y)" % (self.cbox,self.cbox)
        print >>f, "J) %f             # Magnitude photometric zeropoint " % magzpt
        print >>f, "K) %s %s       # Plate scale (dx dy). " % (self.xscale,self.yscale)
        print >>f, "O) regular             # Display type (regular, curses, both) "
        print >>f, "P) 0                   # Create ouput only? (1=yes; 0=optimize) "
        print >>f, "S) 0                   # Modify/create objects interactively? "
        print >>f, ""
        print >>f, "# INITIAL FITTING PARAMETERS "
        print >>f, ""
        print >>f, "#   For object type, allowed functions are: sersic, nuker, "
        print >>f, "#                       expdisk, devauc, moffat, gaussian. "
        print >>f, ""
        print >>f, "# Objtype:      Fit?         Parameters "
        print >>f, ""
        # Write individual component
        for j in range(len(id_fit)):
            self.component(f,j+1,x_fit[j]-xmin+1,y_fit[j]-ymin+1,mag_fit[j],
              Re_fit[j],axratio_fit[j],theta_fit[j])
        # Make sky component
        #bkgnd = skyest('images/obj%d_drz.fits'%objectid,'images/obj%d_seg.fits'%objectid,5.)
        bkgnd = 0.
        self.sky(f,len(id_fit)+1,bkgnd)
        f.close()

        # step 5: run GALFIT
        t1 = time.time()
        galfitrun = subprocess.Popen(["galfit",parfile],stdout=sys.stdout)
        # runs GALFIT
        # ======================================================= #
        # Try to time GALFIT and kill it after spending more than 
        # 10 min on one object
        # ======================================================= #
        timelimit = 10. * 60. # time limit: 10 minutes
        timeint = 10.   # check every 10 seconds
        maxcheck = int(timelimit/timeint)
        kill = 1
        for j in range(maxcheck):
            time.sleep(timeint)  # sleep for check time interval
            if galfitrun.poll() == None:   # if it's still running:
                pass   # haven't reached time limit; let it run
            else:
                kill = 0  # finished before time limit; don't need to kill
                break   # jump out of process checking
        if kill:   # if reached time limit and still running:
            os.kill(galfitrun.pid, signal.SIGKILL)   # kill GALFIT
            print "Killed GALFIT; GALFIT froze"
        retcode = galfitrun.returncode
        t2 = time.time()
        try:
            print >>OUT1, "returned code is %d" % retcode
        except:
            print >>OUT1, "returned code is ", retcode
        dt = t2 - t1
        if retcode != 7:    # not sure about the behavior when GALFIT crashes...!!!
            print >>OUT1, "%d  %s" % (objectid,parfile)
        else:
            print >>OUT2, "%s-wmask-out.fits %.1f seconds" % (outname,dt)
        OUT1.close()
        OUT2.close()
        print "Finishes GALFIT fitting"

    def sky(self,f,ncomp,bkgnd):
        # Make sky component
        print >>f, "# Object number: %d" % ncomp
        print >>f, " 0)    sky            #  Object type"
        print >>f, " 1)    %f   1         #  sky brightness" % bkgnd
        print >>f, " 2)    0.   0         #  sky gradient in x" 
        print >>f, " 3)    0.   0         #  sky gradient in y" 
        print >>f, ""
        print >>f, "==============================================================================="

    def component(self,f,ncomp,x,y,mag,re,axratio,pa):
        print >>f, "# Object number: %d" % ncomp
        print >>f, " 0)     sersic      # Object type"
        print >>f, " 1)  %f  %f  1  1   # position x, y" % (x,y)
        print >>f, " 3)  %f      1      # integrated magnitude" % mag
        print >>f, " 4)  %f      1      # R_e (half-light radius)   [pix]" % re
        print >>f, " 5)  1.5     1      # Sersic index n (de Vaucouleurs n=4)"
        print >>f, " 9)  %f      1      # axis ratio (b/a)" % axratio
        print >>f, " 10) %f      1      # position angle (PA) [deg: Up=0, Left=90]" % pa

def skyest(inputimg,segimg,gsigma):
    #if len(glob.glob("cutout.fits")):
    #    os.system("rm cutout.fits")
    #iraf.imcopy("%s[%d:%d,%d:%d]"%(inputimg,xmin,xmax,ymin,ymax),"cutout.fits")
    #iraf.imcopy(segimg+"[%d:%d,%d:%d]"%(xmin,xmax,ymin,ymax),"segcutout.fits")
    if len(glob.glob('segcutout*.fits')):
        os.system('rm segcutout*.fits')
    iraf.imexpr("a > 0 ? 100 : 0","segcutout2.fits",a=segimg)
    # grow mask by Gaussian wing
    iraf.gauss(input="segcutout2.fits",output="segcutout3.fits",sigma=gsigma)
    h1 = pyfits.open(inputimg)
    h2 = pyfits.open('segcutout3.fits')
    him1 = h1[0].data  # image data
    him2 = h2[0].data
    skyim = compress(him2.ravel()<10.,him1.ravel())
    bkgnd = median(skyim)  # take the median of unmasked sky
    #os.system('rm cutout.fits')
    os.system('rm segcutout*.fits')
    return bkgnd


#if __name__ == "__main__":
#    parfile = sys.argv[1]
#    c = parseconfig(parfile)
#    startnum = c['STARTNUM']
#    stopnum = c['STOPNUM']
#    band = c['BAND']
#    consfile = c['CONSFILE']
#    dropband = c['DROPBAND']
#    dropcat = c['DROPCAT']
#    magzpt = c['MAGZPT']
#    psffile = c['PSFFILE']
#    osectdir_n = c['OSECTDIR_N']
#    osectdir_s = c['OSECTDIR_S']
#    fit_LBG(startnum,stopnum,osectdir_n,osectdir_s,psffile,consfile,magzpt,dropcat,band,dropband)
