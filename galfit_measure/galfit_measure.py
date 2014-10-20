#!/usr/bin/env python

from numpy import *
import numpy as np
from pygoods import *
from pyraf import iraf
from iraf import images, stsdas
import os, glob, subprocess, sys, time, signal
import pyfits
import yaml
import pywcs
from scipy import ndimage
from ImageTools.pyraf_utils import imcopy

re_factor = 2.0   
# automatically include a neighbor if it is within Re*re_factor

def sky(f,ncomp,bkgnd):
    # Make sky component
    
    print >>f, "# Object number: %d" % ncomp
    print >>f, " 0)    sky            #  Object type"
    print >>f, " 1)    %f   1         #  sky brightness" % bkgnd
    print >>f, " 2)    0.   0         #  sky gradient in x" 
    print >>f, " 3)    0.   0         #  sky gradient in y" 
    print >>f, ""
    print >>f, "==============================================================================="

def component(f,ncomp,x,y,mag,re,axratio,pa):
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

def run_galfit(paramfile):
    """
    Runs GALFIT measurement in batch mode over all (or a subset of) objects 
    in the given catalog. Note that the x, y positions of each object need to
    be consistent with the provided input images.
    Provide a parameter file that include the following information:
    - startnum: which object to start (not object ID but the order of entries
                                       in the catalog)
    - stopnum: which object to stop (including this object)
    - drzimage: the image over which to run GALFIT measurement
    - rmsimage: the RMS map
    - segimage: the segmentation map (used for masking neighbors)
    **IMPORTANT: we assume that drzimage, rmsimage, and segimage are registered
    on to the same WCS grid, with the same image size. Otherwise proper masks
    for each object cannot be made.
    - psfimage: the PSF image
    - magzpt: the magnitude zero-point
    - consfile: GALFIT parameter constraint file
    - pixscale: pixel scale of the images in arcsec (does not really matter...?)
    - galfitpath: the path to the executable GALFIT code
    - target_catalog: the SExtractor catalog containing the targets to fit
    - all_catalog: the SExtractor catalog containing ALL objects in the field,
                   in order to include also neighbors in the fit
    *** Make sure that target_catalog and all_catalog have matching IDs! ***
    - band: filter name of the image
    - id_column: the column name for object ID (usually either ID_MOSAIC or NUMBER)
    - ra_column: the column name for right ascension
    - dec_column: the column name for delination
    - re_column: the column name for the SExtractor-measured half-light radius
    - mask: whether to mask neighbors or not (should be true unless in special circumstances)
    - mag_column: the column name for the best magnitude (in the measurement band)
    - mag_column
    - x_column
    - y_column
    - ell_column
    - theta_column
    - isoa_column
    - overwrite: whether to re-run previously measured objects
    - nblocks_x & nblocks_y: number of subpanels in each direction
    - paddiing: the width of border (in pixels) that each subpanel extends
                over its nominal limit, so that galaxies at borders of each
                subpanel will not be cut-off
    - split_image: whether to split the input image, or use existing subpanels
      If split_image == no, it will look for subpanels in the directory 
      subpanels/input_image_[type]_[XX]_[YY].fits, where [type] could be either
      drz or rms, and [XX], [YY] are the subpanel numbers in each dimension.

    Future developments:
    - The ability to use different PSFs for different regions in the mosaic
    """
    print "Entering galfit_measure.run_galfit()..."
    pars = yaml.load(open(paramfile, 'rb'))
    curdir = os.getcwd()
    if pars['target_catalog'].endswith('.fits'):
        c = Ftable(pars['target_catalog'])
        Ncat_target = len(c.d)
    else:
        c = sextractor(pars['target_catalog'])
        # Ncat = len(c)
        c.__getitem__ = c.__getattribute__
        # print "Please give a FITS table."
        # sys.exit()
        Ncat_target = len(c)
    if pars['all_catalog'].endswith('.fits'):
        c_all = Ftable(pars['all_catalog'])
        Ncat = len(c_all.d)
    else:
        c_all = sextractor(pars['all_catalog'])
        c_all.__getitem__ = c_all.__getattribute__
        Ncat = len(c_all)
    if pars['stopnum'] > Ncat_target:
        pars['stopnum'] = Ncat_target
    object_id_target = c.__getitem__(pars['id_column'].lower())
    object_id_all = c_all.__getitem__(pars['id_column'].lower())
    ra = c_all.__getitem__(pars['ra_column'].lower())
    dec = c_all.__getitem__(pars['dec_column'].lower())
    radec = np.array([ra, dec]).swapaxes(0, 1)
    OUT3 = open("galfit.crashes","w")
    OUT4 = open("galfit.fitted","a")
    # Load mask image
    # if pars['mask'] == True:
    #     maskname = pars['segimage']
    #     mask0 = pyfits.getdata(maskname)
    GSIGMA = 10.
    if not os.path.isdir('galfit_output'):
        os.mkdir('galfit_output')
    if not os.path.isdir('galfit_input'):
        os.mkdir('galfit_input')
    mag = c_all.__getitem__(pars['mag_column'])
    Re = c_all.__getitem__(pars['re_column'].lower())
    ellip = c_all.__getitem__(pars['ell_column'])
    axratio = 1.- ellip
    theta_sex = c_all.__getitem__(pars['theta_column'])
    theta = theta_sex + 90.
    # osect = acscat.osect_refnum
    isoa = c_all.__getitem__(pars['isoa_column'])
    # Now split up the large image into several pieces if necessary, and then
    # find which sub-panel each source lies in.
    hdr_input = pyfits.getheader(pars['drzimage'])
    ncols0 = hdr_input['naxis1']
    nlines0 = hdr_input['naxis2']
    # x = c_all.__getitem__(pars['x_column'])
    # y = c_all.__getitem__(pars['y_column'])
    # Calculate (x, y) from (RA, DEC) 
    wcs_input = pywcs.WCS(hdr_input)
    xy0 = wcs_input.wcs_sky2pix(radec, 1)
    x_input = xy0[:,0]
    y_input = xy0[:,1]
    # Now split the images if necessary, and calculate the shifts in x and y 
    # in each subpanel
    # x0_shift = {}
    # y0_shift = {}
    mask0_dic = {}
    x_image = {}   # store the pixel coordinates of all objects in all_catalog
    y_image = {}
    if (pars['nblocks_x'] > 1) or (pars['nblocks_y'] > 1):
        x0_member = np.linspace(1, ncols0, pars['nblocks_x'] + 1)
        x0_member = np.around(x0_member).astype('int')
        y0_member = np.linspace(1, nlines0, pars['nblocks_y'] + 1)
        y0_member = np.around(y0_member).astype('int')
        # ix, iy are the indices of subpanels for each galaxy in all_catalog
        ix = np.searchsorted(x0_member, x_input)
        iy = np.searchsorted(y0_member, y_input)
        ## TO DO: cut images into pieces...
        ## Be careful about how I handle segmentation map... but worry about it
        ## later
        for nx in range(pars['nblocks_x']):
            for ny in range(pars['nblocks_y']):
                x0 = np.maximum(1, x0_member[nx]-pars['padding'])
                x1 = np.minimum(ncols0, x0_member[nx+1]+pars['padding'])
                y0 = np.maximum(1, y0_member[ny]-pars['padding'])
                y1 = np.minimum(nlines0, y0_member[ny+1]+pars['padding'])
                # x0_shift[(nx+1, ny+1)] = x0 - 1
                # y0_shift[(nx+1, ny+1)] = y0 - 1
                if pars['split_image']:
                    if not os.path.isdir('subpanels'):
                        os.mkdir('subpanels')
                    print "Copying subpanel (%d, %d): [%d:%d,%d:%d]" % (nx+1, ny+1, x0, x1, y0, y1)
                    imcopy(pars['drzimage'], [x0, x1, y0, y1],
                          'subpanels/input_drz_%d_%d.fits' % (nx+1, ny+1))
                    imcopy(pars['rmsimage'], [x0, x1, y0, y1],
                          'subpanels/input_rms_%d_%d.fits' % (nx+1, ny+1))
                    imcopy(pars['segimage'], [x0, x1, y0, y1],
                          'subpanels/input_seg_%d_%d.fits' % (nx+1, ny+1))
                mask0_dic[(nx+1, ny+1)] = pyfits.getdata('subpanels/input_seg_%d_%d.fits' % (nx+1, ny+1))
                # Now calculate the pixel coordinates of all objects in this
                # subpanel
                hdr_subpanel = pyfits.getheader('subpanels/input_drz_%d_%d.fits' % (nx+1, ny+1))
                wcs_subpanel = pywcs.WCS(hdr_subpanel)
                print "Calculate (x, y) in subpanel (%d, %d) from WCS..." % (nx+1, ny+1)
                xy_subpanel = wcs_subpanel.wcs_sky2pix(radec, 1)  # 1-based index
                x_image[(nx+1, ny+1)] = xy_subpanel[:,0]
                y_image[(nx+1, ny+1)] = xy_subpanel[:,1]

    else:
        ix = np.ones(len(x_input), 'int')
        iy = np.ones(len(y_input), 'int')
        x_image[(1,1)] = x_input
        y_image[(1,1)] = y_input
        if os.path.exists('subpanels'):
            os.system('rm subpanels/*.fits')
        else:
            os.mkdir('subpanels')
        os.symlink(os.path.join(curdir, pars['drzimage']), 'subpanels/input_drz_1_1.fits')
        os.symlink(os.path.join(curdir, pars['rmsimage']), 'subpanels/input_rms_1_1.fits')
        os.symlink(os.path.join(curdir, pars['segimage']), 'subpanels/input_seg_1_1.fits')
        # x0_shift[(1,1)] = 0
        # y0_shift[(1,1)] = 0
        mask0_dic[(1,1)] = pyfits.getdata('subpanels/input_seg_1_1.fits')

    for i in range(pars['startnum'], pars['stopnum']):
        # i: index in c
        objid = object_id_target[i]
        ixy = (ix[object_id_all==objid][0], iy[object_id_all==objid][0])
        # print "ixy", ixy
        input_drz = 'subpanels/input_drz_%d_%d.fits' % (ixy[0], ixy[1])
        input_rms = 'subpanels/input_rms_%d_%d.fits' % (ixy[0], ixy[1])
        hdr = pyfits.getheader(input_drz)
        ncols = hdr['naxis1']
        nlines = hdr['naxis2']
        if pars['overwrite'] == False:
            if os.path.exists('galfit_output/obj%d_out.fits' % objid):
                print "Skipping object ID %d" % objid
                continue
        print "Fitting object ID %d" % objid
        ## TO DO: figure out which subpanel each target is in, and calculate
        ## the pixel coordinates in the subpanel
        # ii: index of main object in acscat        
        sidarray = [object_id_target[i]]
        magarray = [mag[object_id_all==objid][0]]
        rearray = [Re[object_id_all==objid][0]]
        xposarray = [x_image[ixy][object_id_all==objid][0]]
        if (xposarray[0] < 0) | (xposarray[0] > ncols):
            print "Object %d is outside the footprint of this image." % objid
            continue
        yposarray = [y_image[ixy][object_id_all==objid][0]]
        if (yposarray[0] < 0) | (yposarray[0] > nlines):
            print "Object %d is outside the footprint of this image." % objid
            continue
        theta_array = [theta[object_id_all==objid][0]]
        axrat_array = [axratio[object_id_all==objid][0]]
        segnum_array = [object_id_target[i]]
        ncomp = 1

        # step 1: determine temporary image size   
        # Experimental image size
        majx = abs(sqrt(isoa[object_id_all==objid][0]) * 4. * sin(theta[object_id_all==objid][0]/180.*pi))
        # convert degree into radian
        majy = abs(sqrt(isoa[object_id_all==objid][0]) * 4. * cos(theta[object_id_all==objid][0]/180.*pi))
        minx = axratio[object_id_all==objid][0] * sqrt(isoa[object_id_all==objid][0]) * 4. * abs(sin((theta[object_id_all==objid][0]+pi/2.)/180.*pi))
        miny = axratio[object_id_all==objid][0] * sqrt(isoa[object_id_all==objid][0]) * 4. * abs(sin((theta[object_id_all==objid][0]+pi/2.)/180.*pi))
        Re_max = 20.
        if Re[object_id_all==objid][0] > Re_max:
            Re_max = Re[object_id_all==objid][0]

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
        # print "xsize", xsize
        # print "ysize", ysize
        xmin = int(xposarray[0] - xsize)
        if xmin <= 0:
            xmin = 1
        xmax = int(xposarray[0] + xsize)
        ymin = int(yposarray[0] - ysize)
        if ymin <= 0:
            ymin = 1
        ymax = int(yposarray[0] + ysize)
        print "xmin, xmax, ymin, ymax", xmin, xmax, ymin, ymax
        # construct postage stamp name
        outname = 'obj%d_out.fits' % objid
        if os.path.exists('galfit_output/%s' % outname):
            os.remove('galfit_output/%s' % outname)

        # step 2: determine which objects to be included
        for j in range(Ncat):  # iterate through other objects
            # if (j != i):
            if object_id_all[j] != objid:
                dx = xposarray[0] - x_image[ixy][j]
                dy = yposarray[0] - y_image[ixy][j]
                dist = sqrt(dx**2 + dy**2)
 
                # Fit all objects that fall within a certain 
                #  radius of the central obj.:
                radck = sqrt(isoa[j])
                angle = arctan2(dy,dx)
                phi = pi/2. - (theta[j]/180.*pi + pi/2. - angle)
                thetap = arctan2(axratio[j] * tan(phi), 1.)
                isorad = radck * sqrt((axratio[j]*cos(thetap))**2 + sin(thetap)**2)
 
                # If object is within the fitted image region...
                if (((x_image[ixy][j] >= xmin) and (x_image[ixy][j] <= xmax) and
                  (y_image[ixy][j] >= ymin) and (y_image[ixy][j] <= ymax) and
                  # but is greater than 2 Re away, and is not more than 3 mag fainter, or
                  #((dist > 2 * Re[object_id_all==objid][0] and mag[j]-3. <= mag[object_id_all==objid][0]) or
                  ((dist > re_factor * rearray[0] and mag[j]-3. <= magarray[0]) or
                  # is less than 2 Re away and not more than 5 mag fainter, or
                  #(dist < 2. * Re[object_id_all==objid][0] and mag[j]-5. <= mag[object_id_all==objid][0]))) or
                  (dist < re_factor * Re[object_id_all==objid][0]))) or 
                  # is outside image, but is about 3 Re away and 1 mag brighter
                  ((x_image[ixy][j] >= xmin - isorad) and (x_image[ixy][j] <= xmax + isorad) and
                  (y_image[ixy][j] >= ymin - isorad) and (y_image[ixy][j] <= ymax + isorad) and
                  mag[j] + 1. <= magarray[0])):
                    xposarray += [x_image[ixy][j]]
                    yposarray += [y_image[ixy][j]]
                    sidarray += [object_id_all[j]]
                    magarray += [mag[j]]
                    rearray += [Re[j]]
                    axrat_array += [axratio[j]]
                    theta_array += [theta[j]]
                    segnum_array += [object_id_all[j]]
                    if Re[j] > Re_max:
                        Re_max = Re[j]
    
        # step 3: determine fitting size & copy images
        # Determine image fitting size and convolution box info
        # into the galfit template file
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
        imgxsize = xmax-xmin+1
        imgysize = ymax-ymin+1

        ## No need to make cutouts myself... Just specify the image region to fit?
        # copy drizzled images
        # if iraf.access('images/obj%d_drz.fits' % objid):
        #     os.system('rm images/obj%d_drz.fits' % objid)
        # if (field=='s') & (band=='i') & (obj_osect==12):
        #     drzimg = osectdir_s+'/%s%s_v1.9_sc03_osect%d_drz.fits[0][%d:%d,%d:%d]'%(
        #        field,band,obj_osect,xmin,xmax,ymin,ymax)
        # elif field=='s':
        #     drzimg = osectdir_s+'/%s%s_v1.9_sc03_osect%d_drz.fits[%d:%d,%d:%d]' % (field,
        #        band,obj_osect,xmin,xmax,ymin,ymax)
        # else:
        #     drzimg = osectdir_n+'/%s%s_v1.9_sc03_osect%d_drz.fits[%d:%d,%d:%d]' % (field,
        #        band,obj_osect,xmin,xmax,ymin,ymax)
        # drzimg_obj = 'images/obj%d_drz.fits' % objid
        # if os.path.exists(drzimg_obj):
        #     os.remove(drzimg_obj)
        # iraf.imcopy(c['drzimage'], drzimg_obj)
        # h = pyfits.open(drzimg_obj, mode='update')
        # h[0].header.update('xmin', xmin)
        # h[0].header.update('xmax', xmax)
        # h[0].header.update('ymin', ymin)
        # h[0].header.update('ymax', ymax) 
        # h.flush()
        # h.close()
        # copy sigma image
        # rmsimg_obj = 'images/obj%d_rms.fits' % objid
        # if os.path.exists(rmsimg_obj):
        #    os.system(rmsimg_obj)
        # iraf.imcopy(c['rmsimage'])
        #if (field=='s') & (band=='i') & (obj_osect==12):
        #    iraf.imcopy(osectdir+'/%s%s_v1.9_sc03_osect%d_rms.fits[0][%d:%d,%d:%d]' % (
        #      field,band,obj_osect,xmin,xmax,ymin,ymax),'images/obj%d_rms.fits' % objid)
        #else:
        #    iraf.imcopy(osectdir+'/%s%s_v1.9_sc03_osect%d_rms.fits[%d:%d,%d:%d]' % (
        #      field,band,obj_osect,xmin,xmax,ymin,ymax),'images/obj%d_rms.fits' % objid)
        # copy segmentation images... have to use z-band seg image?
        #if iraf.access('images/obj%d_seg.fits' % objid):
        #    os.system('rm images/obj%d_seg.fits' % objid)
        #iraf.imcopy('%sz_v1.9_sc03_osect%d_seg.fits[%d:%d,%d:%d]' % (field,obj_osect,
        #  xmin,xmax,ymin,ymax),'images/obj%d_seg.fits' % objid)
  
        # make mask image -- zero-out the included components in the mask image
        if pars['mask'] == True:
            if os.path.exists('obj%d_mask.fits' % objid):
                os.remove('obj%d_mask.fits' % objid)
            print "Making mask..."
            mask0 = mask0_dic[ixy]
            mask0_obj = mask0[ymin-1:ymax, xmin-1:xmax].copy()
            mask0_flat = mask0_obj.flatten()
            mask2_flat = where(in1d(mask0_flat, sidarray), 0, 
                                  mask0_flat)
            mask2 = mask2_flat.reshape(mask0_obj.shape)
            mask3 = where(mask2 > 0, 100., 0.)
            # normalize mask value to 100
            mask4 = ndimage.filters.gaussian_filter(mask3, GSIGMA, order=0)
            mask4 = where(mask4 <= 10., 0, 1).astype('int')
            hdu = pyfits.PrimaryHDU(mask4)
            hdu.writeto('obj%d_mask.fits' % objid)
            # make a cutout of the mask, because GALFIT won't do it
            # iraf.imcopy('%s[%d:%d,%d:%d]' % (maskname, xmin, xmax, ymin, ymax),
            #             'obj%d_mask.fits' % objid)

        # step 4: write GALFIT input file
    
        parfile = "obj%d" % objid
        f = open(parfile,'wb')
        # Write image parameters
        print >>f, "==============================================================================="
        print >>f, "# IMAGE PARAMETERS"
        print >>f, "A) %s   # Input Data image (FITS file) " % input_drz 
        print >>f, "B) galfit_output/%s   # Output data image block " % outname
        print >>f, 'C) %s    # Noise image name (made from data if blank or "none")' % input_rms
        print >>f, "D) %s         # Input PSF image for convolution (FITS file) " % pars['psfimage']
        print >>f, "E) 1                   # PSF oversampling factor relative to data "
        print >>f, "F) obj%d_mask.fits  # Bad pixel mask (FITS image or ASCi coord list; copied to obj%d_mask.fits) " % (objid, objid)
        print >>f, "G) %s           # Parameter constraint file " % pars['consfile']
        print >>f, "H) %d %d %d %d  # Image region to fit " % (xmin, xmax, ymin, ymax)
        print >>f, "I) %d  %d       # Size of the convolution box (x y)" % (cbox,cbox)
        print >>f, "J) %f             # Magnitude photometric zeropoint " % pars['magzpt']
        print >>f, "K) %.2f %.2f      # Plate scale (dx dy). " % (pars['pixscale'], pars['pixscale'])
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

        for jj in range(len(sidarray)):
            # component(f,j+1,xposarray[j]-xmin+1,yposarray[j]-ymin+1,magarray[j],
            #   rearray[j],axrat_array[j],theta_array[j])
            component(f, jj+1, xposarray[jj], yposarray[jj], magarray[jj],
                      rearray[jj], axrat_array[jj], theta_array[jj])
    
        # Make sky component
        #bkgnd = skyest('images/obj%d_drz.fits'%objid,'images/obj%d_seg.fits'%objid,5.)
        bkgnd = 0.
        sky(f,len(sidarray)+1,bkgnd)
        ncomp = len(sidarray)
        f.close()

        # step 5: run GALFIT
        t1 = time.time()
        galfitrun = subprocess.Popen([os.path.join(pars['galfitpath'],'galfit'), parfile],
                                     stdout=sys.stdout)
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
            print "Killed GALFIT, GALFIT froze"

        retcode = galfitrun.returncode
        t2 = time.time()
        try:
            print >>OUT3, "returned code is %d" % retcode
        except:
            print >>OUT3, "returned code is ", retcode
        dt = t2 - t1
        if retcode != 7:    # not sure about the behavior when GALFIT crashes...!!!
            print >>OUT3, "%d  %s" % (objid,parfile)
        else:
            print >>OUT4, "%s %.1f seconds" % (outname,dt)
        # if (dropband=='b') & (field=='n'):
        #     os.system('mv obj%d objfiles/bdrops/north' % objid)
        #     os.system('mv images/*obj%d*.fits images/bdrops/north' % objid)
        # elif (dropband=='b') & (field=='s'):
        #     os.system('mv obj%d objfiles/bdrops/south' % objid)
        #     os.system('mv images/*obj%d*.fits images/bdrops/south' % objid)
        # elif (dropband=='v') & (field=='n'):
        #     os.system('mv obj%d objfiles/vdrops/north' % objid)
        #     os.system('mv images/*obj%d*.fits images/vdrops/north' % objid)
        # else: 
        #     os.system('mv obj%d objfiles/vdrops/south' % objid)
        #     os.system('mv images/*obj%d*.fits images/vdrops/south' % objid)
        # Update the output image header
        if os.path.exists('galfit_output/'+outname):
            hout = pyfits.open('galfit_output/'+outname, mode='update')
            hout[2].header['objnum'] = objid
            hout[2].header['datain'] = pars['drzimage']
            hout[2].header['noise'] = pars['rmsimage']
            hout[2].header['xmin'] = xmin
            hout[2].header['xmax'] = xmax
            hout[2].header['ymin'] = ymin
            hout[2].header['ymax'] = ymax
            hout[2].header['ncomp'] = ncomp
            hout.flush()
            hout.close()
        # Move files to their places
        # os.system('mv %s galfit_output/' % outname)
        os.system('mv obj%d galfit_input/' % objid) 
        os.system('mv obj%d_mask.fits galfit_input/' % objid)

    OUT3.close()
    OUT4.close()


    # Finishes GALFIT fitting
    print "Finishes GALFIT fitting"

if __name__ == "__main__":
    parfile = sys.argv[1]
    run_galfit(parfile)