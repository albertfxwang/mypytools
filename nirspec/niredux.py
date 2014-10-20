import pyfits,numpy
from scipy import ndimage,signal,interpolate,optimize
from mostools import specfit as spf
from mostools import spectools as st
from mostools import ycorrect
import indexTricks as iT
import special_functions as sf
from nirspec_filters import nirspec_filters
from nirspec.nirspecTools import *
import nirspec,crsub
import matplotlib.pyplot as plt


def niredux(flats,stars,imgpairs,oname,isstar=False,
            docrrej=True,userlow=0.,userhigh=0.):

    #
    # Set detector characteristics. Numbers from NIRSPEC specs web site:  
    #  http://www2.keck.hawaii.edu/inst/nirspec/Specifications.html
    #
    # First, set plt.interactive to False, otherwise it won't run.
    plt.interactive(False)
    RN = 25.0  # Read noise in electrons
    gain = 4.0 # Gain in e-/ADU
    
    #
    # Get meta-data
    #
    try:
        tmp = nirspecOpen(imgpairs[0][0])[0]
    except:
        print ""
        print "ERROR. Cannot open first data file."
        print ""
        return
    ny,nx = tmp.data.shape
    tmphdr = tmp.header.copy()
    wstart,wend,disp,smoothKernel = nirspec_filters(tmp.header['FILNAME'])
    wmodel = {'blue':wstart,'red':wend,'scale':disp}
    if userlow>0.:
        wmodel['blue'] = userlow
    if userhigh>0.:
        wmodel['red'] = userhigh
    del tmp
    if flats is None:
        doflat = False
    else:
        doflat = True

    #
    # Make skymodel
    #
    skymodel_file = nirspec.__path__[0]+'/nirspec_skymodel.dat'
    skymodel = numpy.load(skymodel_file)

    # This might need to change depending on the filter, particularly
    #   the amount of smoothing
    wave = numpy.arange(wstart-150.,wend+150.,0.2)
    skymod = interpolate.splev(wave,skymodel)
    skymod = ndimage.gaussian_filter(skymod,smoothKernel)
    model = interpolate.splrep(wave,skymod)
    wmodel['model'] = model

    wave = numpy.arange(wstart,wend,disp)
    mkskymod = interpolate.splev(wave,model)

    #
    # Make flat by median-stacking and use the flat to do edge detection
    #
    if doflat:
        print ""
        print "Making initial flat"
        print "-------------------"
        nflat = len(flats)
        flat = numpy.empty((nflat,ny,nx))
        count = 0
        for img in flats:
            flat[count] = nirspecOpen(img,verbose=True)[0].data.copy()
            count += 1
        flat = numpy.median(flat,0)

        #
        # Get good edges
        #
        low,high,left,right = illumination(flat)

    #
    # If no flats, then use pre-determined edges
    #
    else:
        print ""
        print "** Not using flat-field for this data reduction **"
        print ""
        low = 267
        high = 500
        left = 76
        right = 1148

    #
    # Obtain the y-distortion solution
    #

    star = makeTraceIm(stars,nx,ny)
    ysoln = getYsoln(star,len(stars),low,high)
    yforw,yback = ysoln[0],ysoln[1]

    #
    # Make the normalized flat
    # DOES NOT ACCOUNT FOR TILTED WAVELENGTH SOLUTION
    #
    if doflat:
        flat = getFlatNorm(flat,ysoln,low,high)
        flat[flat<=0] = 1.
        pyfits.PrimaryHDU(flat).writeto('flat.fits',clobber=True)

    #
    # Define the output grid
    #
    scale = wmodel['scale']
    owgrid = numpy.arange(wmodel['blue'],wmodel['red']+0.1,scale/2.)
    if owgrid.size%2==1:
        owgrid = owgrid[:-1]

    #
    # Reduce the science images
    #
    print ""
    print "Processing science images, in pairs"
    print "-----------------------------------"
    wsolution = None
    wmodel['scale'] = scale/2.
    count = 0
    done = []
    for files in imgpairs:
        if len(files)==2:
            pos,neg = files
        else:
            pos = neg = files
            isstar = True
        pimg = nirspecOpen(pos,verbose=True)[0].data.astype(numpy.float32)
        nimg = nirspecOpen(neg,verbose=True)[0].data.astype(numpy.float32)
        pimg[pimg<=0.] = 1.
        nimg[nimg<=0.] = 1.
        time = nirspecOpen(neg)[0].header['ITIME']

        # Two-pass `cosmic ray' detection (and hot pixels)
        if docrrej:
            print "Doing cosmic ray rejection"
            for i in range(2):
                if isstar:
                    continue
                pcleanimg = crsub.crsub(pimg[:,500:]-nimg[:,500:],nimg[:,500:],
                                        niter=1)+nimg[:,500:]
                ncleanimg = crsub.crsub(nimg[:,500:]-pimg[:,500:],pimg[:,500:],
                                        niter=1)+pimg[:,500:]
                pimg[:,500:] = pcleanimg.copy()
                nimg[:,500:] = ncleanimg.copy()

        #
        # Find line-straigtening solution
        #
        print "Finding x-distortion so sky lines can be straightened"
        sky = getStraight(pimg,ysoln,low,False)
        if wsolution is None:
            xsoln = spf.new_xsoln(sky,xord=3,yord=3,tol=2.)
        else:
            xsoln = spf.update_xsoln(sky,xsoln,xord=3,yord=3)

        # Apply x-distortion correction (ie straighten the lines)
        skymod = st.resample1d(sky,xsoln['forw'],'x')

        # Sky is the median along the slit
        sky = numpy.median(skymod,0)
        sky = sky[left*2:right*2].copy()
        sky /= sky.mean()

        #
        # Find wavelength solution
        #
        print "Finding wavelength solution"
        print "wmodel", wmodel
        wsolution,rwsoln = getWsoln(sky,numpy.arange(sky.size)+left*2,
                                    wsolution,wmodel)
        print wsolution
        print rwsoln

        # Create difference image
        img = pimg-nimg
        if doflat:
            img /= flat

        # Resample to output grid
        print "Resampling to output wavelength grid"
        img = resampleData(img,owgrid,low,high,ysoln,xsoln,rwsoln)
        pimg = resampleData(pimg,owgrid,low,high,ysoln,xsoln,rwsoln)
        nimg = resampleData(nimg,owgrid,low,high,ysoln,xsoln,rwsoln)

        # Subtraction isn't perfect, so do one more
        bg = numpy.median(img,0)
        img -= bg

        # Create a variance image.  For now use var(pix_i) = N_e(pix_i) + RN**2
        # Since the output is the difference of two input images (pimg - nimg)
        #  the variance of the output image is 
        #     var(img) = var(pimg) + var(nimg)
        #              = N_e(pix_i,pimg) + RN**2 + N_e(pix_i,nimg) + RN**2
        # Alternative method is commented out below
        # s = numpy.sort(img,0)
        # s = s[15:-15].std(0)
        # s = s.repeat(img.shape[0]).reshape(img.shape[::-1]).T
        # var = s**2
        img *= gain
        var = gain * (pimg + nimg) + 2.0*RN**2
        vhdu = pyfits.ImageHDU(var)
        vhdu.header.update('CTYPE1','LINEAR')
        vhdu.header.update('CRVAL1',wmodel['blue'])
        vhdu.header.update('CRPIX1',1.)
        vhdu.header.update('CD1_1',scale)
        vhdu.header.update('CRVAL2',1.)
        vhdu.header.update('CRPIX2',1.)
        vhdu.header.update('CD2_2',1.)
        vhdu.header.update('CD2_1',0.)
        vhdu.header.update('CD1_2',0.)

        if pos not in done:
            hdu = pyfits.PrimaryHDU(img)
            hdu.header.update('FILENAME',pos)
            hdu.header.update('CTYPE1','LINEAR')
            hdu.header.update('CRVAL1',wmodel['blue'])
            hdu.header.update('CRPIX1',1.)
            hdu.header.update('CD1_1',scale)
            hdu.header.update('CRVAL2',1.)
            hdu.header.update('CRPIX2',1.)
            hdu.header.update('CD2_2',1.)
            hdu.header.update('CD2_1',0.)
            hdu.header.update('CD1_2',0.)
            hdu = pyfits.HDUList([hdu,vhdu])#,whdu,shdu])
            hdu.writeto('%s_%02d.fits'%(oname,count),clobber=True)
            print "Writing output file %s_%02d.fits" % (oname,count)
            done.append(pos)
            count += 1
        if neg not in done:
            hdu = pyfits.PrimaryHDU(img*-1)
            hdu.header.update('FILENAME',neg)
            hdu.header.update('CTYPE1','LINEAR')
            hdu.header.update('CRVAL1',wmodel['blue'])
            hdu.header.update('CRPIX1',1.)
            hdu.header.update('CD1_1',scale)
            hdu.header.update('CRVAL2',1.)
            hdu.header.update('CRPIX2',1.)
            hdu.header.update('CD2_2',1.)
            hdu.header.update('CD2_1',0.)
            hdu.header.update('CD1_2',0.)
            hdu = pyfits.HDUList([hdu,vhdu])#,whdu,shdu])
            hdu.writeto('%s_%02d.fits'%(oname,count),clobber=True)
            print "Writing output file %s_%02d.fits" % (oname,count)
            done.append(neg)
            count += 1
        print ""
    # Now get back plt.interactive
    plt.interactive(True)

def extract(filename,outname,owave=None):
    d = nirspecOpen(filename)
    w = d[1].data.copy()
    s = d[2].data.copy()
    d = d[0].data[:-20].copy()

    if owave is None:
        w0 = w[0]
        w1 = w[-1]
        owave = numpy.linspace(w0,w1,w.size)

    trace = d.sum(1)
    maxd = ndimage.maximum_filter(trace,3)
    peaks = numpy.where(maxd==trace)[0]
    peak = peaks[maxd[peaks].argmax()-1]

    spec = d[peak-3:peak+4].sum(0)

    specmod = interpolate.splrep(w,spec)
    outspec = interpolate.splev(owave,specmod)

    skymod = interpolate.splrep(w,s)
    sky = interpolate.splev(owave,skymod)

    specmean = outspec.mean()
    outspec /= specmean

    var = ((sky*7+outspec)*5. + 2*25**2)/(5.*specmean)**2

    st.make_spec(outspec,var,owave,outname,clobber=True)
